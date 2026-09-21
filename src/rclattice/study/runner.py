"""`run.py`'s body: one run = one timestamped directory (Aldemir PLAN.md §9, D70; shared in D103).

    runner.main(SPEC)            # from a study's own run.py

Every parameter comes from the study's registry, so the CLI, the run-directory name, `params.json`
and the master report never drift apart. One run = one timestamped directory containing
`command.txt`, `params.json` (EVERY parameter, defaults included, plus a schema version), a teed
`console.log`, `data.json`, `figures/` and a self-contained `report.md`.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from . import metrics
from . import report as report_mod
from .spec import StudySpec


class _Tee:
    """Duplicate a stream into the run's `console.log`, writing through immediately."""

    def __init__(self, stream, log):
        self._stream, self._log = stream, log

    def write(self, s: str) -> int:
        self._stream.write(s)
        self._log.write(s)
        return len(s)

    def flush(self) -> None:
        self._stream.flush()
        self._log.flush()

    def isatty(self) -> bool:
        return self._stream.isatty()

    def __getattr__(self, name):
        return getattr(self._stream, name)


def open_run_dir(spec: StudySpec, params: dict, root: Path) -> Path:
    """A fresh timestamped directory, with `command.txt` and a teed `console.log` (D70)."""
    name = spec.registry.run_name(params)
    d = root / name
    (d / "figures").mkdir(parents=True, exist_ok=True)
    (d / "command.txt").write_text(" ".join(sys.argv) + "\n")
    if not isinstance(sys.stdout, _Tee):
        log = open(d / "console.log", "a", buffering=1, encoding="utf-8")
        sys.stdout, sys.stderr = _Tee(sys.stdout, log), _Tee(sys.stderr, log)
    return d


def gf_multiple(model, params) -> tuple[float, str]:
    """How much energy the tension tail ACTUALLY dissipates, as a multiple of Gf (PLAN §5).

    Printed on every run rather than buried: a printed a2/a3 fitted at one grid transplanted to
    another is not neutral, since the tail's energy scales with strut length. `tail=solved` is 1.0
    by construction and this is then a self-check.
    """
    mats = {m.id: m for m in model.uniaxial_materials}
    lengths, worst = [], None
    for e in model.elements:
        if e.kind != "concrete" or len(e.nodes) != 2:
            continue
        (xa, ya), (xb, yb) = model.nodes[e.nodes[0]].coords, model.nodes[e.nodes[1]].coords
        lengths.append(((xb - xa) ** 2 + (yb - ya) ** 2) ** 0.5)
        if worst is None:
            worst = (e, lengths[-1])
    if worst is None:
        return float("nan"), "no concrete struts"
    e, L = worst
    m = mats.get(e.args[-1])
    if m is None:
        return float("nan"), "no material found for a concrete strut"
    if m.mtype != "ElasticMultiLinear":
        # Concrete02's tail is bilinear and `concrete_uniaxial_regularized` sets its softening
        # slope FROM Gf and the strut length (crack-band, D20), so the multiple is 1 by
        # construction. Stated rather than left blank: PLAN §5 asks every run to print it.
        return 1.0, "1.00 x Gf by construction (Concrete02, crack-band regularized per strut)"
    a = list(m.args)
    st = a[a.index("-strain") + 1:a.index("-stress")]
    sg = a[a.index("-stress") + 1:]
    area = 0.0
    for (e0, s0), (e1, s1) in zip(zip(st, sg), list(zip(st, sg))[1:]):
        if e0 >= 0.0:                                   # tension branch only
            area += 0.5 * (s0 + s1) * (e1 - e0)
    gf_eff = area * L                                   # crack band: energy per unit crack area
    return gf_eff / float(params["gf"]), (
        f"{gf_eff / float(params['gf']):.2f} x Gf per cracking strut at L = {L:.1f} mm")


def explicit_steps_per_period(spec: StudySpec, model, *, safety: float = 0.8) -> tuple[int, float, float]:
    """Size an EXPLICIT run: `(steps_per_period, dt_crit, T1)` — D74.

    The stable step comes from the stiffest, lightest ELEMENT (`2/w_max`), never from T1. The
    runners take their step as a fraction of a period, so T1 is fetched only to express the answer
    in the units they want. Sizing off T1 directly is the implicit habit and diverges.
    """
    from rclattice.builders import critical_time_step
    from rclattice.opensees import run_modal

    if spec.element_modulus is None:
        raise SystemExit("an explicit run needs spec.element_modulus to size its step (D74)")
    dt_crit, _w, _worst = critical_time_step(model, spec.element_modulus)
    t1 = run_modal(model, 1)["periods"][0]
    return int(t1 / (safety * dt_crit)) + 1, dt_crit, t1


def bond_artefact(spec: StudySpec, params: dict, height: float, run_pushover) -> dict | None:
    """`K/K_perfect` for a bonded run, from its OWN elastic twin (Aldemir PLAN §4).

    A full-area bond ring is a parallel load path beside an intact lattice, so a bonded model is
    STIFFER than perfect bond. Measuring it per configuration is the point; calibrating it away is
    not.
    """
    if not spec.is_bonded(params):
        return None
    sel = spec.selectors
    out = {}
    for label, bond in (("bond", params[spec.bond_axis]), ("perfect", spec.perfect_bond)):
        twin = dict(params, analysis="elastic")
        twin[spec.bond_axis] = bond
        model, _cal, _meta = spec.build(twin)
        target = spec.elastic_target_drift * height
        res = run_pushover(model, lateral_loads=sel.lateral_loads(model, height=height),
                           control_node=sel.control_node(model, height=height), control_dof=1,
                           dU=target / 40.0, target=target, base_nodes=sel.base_nodes(model))
        out[label] = res["shear"][-1] / res["disp"][-1]
    out["ratio"] = out["bond"] / out["perfect"] if out["perfect"] else float("nan")
    return out


def main(spec: StudySpec, argv=None) -> None:
    R = spec.registry
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    R.add_arguments(ap)
    a = ap.parse_args(argv)
    p = R.resolve(a)
    if spec.prepare:
        spec.prepare(p)

    from rclattice import viz
    from rclattice.opensees import run_cyclic_dynamic, run_pushover, run_pushover_dynamic

    sel = spec.selectors
    refs = spec.references
    root = spec.out_root
    out = open_run_dir(spec, p, root)
    (out / "params.json").write_text(json.dumps(p, indent=2, sort_keys=True))
    print(f"run directory: {out}")
    print("parameters: " + "  ".join(f"{k}={p[k]}" for k in R.stem))

    length_hw = spec.panel(p)
    height = spec.drift_denominator(p)
    print(f"drift denominator {height:.0f} mm — cross-panel comparisons are made in MILLIMETRES, "
          f"never in drift (PLAN §2)")

    model, cal, meta = spec.build(p)
    if spec.report_calibration:
        spec.report_calibration(cal, p)
    print(f"\n{spec.specimen}  {meta['panel_mm'][0]:.0f} x {meta['panel_mm'][1]:.0f} x "
          f"{meta['thickness_mm']:.0f}, mesh {p['mesh']:g}, horizon {p['horizon']:g}")
    print(f"  {meta['describe']}")
    print(f"  strut life: orthogonal {meta['strut_life_orthogonal']:.1f}, "
          f"diagonal {meta['strut_life_diagonal']:.1f}")
    print(f"  first cracking V_cr = {meta['V_cr_N'] / 1e3:.0f} kN at {meta['drift_cr']:.4%} drift")

    gfm, gf_note = gf_multiple(model, p)
    print(f"  tension tail dissipates {gf_note}")

    data = {"run_name": out.name, "meta": meta, "gf_multiple": gfm,
            "gf_multiple_note": gf_note, "figures": []}

    # PLAN §4: every bonded run measures its own artefact before anything else is read from it.
    art = bond_artefact(spec, p, height, run_pushover)
    if art:
        data["k_over_perfect"] = art["ratio"]
        data["k_bond"], data["k_perfect"] = art["bond"], art["perfect"]
        print(f"\n  BOND ARTEFACT (elastic twins): K_bond {art['bond'] / 1e3:,.1f} kN/mm, "
              f"K_perfect {art['perfect'] / 1e3:,.1f} kN/mm  ->  K/K_perfect = {art['ratio']:.3f}")
        # MEASURE the density, never quote it (a hardcoded mesh-50 share once appeared on a
        # mesh-25 run whose real figure was 43.0%).
        hosts = {e.nodes[1] for e in model.elements if e.kind == "bond"}
        concrete = len(model.nodes) - len(model.steel_nodes)
        share = len(hosts) / concrete if concrete else float("nan")
        print("  a full-area ring is a parallel load path, so >1 is expected and is NOT calibrated")
        print(f"  away; the bar density is the variable — {share:.1%} of this model's {concrete:,} "
              f"concrete nodes carry a steel node.")

    eg = spec.load_path_groups(model, p) if spec.load_path_groups else None
    if eg:
        print(f"  load-path probe: {len(eg)} groups across a cut at y = {0.5 * float(p['mesh']):.0f} mm")

    explicit = p["integrator"] == "explicit"
    integrator = ("CentralDifference",) if explicit else ("Newmark", 0.5, 0.25)
    spp = int(p["steps_per_period"])
    if explicit and spp == 0:
        spp, dt_crit, t1 = explicit_steps_per_period(spec, model)
        print(f"\n  EXPLICIT sizing (D74): dt_crit = {dt_crit * 1e6:.3f} us, T1 = {t1 * 1e3:.3f} ms"
              f"  ->  steps_per_period {spp:,}  (dt = {t1 / spp * 1e6:.3f} us)")
        print("  sizing off T1 instead would be 15-20x too large and would diverge.")
    elif not explicit and spp == 0:
        spp = 40
    data["steps_per_period_used"] = spp

    t0 = time.time()

    def progress(i, n, u, s_):
        print(f"    step {i:8,d}/{n:,}  drift {u / height:+8.4%}  shear {s_ / 1e3:+9.1f} kN"
              f"  [{time.time() - t0:7.0f}s]", flush=True)

    analysis = p["analysis"]
    if analysis == "elastic":
        target = spec.elastic_target_drift * height
        res = run_pushover(model, lateral_loads=sel.lateral_loads(model, height=height),
                           control_node=sel.control_node(model, height=height), control_dof=1,
                           dU=target / 40.0, target=target, base_nodes=sel.base_nodes(model))
        k = res["shear"][-1] / res["disp"][-1]
        data.update({"k_lattice": k, "disp": res["disp"], "shear": res["shear"],
                     "converged": res["converged"], "steps": len(res["disp"])})
        print(f"\nELASTIC")
        print(f"  K_lattice             {k / 1e3:9,.1f} kN/mm")
        # THE STAGE 2 GATE IS K_lattice/K_continuum, NOT K/K_measured (PLAN §8). A plane-stress
        # continuum on the SAME grid with the same rebar and BCs is what the energy balance claims
        # to reproduce; the measured stiffness is a different question, entangled with the geometry
        # (D73). Judge the calibration on the continuum ratio.
        if spec.continuum:
            con = spec.continuum(p)
            res_c = run_pushover(con, lateral_loads=sel.lateral_loads(con, height=height),
                                 control_node=sel.control_node(con, height=height), control_dof=1,
                                 dU=target / 40.0, target=target, base_nodes=sel.base_nodes(con))
            k_con = res_c["shear"][-1] / res_c["disp"][-1]
            print(f"  K_continuum           {k_con / 1e3:9,.1f} kN/mm   (same grid, same rebar, same BCs)")
            print(f"  K_lattice / K_cont.   {k / k_con:9.4f}   <- THE GATE (PLAN §8 Stage 2)")
            data.update({"k_continuum": k_con, "k_over_continuum": k / k_con})
        print(f"  K_lattice / measured  {k / 1e3 / spec.measured_K_kNmm:9.4f}   "
              f"(a different question — entangled with the geometry, D73)")
    elif analysis == "static":
        target = float(p["drift"]) * height
        res = run_pushover(model, lateral_loads=sel.lateral_loads(model, height=height),
                           control_node=sel.control_node(model, height=height), control_dof=1,
                           dU=target / int(p["steps"]), target=target,
                           base_nodes=sel.base_nodes(model),
                           algorithm=("ModifiedNewton", "-initial"), element_groups=eg,
                           capture=bool(p["capture"]))
    elif analysis == "pushover":
        target = float(p["drift"]) * height
        print(f"\n  QUASI-STATIC (dynamic relaxation): {p['rate']:g} mm/s, damping "
              f"{float(p['damping']):.0%}, {'explicit CentralDifference' if explicit else 'Newmark'}")
        res = run_pushover_dynamic(model, control_node=sel.control_node(model, height=height),
                                   control_dof=1, target=target,
                                   drive_nodes=sel.drive_nodes(model, height=height),
                                   base_nodes=sel.base_nodes(model), rate=float(p["rate"]),
                                   steps_per_period=spp, damping_ratio=float(p["damping"]),
                                   quasi_static=True, element_groups=eg,
                                   capture=bool(p["capture"]), progress=progress,
                                   progress_every=int(p["progress_every"]),
                                   integrator=integrator)
    elif analysis == "cyclic":
        target = float(p["drift"])
        cycles = int(p["cycles"])
        pro = spec.protocols
        lv = pro.levels(p["proto"], target_drift=target)
        hist = pro.history(p["proto"], height=height, target_drift=target,
                           cycles_per_level=cycles)
        if not hist:
            raise SystemExit("a cyclic run needs --proto: a preset, 'ladder' with --drift, or an "
                             "explicit list of drifts." + (" And note THE PROTOCOL IS INVENTED"
                                                          if pro.invented else ""))
        path = pro.drive_path_mm(p["proto"], height=height, target_drift=target,
                                 cycles_per_level=cycles)
        cost = pro.cost_hours(p["proto"], height=height, bond=spec.is_bonded(p),
                              target_drift=target, cycles_per_level=cycles)
        print(f"\n  CYCLIC, protocol {p['proto']!r}: {len(lv)} levels x {cycles} cycle(s), "
              f"largest {lv[-1]:.3%} drift = {lv[-1] * height:.2f} mm")
        print(f"    levels (% drift): " + ", ".join(f"{v:.3%}" for v in lv))
        print(f"    drive path {path:,.0f} mm  ->  estimated {cost:,.1f} h at the measured "
              f"step cost")
        if pro.invented_note:
            for line in pro.invented_note.splitlines():
                print(f"    {line}")
        data["protocol"] = {"spec": p["proto"], "levels": list(lv), "cycles_per_level": cycles,
                            "largest_mm": lv[-1] * height, "drive_path_mm": path,
                            "estimated_hours": cost, "invented": pro.invented}
        res = run_cyclic_dynamic(model, control_node=sel.control_node(model, height=height),
                                 control_dof=1, history=hist,
                                 drive_nodes=sel.drive_nodes(model, height=height),
                                 base_nodes=sel.base_nodes(model), rate=float(p["rate"]),
                                 steps_per_period=spp, damping_ratio=float(p["damping"]),
                                 # NOTE: `run_cyclic_dynamic` takes no `element_groups` — the
                                 # load-path probe is available on the pushover runners only.
                                 quasi_static=True, capture=bool(p["capture"]),
                                 progress=progress, progress_every=int(p["progress_every"]),
                                 integrator=integrator)
    else:
        raise SystemExit(f"unknown analysis {analysis!r}")

    elapsed = time.time() - t0
    if analysis != "elastic":
        peak = max(res["shear"])
        at_peak = res["disp"][res["shear"].index(peak)] / height
        end_drift = res["disp"][-1] / height
        k0 = res["shear"][1] / res["disp"][1] if len(res["disp"]) > 1 and res["disp"][1] else float("nan")
        # A SAMPLE MAXIMUM IS NOT RESISTANCE once the run carries failure switches (D92): the
        # release of a ruptured bar rings in a few elements at a local period well below T1, and
        # `max` lands on a crest of it. Smooth over a physical window first; both are recorded, and
        # `peak_ringing_ratio` says whether the two ever differed.
        dt_marched = res.get("dt") or (res.get("T1", float("nan")) / spp)
        data.update(metrics.response_metrics(res["shear"], res["disp"], dt=dt_marched, height=height,
                                             window_s=spec.peak_window_s or metrics.WINDOW_S))
        data.update({"peak_shear": peak, "drift_at_peak": at_peak, "end_drift": end_drift,
                     "k_initial": k0, "converged": res["converged"], "steps": len(res["disp"]),
                     "disp": res["disp"], "shear": res["shear"], "groups": res.get("groups"),
                     "rate": res.get("rate"), "T1": res.get("T1"),
                     "integrator": list(res.get("integrator", integrator)),
                     "dynamic_residual": res.get("dynamic"),
                     "disps_peak": res.get("disps_peak"), "disps_final": res.get("disps_final")})
        print(f"\ntraced to {end_drift:.4%} drift (converged={res['converged']}, "
              f"{len(res['disp']):,} steps, {elapsed:,.0f}s)")
        pk, pk_drift, smoothed = metrics.headline_peak(data)
        print(f"  peak base shear   {pk / 1e3:9,.1f} kN at {pk_drift:.4%} drift"
              + (f"  ({data['peak_smoothing_window_s'] * 1e3:g} ms moving average)" if smoothed else ""))
        print(f"  measured ({refs.f_source}){spec.measured_F_kN:9,.1f} kN  ->  model/test = "
              f"{pk / 1e3 / spec.measured_F_kN:.3f}")
        note = metrics.ringing_note(data)
        if note:
            print(f"  {note}")

        # The residual on the ASCENDING branch is what says whether the peak is real resistance;
        # a headline maximum is usually the collapse transient (D75).
        dyn = res.get("dynamic")
        if dyn:
            i_peak = res["shear"].index(peak)
            asc = sorted(abs(d) for d in dyn[:max(1, i_peak)])
            data["residual_ascending_p95"] = asc[int(0.95 * (len(asc) - 1))]
            print(f"  residual, ascending branch p95: "
                  f"{data['residual_ascending_p95'] / 1e3:,.1f} kN = "
                  f"{data['residual_ascending_p95'] / abs(pk):.1%} of peak")

        drift_pct = [u / height * 100.0 for u in res["disp"]]
        curves = [{"disp": drift_pct, "shear": [s / 1e3 for s in res["shear"]],
                   "label": f"lattice ({spec.label(p)}) — peak {peak / 1e3:.0f} kN",
                   "style": {"color": "C0", "lw": 2}},
                  {"disp": [min(drift_pct), max(drift_pct)],
                   "shear": [spec.measured_F_kN] * 2,
                   "label": f"test peak, {refs.f_source} ({spec.measured_F_kN:.0f} kN)",
                   "style": {"color": "C3", "ls": "--", "lw": 1.5}}]
        viz.figure_pushover(curves, savepath=str(out / "figures" / "backbone.png"),
                            xlabel="drift ratio (%)", ylabel="base shear (kN)",
                            title=f"{spec.title} — {analysis}, {spec.label(p)}")
        data["figures"].append("backbone.png")

        # LOAD PATH (PLAN §9/§10): the share of the overturning moment carried by the flexural
        # couple (vertical struts + vertical rebar) against truss action (diagonals), step by step.
        if res.get("groups"):
            import matplotlib.pyplot as plt
            g = res["groups"]
            fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=170)
            n = min(len(drift_pct), min(len(v) for v in g.values()))
            tot = [sum(abs(g[k][i]) for k in g if k.startswith("M:")) or float("nan")
                   for i in range(n)]
            for key in sorted(k for k in g if k.startswith("M:")):
                share = [abs(g[key][i]) / tot[i] * 100.0 for i in range(n)]
                ax.plot(drift_pct[:n], share, lw=1.8, label=key[2:])
            ax.set_xlabel("drift ratio (%)")
            ax.set_ylabel("share of overturning moment (%)")
            ax.set_title("Load path across a base cut — flexural couple vs truss action")
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(out / "figures" / "load_path.png")
            plt.close(fig)
            data["figures"].append("load_path.png")

        # COMPARISON against the digitized record and the published peak, at MATCHED DISPLACEMENT.
        rows = refs.comparison_points(res["disp"], res["shear"], horizon=float(p["horizon"]))
        data["comparison"] = rows
        ref = refs.load()
        if ref is not None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(7.6, 4.8), dpi=170)
            clip = f" (clipped at +{refs.clip_mm:g} mm)" if refs.clip_mm < float("inf") else ""
            ax.scatter(ref["cloud"][:, 0], ref["cloud"][:, 1], s=0.5, color="0.75",
                       label=f"test, {refs.figure}{clip}")
            ax.plot(ref["mono_x"], ref["mono_y"], color="C3", lw=1.4, ls="--",
                    label="test envelope (UPPER BOUND)")
            key = refs.author_key(float(p["horizon"]))
            if key and key + "_x" in ref:
                ax.plot(ref[key + "_x"], ref[key + "_y"], color="C2", lw=0.9, alpha=0.7,
                        label=f"{refs.author_label}, horizon {p['horizon']:g}")
            ax.plot(res["disp"], [sh / 1e3 for sh in res["shear"]], color="C0", lw=2.0,
                    label=f"this run ({spec.label(p)})")
            if refs.clip_mm < float("inf"):
                ax.axvline(refs.clip_mm, color="0.5", lw=0.8, ls=":")
            ax.set_xlabel("top displacement (mm)")
            ax.set_ylabel("base shear (kN)")
            ax.secondary_xaxis("top", functions=(lambda mm: mm / height * 100.0,
                                                 lambda d: d * height / 100.0)
                               ).set_xlabel("drift (%)")
            ax.set_title(f"Against the test record and {refs.author_label}'s own curve")
            ax.grid(alpha=0.3)
            ax.legend(fontsize=7.5, loc="lower right")
            fig.tight_layout()
            fig.savefig(out / "figures" / "comparison.png")
            plt.close(fig)
            data["figures"].append("comparison.png")
            if rows:
                print("\n  AT MATCHED DISPLACEMENT (the peak-to-peak ratio compares two "
                      "different states, D77):")
                print(f"    {'u (mm)':>8s}{'model':>10s}{'test env':>10s}{refs.author_label:>9s}"
                      f"{'m/test':>9s}{'m/auth':>9s}")
                for r in rows:
                    print(f"    {r['at_mm']:>8.2f}{r['model_kN']:>10.1f}"
                          f"{r.get('test_envelope_kN') or float('nan'):>10.1f}"
                          f"{r.get('aydin_kN') or float('nan'):>9.1f}"
                          f"{r.get('model_over_test') or float('nan'):>9.3f}"
                          f"{r.get('model_over_aydin') or float('nan'):>9.3f}"
                          + (f"   [past the +{refs.clip_mm:g} mm clip]" if r["clipped"] else ""))

        if p["capture"] and res.get("disps_peak") is not None and spec.damage_thresholds:
            from rclattice.viz import figure_damage
            panels = [(f"at peak shear ({peak / 1e3:.0f} kN)", model, res["disps_peak"]),
                      (f"at the last step ({end_drift:.4%} drift)", model, res["disps_final"])]
            # eps_crush is a POSITIVE magnitude — passing -EPSC0 flags ~99% of struts as crushed.
            eps_crack, eps_crush, crush_label = spec.damage_thresholds(p, meta)
            figure_damage(panels, eps_crack=eps_crack, eps_crush=eps_crush,
                          crush_label=crush_label,
                          savepath=str(out / "figures" / "damage.png"))
            data["figures"].append("damage.png")

    data["elapsed_s"] = elapsed
    data["nodes"], data["elements"] = meta["nodes"], meta["elements"]
    (out / "data.json").write_text(json.dumps(data, indent=2))
    report_mod.write(spec, out, p, data)
    print(f"\nsaved to {out}\n  report.md, data.json, params.json, console.log, "
          f"figures/ ({len(data['figures'])})")
