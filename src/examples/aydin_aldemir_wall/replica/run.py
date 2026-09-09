"""Run the Aydin replica and compare against his Table 4.

    uv run python examples/aydin_aldemir_wall/replica/run.py --elastic
    uv run python examples/aydin_aldemir_wall/replica/run.py --drift 0.0025 [--damping 0.5]

`--elastic` is the cheap first check: his Table 4 reports K_sim = 943.16 kN/mm, so if the replica's
elastic stiffness misses that, nothing downstream is worth running.
"""

from __future__ import annotations

import argparse, json, os, time
from datetime import datetime
from pathlib import Path

import numpy as np

from rclattice.builders import select_nodes
from rclattice.model import Load
from rclattice.opensees import (
    cyclic_protocol, nodal_displacements, run_cyclic_dynamic, run_pushover,
    run_pushover_dynamic,
)

from build import (
    calibrate, describe, explicit_steps_per_period, report_bond, report_calibration, wall_lattice,
)
from specimen import A_SHEAR, HORIZON, HW, LW, MESH, SPECIMEN, EPS, not_replicated, paper_for

OUT = Path(__file__).resolve().parent.parent.parent / "output" / "aydin_aldemir_wall" / "replica"


def nodes(model):
    top = select_nodes(model, (-1.0, LW + 1.0, HW - EPS, HW + EPS))
    ctrl = select_nodes(model, (-EPS, EPS, HW - EPS, HW + EPS))[0]
    base = select_nodes(model, (-1.0, LW + 1.0, -EPS, EPS))
    return top, ctrl, base


def hysteresis_logger(out: Path, a_shear: float, every: int = 1):
    """Accumulate and dump the load-displacement history AS THE RUN GOES.

    The full per-step history only reaches disk in `data.json` when a run RETURNS, which for a
    multi-day cyclic run means no loops are visible until the end — and none at all if it is
    stopped. This costs nothing: `(u, shear)` are already handed to the progress callback, so the
    logger just accumulates them and rewrites a small npz atomically. Field snapshots are dumped
    far more rarely by `snapshotter`, because those DO cost a pass over every node.
    """
    us: list[float] = []
    ss: list[float] = []
    ts: list[int] = []
    tmp, dst = out / "_hyst.tmp.npz", out / "hysteresis.npz"
    n = {"k": 0}

    def cb(i, _ntot, u, sh):
        n["k"] += 1
        if n["k"] % max(1, every):
            return
        ts.append(i); us.append(u); ss.append(sh)
        np.savez(tmp, step=np.array(ts), disp=np.array(us), shear=np.array(ss),
                 drift=np.array(us) / a_shear)
        os.replace(tmp, dst)

    return cb


def snapshotter(model, out: Path, area: float, ring_size: int = 10):
    """Progress callback that also DUMPS the displacement field to disk as the run goes.

    `capture=True` on the runners returns the peak and final fields, but only when the run RETURNS.
    Every interesting run in this study so far has been killed mid-way — after a collapse, after a
    hypothesis was settled, after a stall — and each left no field behind, which is exactly when a
    damage picture is most wanted. So the field is written on every progress tick, overwriting
    `snapshots.npz` atomically, and the PEAK field is kept alongside the latest.

    Cost is one `ops.nodeDisp` per node per DOF per tick — tens of milliseconds every few thousand
    steps, i.e. nothing.
    """
    ids = np.array(sorted(model.nodes), dtype=int)
    st = {"peak_abs": -1.0, "peak_u": None, "peak_meta": None}
    ring: list = []                # the last `ring_size` ticks, oldest first
    tmp, dst = out / "_snapshots.tmp.npz", out / "snapshots.npz"

    def field():
        d = nodal_displacements(model)
        return np.array([d[i] for i in ids], dtype=float)

    def cb(i, n, u, sh):
        cur = field()
        if abs(sh) > st["peak_abs"]:
            st.update(peak_abs=abs(sh), peak_u=cur, peak_meta=(i, u, sh))
        # A RING, not just the latest. The h=3.01 run showed peak at one tick and a fully collapsed
        # wall at the next, with the whole event in the gap — so "latest" recorded the aftermath and
        # could not separate a top-row tear that CAUSED the collapse from one that FOLLOWED it.
        ring.append((cur, np.array([i, u, sh], float)))
        del ring[:-ring_size]
        np.savez(tmp, ids=ids, u_latest=cur, meta_latest=np.array([i, u, sh], float),
                 u_peak=st["peak_u"], meta_peak=np.array(st["peak_meta"], float), area=area,
                 ring_u=np.array([r[0] for r in ring]),
                 ring_meta=np.array([r[1] for r in ring]))
        os.replace(tmp, dst)       # atomic: a kill mid-write cannot corrupt the good copy
        return cur

    cb.state = st
    return cb


def save_capture(out: Path, model, res: dict, area: float) -> None:
    """Persist the runner's own end-of-run capture, for runs that actually finish."""
    if "disps_final" not in res:
        return
    ids = np.array(sorted(model.nodes), dtype=int)
    pack = lambda d: np.array([d[i] for i in ids], dtype=float)   # noqa: E731
    np.savez(out / "capture.npz", ids=ids, u_final=pack(res["disps_final"]),
             u_peak=pack(res["disps_peak"]) if res.get("disps_peak") else pack(res["disps_final"]),
             area=area)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--elastic", action="store_true", help="elastic stiffness check only")
    ap.add_argument("--drift", type=float, default=0.0025)
    ap.add_argument("--mesh", type=float, default=MESH)
    ap.add_argument("--horizon", type=float, default=HORIZON,
                    help="strut connectivity radius in mesh sizes. HIS two Table 4 rows are 1.5 and "
                         "3.01. Changing it also moves the calibration (C = 0.621 -> 0.102) and the "
                         "trilinear a1 (Table 1 footnote d), both handled automatically.")
    ap.add_argument("--damping", type=float, default=0.5)
    ap.add_argument("--rate", type=float, default=7.6)
    ap.add_argument("--progress-every", type=int, default=250)
    ap.add_argument("--explicit", nargs="?", const="CentralDifference", default=None,
                    help="EXPLICIT integrator (D74). Required in practice: his tension-only "
                         "ElasticMultiLinear law is path-independent, so Newton flips branches and "
                         "the implicit solve dies at 0.006%% drift.")
    ap.add_argument("--steps-per-period", type=int, default=None,
                    help="explicit: sized from critical_time_step() when omitted")
    ap.add_argument("--cyclic", type=str, default=None, metavar="LEVELS",
                    help="REVERSED-CYCLIC run instead of a pushover. Comma-separated drift levels, "
                         "e.g. '0.0005,0.001,0.0015'. Cost scales with total PATH, not peak drift, "
                         "so it is ~4x a monotonic run to the same drift per cycle. Switches the "
                         "concrete to a HystereticSM (his ElasticMultiLinear is path-independent "
                         "and meaningless under reversal, D60) — see not_replicated.")
    ap.add_argument("--cycles-per-level", type=int, default=1)
    ap.add_argument("--cyclic-law", choices=("concrete02", "aydin"), default="concrete02",
                    help="concrete02 PINCHES (0.29 MPa left at zero strain); aydin keeps his "
                         "trilinear tail but holds 4.4 MPa there, so cracked struts never release "
                         "and the loops stay fat. Measured at the strut level.")
    ap.add_argument("--gf-factor", type=float, default=1.0,
                    help="scale Gf. Concrete02's shorter tail fails earlier than his trilinear at "
                         "the same Gf (3x in strain on VK3) — raise this if it collapses early.")
    ap.add_argument("--hyst-every", type=int, default=1,
                    help="log the load-displacement point every N progress ticks (1 = every tick). "
                         "Cheap: it is the (u, shear) the callback already receives.")
    ap.add_argument("--field-every", type=int, default=25,
                    help="dump the full displacement FIELD every N progress ticks. Expensive "
                         "(one ops.nodeDisp per node per DOF), so much rarer than the hysteresis.")
    ap.add_argument("--rebar-to-top", action="store_true",
                    help="run the LONGITUDINAL bars to y = HW. By default the bar grid tops out at "
                         "y = 2580, leaving a 100 mm unreinforced band that the drive is applied "
                         "through — and all three measured collapses tore inside it. The paper does "
                         "not state bar termination, so both are inferences.")
    ap.add_argument("--cap", type=int, default=0, metavar="ROWS",
                    help="make the top N rows ELASTIC — the stiff loading beam the specimen had. "
                         "Tests whether the collapse is a load-introduction artefact: the h=3.01 "
                         "run tore through the row below the drive (93%% of post-collapse strain "
                         "there, at 25%% strain, with the wall below fully unloaded). NOTE equalDOF "
                         "across the driven row would NOT test this — the drive already imposes the "
                         "same ops.sp on every top node (D53).")
    ap.add_argument("--bond", action="store_true",
                    help="HIS bond scheme (D72/D75) instead of shared nodes: each bar gets its own "
                         "steel nodes, tied to the concrete by a horizon ring of elastic-brittle "
                         "links. This is the largest single departure the replica otherwise carries. "
                         "Roughly doubles the model — see preflight.py --bond.")
    ap.add_argument("--bond-area-ratio", type=float, default=0.01,
                    help="link area as a fraction of the calibrated strut area A_t. 0.01 is the "
                         "ratio calibrated against perfect bond in the PARENT study at mesh 50; it "
                         "is not calibrated at this mesh. Check with --elastic --bond.")
    a = ap.parse_args()

    stamp = f"{datetime.now():%Y-%m-%d_%H%M%S}"
    out = OUT / (f"{stamp}_" + ("elastic" if a.elastic else f"pushover_d{a.drift:g}")
                 + (f"_h{a.horizon:g}" if a.horizon != HORIZON else "")
                 + (f"_cap{a.cap}" if a.cap else "")
                 + ("_rebartop" if a.rebar_to_top else "")
                 + (f"_cyclic-{a.cyclic_law}" if a.cyclic else "")
                 + (f"_bond{a.bond_area_ratio:g}" if a.bond else "")
                 + (f"_{a.explicit}" if getattr(a, "explicit", None) else ""))
    out.mkdir(parents=True, exist_ok=True)
    print(f"{SPECIMEN}\noutput: {out}\n")

    (out / "params.json").write_text(json.dumps({
        "mesh": a.mesh, "horizon": a.horizon, "drift": a.drift, "bond": a.bond,
        "bond_area_ratio": a.bond_area_ratio if a.bond else None,
        "thickness": float(os.environ.get("ALDEMIR_TW", 120.0)), "elastic": a.elastic,
        "explicit": a.explicit, "rate": a.rate, "damping": a.damping,
        "cap_rows": a.cap, "full_height_rebar": a.rebar_to_top,
        "cyclic": a.cyclic, "cycles_per_level": a.cycles_per_level,
        "cyclic_law": a.cyclic_law, "gf_factor": a.gf_factor,
    }, indent=2))

    PAPER = paper_for(a.horizon)      # his Table 4 row for THIS horizon, not the 1.5 one
    cal = calibrate(mesh_size=a.mesh, horizon=a.horizon)
    report_calibration(cal, mesh_size=a.mesh, horizon=a.horizon)

    bond_area = a.bond_area_ratio * cal.area if a.bond else None
    if a.bond:
        report_bond(cal.area, bond_area, mesh_size=a.mesh)

    t0 = time.time()
    model, _e = wall_lattice(cal.area, mesh_size=a.mesh, horizon=a.horizon,
                             nonlinear=not a.elastic, bond=a.bond, bond_area=bond_area,
                             cap_rows=a.cap, full_height_rebar=a.rebar_to_top,
                             cyclic=bool(a.cyclic), cyclic_law=a.cyclic_law,
                             gf_factor=a.gf_factor)
    print(f"\n{describe(model)}   [built in {time.time() - t0:.0f}s]")
    print(f"  his Table 2: {PAPER['nodes']:,} nodes, {PAPER['elements']:,} elements "
          f"(concrete only, horizon {PAPER['horizon']:g})")
    if a.bond:
        print("  the concrete counts are unchanged by --bond; steel nodes and bond links are extra, "
              "and his Table 2 does not count them either")
    top, ctrl, base = nodes(model)

    if a.cyclic:
        levels = [float(v) for v in a.cyclic.split(",")]
        peaks = [d * A_SHEAR for d in levels]
        hist = cyclic_protocol(peaks, cycles_per_level=a.cycles_per_level)
        path = sum(abs(b - x) for x, b in zip([0.0] + hist, hist))
        print(f"\n  CYCLIC: levels {[f'{d:.4%}' for d in levels]}, "
              f"{a.cycles_per_level} cycle(s) each -> {len(hist)} reversals")
        print(f"  total drive PATH {path:.1f} mm (a monotonic run to the top level is "
              f"{peaks[-1]:.1f} mm) — cost follows PATH, not peak drift")

        spp = a.steps_per_period
        if a.explicit and spp is None:
            spp, dtc, t1 = explicit_steps_per_period(model)
            print(f"  explicit: dt_crit = {dtc * 1e6:.3f} us, T1 = {t1 * 1e3:.3f} ms "
                  f"-> steps_per_period {spp:,} (dt = {t1 / spp * 1e6:.3f} us)", flush=True)
            print(f"  estimated {path / a.rate / (t1 / spp):,.0f} steps", flush=True)

        hyst = hysteresis_logger(out, A_SHEAR, every=a.hyst_every)
        snap = snapshotter(model, out, cal.area)
        tick = {"n": 0}

        def report(i, n, u, sh):
            hyst(i, n, u, sh)                       # cheap: every tick
            tick["n"] += 1
            if tick["n"] % max(1, a.field_every) == 0:
                snap(i, n, u, sh)                   # expensive: rarely
            print(f"    step {i:8d}/{n}  drift {u / A_SHEAR:+8.4%}  shear {sh / 1e3:+8.1f} kN"
                  f"  [{time.time() - t0:7.0f}s]", flush=True)

        res = run_cyclic_dynamic(model, control_node=ctrl, control_dof=1, history=hist,
                                 drive_nodes=top, base_nodes=base, rate=a.rate,
                                 damping_ratio=a.damping, quasi_static=True, capture=True,
                                 progress=report, progress_every=a.progress_every,
                                 integrator=(a.explicit,) if a.explicit else ("HHT", 0.7),
                                 **({"steps_per_period": spp} if spp else {}))
        save_capture(out, model, res, cal.area)
        sh = res["shear"]
        print(f"\nCYCLIC  peak +{max(sh) / 1e3:,.1f} / {min(sh) / 1e3:,.1f} kN"
              f"   (reached {max(abs(d) for d in res['disp']) / A_SHEAR:.4%} drift, "
              f"converged={res['converged']})")
        print(f"  measured   = {PAPER['F_exp_kN']:,.1f} kN  -> peak/test = "
              f"{max(sh) / 1e3 / PAPER['F_exp_kN']:.3f}")
        data = {"mode": "cyclic", "levels": levels, "cycles_per_level": a.cycles_per_level,
                "mesh": a.mesh, "horizon": a.horizon, "area": cal.area, "rate": a.rate,
                "damping": a.damping, "converged": res["converged"],
                "disp": res["disp"], "shear": res["shear"]}
    elif a.elastic:
        target = 0.0002 * A_SHEAR
        res = run_pushover(model, lateral_loads=[Load(n, (1.0 / len(top), 0.0)) for n in top],
                           control_node=ctrl, control_dof=1, dU=target / 40.0, target=target,
                           base_nodes=base, capture=True)
        save_capture(out, model, res, cal.area)
        k = res["shear"][-1] / res["disp"][-1]
        print(f"\nELASTIC  K = {k / 1e3:,.1f} kN/mm")
        print(f"  his K_sim  = {PAPER['K_sim_kN_per_mm']:,.1f} kN/mm   -> replica/his = "
              f"{k / 1e3 / PAPER['K_sim_kN_per_mm']:.3f}")
        print(f"  measured   = {PAPER['K_exp_kN_per_mm']:,.1f} kN/mm   -> replica/test = "
              f"{k / 1e3 / PAPER['K_exp_kN_per_mm']:.3f}")
        if a.bond:
            print("  READ AS A RATIO AGAINST A --elastic RUN WITHOUT --bond, not on its own: bond "
                  "can only ADD")
            print("  flexibility, so K/K_perfect > 1 is the D72 artefact and sizes it. NOTE the "
                  "links keep their")
            print("  brittle law here — if any of them cracks at this drift the comparison stops "
                  "being linear.")
        data = {"mode": "elastic", "K": k, "area": cal.area, "mesh": a.mesh,
                "horizon": a.horizon,
                "bond": a.bond, "bond_area": bond_area}
    else:
        target = a.drift * A_SHEAR
        snap = snapshotter(model, out, cal.area)

        def report(i, n, u, s):
            snap(i, n, u, s)
            print(f"    step {i:6d}/{n}  drift {u / A_SHEAR:+8.4%}  shear {s / 1e3:+8.1f} kN"
                  f"  [{time.time() - t0:7.0f}s]", flush=True)
        integ = (a.explicit,) if a.explicit else ("Newmark", 0.5, 0.25)
        spp = a.steps_per_period
        if a.explicit and spp is None:
            spp, dtc, t1 = explicit_steps_per_period(model)
            print(f"  explicit: dt_crit = {dtc * 1e6:.3f} us, T1 = {t1 * 1e3:.3f} ms "
                  f"-> steps_per_period {spp:,} (dt = {t1 / spp * 1e6:.3f} us)", flush=True)
        # Rigid loading beam: every driven node shares the control node's horizontal DOF, so the
        # top row must translate as one rather than each node being pushed independently.
        res = run_pushover_dynamic(model, control_node=ctrl, control_dof=1, target=target,
                                   drive_nodes=top, base_nodes=base, rate=a.rate,
                                   damping_ratio=a.damping, quasi_static=True,
                                   progress=report, progress_every=a.progress_every,
                                   integrator=integ, capture=True,
                                   **({"steps_per_period": spp} if spp else {}))
        save_capture(out, model, res, cal.area)
        peak = max(res["shear"]); ip = res["shear"].index(peak)
        print(f"\nPUSHOVER  peak {peak / 1e3:,.1f} kN at {res['disp'][ip] / A_SHEAR:.4%} drift"
              f"   (reached {res['disp'][-1] / A_SHEAR:.4%}, converged={res['converged']})")
        print(f"  his F_sim  = {PAPER['F_sim_kN']:,.1f} kN  -> replica/his  = {peak / 1e3 / PAPER['F_sim_kN']:.3f}")
        print(f"  measured   = {PAPER['F_exp_kN']:,.1f} kN  -> replica/test = {peak / 1e3 / PAPER['F_exp_kN']:.3f}")
        data = {"mode": "pushover", "drift": a.drift, "peak": peak, "area": cal.area,
                "mesh": a.mesh, "horizon": a.horizon, "damping": a.damping, "rate": a.rate, "integrator": list(integ),
                "bond": a.bond, "bond_area": bond_area,
                "converged": res["converged"], "disp": res["disp"], "shear": res["shear"]}

    (out / "command.txt").write_text(" ".join(__import__("sys").argv) + "\n")
    (out / "data.json").write_text(json.dumps(data))
    print("\nNOT REPLICATED:")
    for n in not_replicated(bond=a.bond, explicit=bool(a.explicit), cyclic=bool(a.cyclic)):
        print(f"  - {n}")
    print(f"\nsaved to {out}")


if __name__ == "__main__":
    main()
