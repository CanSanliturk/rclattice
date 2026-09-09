"""Resolved parameters -> the model to analyse (PLAN.md §9).

This module owns the MAPPING ONLY. Every number, grade and helper comes from the parent package
(`specimen.py`, `build.py`), so the study cannot drift away from the specimen it is about: change
the wall in one place and every cell of the matrix moves with it.

An axis value that is not implemented yet raises with the exact change it needs, rather than
silently falling back to something adjacent — a matrix whose cells quietly mean different things is
worse than a matrix with holes in it.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import specimen                                                    # noqa: E402
from build import (calibrate, cracking_shear, describe, strut_groups,   # noqa: E402
                   strut_life, wall_lattice)

# Panel options, PLAN §2. `fig10a` is the specimen; the other two are the Table 2 inversion and its
# transposition, and both are ELASTIC-ONLY here — they force mesh 10, which is 125x the cost.
PANELS = {
    "fig10a": (specimen.LW, specimen.HW),
    "table2": specimen.PAPER_GRID,
    "table2t": specimen.PAPER_GRID_T,
}

# The paper's printed tension tail (Table 1). Fitted at his 20 mm grid — see `wall_lattice`.
PAPER_TAIL = (70.0, 360.0)


def panel(params: dict) -> tuple[float, float]:
    try:
        return PANELS[params["panel"]]
    except KeyError:
        raise SystemExit(f"unknown panel {params['panel']!r}; have {sorted(PANELS)}")


def drift_denominator(params: dict) -> float:
    """The height drift is measured against — it MOVES WITH THE PANEL (PLAN §2).

    The test's 20 mm is 0.89% on `fig10a` and 0.75% on `table2`, so cross-panel comparisons are made
    in millimetres and never in drift. Every run prints this.
    """
    return panel(params)[1]


def material_choice(params: dict) -> tuple[str, tuple[float, float] | None, float | None]:
    """(`wall_lattice` material name, explicit tail, compression cap) for the comp x tail cell.

    `lincomp` is his published law — compression LINEAR AT E FOREVER, crushing meant to emerge as
    indirect tensile splitting. `eppcomp` caps it at `fc * fcx` and goes perfectly plastic.
    `c02` is the repo's own control law (Concrete02), which has a real softening compression
    branch and is the only one of the three that can be cycled as it stands (D79, PLAN §3).
    """
    comp = params["comp"]
    if comp not in ("linear", "crushing", "capped"):
        raise SystemExit(f"unknown comp {comp!r}")
    material = "concrete02" if comp == "crushing" else "aydin"
    cap = None
    if comp == "capped":
        cap = float(params["fc"]) * float(params["fcx"])
    tail = params["tail"]
    if tail == "paper":
        if comp == "crushing":
            raise NotImplementedError(
                "tail='paper' means Aydin's printed a2/a3, which only exist on his trilinear "
                "backbone: it is meaningful with comp='lincomp' or 'eppcomp', not with "
                "Concrete02's bilinear tail.")
        return material, PAPER_TAIL, cap
    if tail != "solved":
        raise SystemExit(f"unknown tail {tail!r}")
    return material, None, cap


def bond_kwargs(params: dict) -> dict:
    """PLAN §4's bond law, or none.

    `bondNN` sets the residual plateau, Aydin's `a`. BOTH AVAILABLE VALUES ARE HIS AND THEY
    DISAGREE (D96): `bond70` is the PUBLISHED one (2019 paper, from Aydin 2017, for deformed bars)
    and is the library default from 2026-09-08; `bond60` is what the author said directly to this
    project on 2026-09-05 — later than the paper, specific to this replication. Every bonded run
    before 2026-09-08 used `bond60`, so a switch is a difference to measure, not a correction. The links take
    the FULL concrete-strut area — that IS the law, not an oversight — so a bonded model comes out
    stiffer than perfect bond, and every bonded run reports `K/K_perfect` so the artefact is
    measured rather than calibrated away.
    """
    b = params["bond"]
    if b in ("perfect", "nobond"):          # "nobond" is the pre-D94 spelling of perfect bond
        return {}
    if b.startswith("bond") and b[4:].isdigit():
        residual = int(b[4:]) / 100.0       # bond60 -> 0.60 of f_t held after the brittle drop
        if not 0.0 < residual <= 1.0:
            raise SystemExit(f"bond residual out of range in {b!r} (got {residual})")
        return {"bond": True, "bond_law": "aydin", "bond_residual": residual,
                "bond_damage": bool(params.get("bond_damage"))}
    if b.startswith("bond-a") and b[6:].isdigit():      # legacy bond-a06 / bond-a07
        return {"bond": True, "bond_law": "aydin", "bond_residual": int(b[6:]) / 10.0}
    raise SystemExit(f"unknown bond {b!r}; use 'perfect' or 'bondNN' (e.g. bond60)")


def gf_factor(params: dict) -> float:
    """`gf` is given absolutely; the builder wants it relative to the specimen's own Gf."""
    return float(params["gf"]) / specimen.GF


def apply_material_overrides(params: dict) -> dict | None:
    """`fc` and `ft` as parameters (PLAN §7), applied to the specimen's own grade.

    **`epsc0` IS RE-DERIVED AS 2*fc/E, never carried over (D56).** Concrete02's compressive tangent
    is `2*fc/epsc0` whatever the grade says, so holding a published `epsc0` against a changed `fc`
    silently multiplies every strut's stiffness — the trap the compression cube was built to find.
    `fcu` and `epsU` follow the same grade's own conventions.

    Mutates `specimen.GRADES` in place, which is what `build.py` reads (same dict object). Returns
    the change so the run records it, or None when nothing moved.
    """
    import dataclasses

    g = specimen.GRADES["wall"]
    fc, ft = float(params["fc"]), float(params["ft"])
    if fc == g.fc and ft == g.ft:
        return None
    epsc0 = 2.0 * fc / g.E
    specimen.GRADES["wall"] = dataclasses.replace(
        g, fc=fc, ft=ft, epsc0=epsc0, fcu=0.2 * fc, epsU=8.0 * epsc0)
    return {"fc": [g.fc, fc], "ft": [g.ft, ft], "epsc0": [g.epsc0, epsc0],
            "note": "epsc0 re-derived as 2fc/E (D56)"}


def build(params: dict):
    """(model, calibration, meta) for one cell. No solver."""
    if params["fcx"] != 1.0 and params["comp"] != "capped":
        raise SystemExit("fcx scales the strut compressive strength and only means anything under "
                         "comp='eppcomp' (PLAN §5)")
    length, height = panel(params)
    mesh, horizon = float(params["mesh"]), float(params["horizon"])
    specimen.check_mesh_alignment(mesh, length=length, height=height)

    material, tail, fc_cap = material_choice(params)
    overrides = apply_material_overrides(params)
    cal = calibrate(mesh_size=mesh, horizon=horizon)
    model, _edges = wall_lattice(
        cal.area, mesh_size=mesh, horizon=horizon,
        nonlinear=params["analysis"] != "elastic",
        gf_factor=gf_factor(params), material=material, tail_a2a3=tail, fc_cap=fc_cap,
        reinforced=bool(params["rebar"]), full_height_rebar=bool(params["rebar_top"]),
        steel_rupture=(float(params["steel_rupture"]) or None),
        concrete_residual=float(params["concrete_residual"]),
        length=length, height=height, **bond_kwargs(params))

    meta = {
        "panel_mm": [length, height],
        "thickness_mm": specimen.TW,
        "drift_denominator_mm": height,
        "material": material,
        "tail_a2a3": list(tail) if tail else None,
        "fc_cap_MPa": fc_cap,
        "steel_rupture": float(params["steel_rupture"]) or None,
        "concrete_residual": float(params["concrete_residual"]),
        "gf_factor": gf_factor(params),
        "material_overrides": overrides,
        "area_mm2": cal.area,
        "EA_N": cal.area * specimen.EC,
        "nodes": len(model.nodes),
        "elements": len(model.elements),
        "steel_nodes": len(model.steel_nodes),
        "describe": describe(model),
        # Strut life = eps_ult/eps_cr. Below ~10 a strut fails almost immediately after cracking and
        # the lattice cannot redistribute — the reason WSH3 needed --gf-factor 2 (D67).
        "strut_life_orthogonal": strut_life(gf_factor(params), mesh),
        "strut_life_diagonal": strut_life(gf_factor(params), mesh, diagonal=True),
        "V_cr_N": cracking_shear(length=length, height=height)[0],
        "drift_cr": cracking_shear(length=length, height=height)[1],
    }
    return model, cal, meta


def load_path_groups(model, params: dict) -> dict | None:
    """The `element_groups` probe: shear and moment decomposition across a base cut."""
    if not params["groups"]:
        return None
    mesh = float(params["mesh"])
    eg = {f"V:{k}": v for k, v in strut_groups(model, quantity="shear", mesh_size=mesh).items()}
    eg.update({f"M:{k}": v for k, v in strut_groups(model, quantity="moment", mesh_size=mesh).items()})
    return eg


def element_modulus(e) -> float:
    """E for `critical_time_step`, which sees only material TAGS and so has to be told (D74).

    Bond links take the CONCRETE modulus, matching the replica's convention.
    """
    return specimen.STEEL.E0 if e.kind in ("longitudinal", "stirrup") else specimen.EC


def explicit_steps_per_period(model, *, safety: float = 0.8) -> tuple[int, float, float]:
    """Size an EXPLICIT run: `(steps_per_period, dt_crit, T1)` — D74.

    The stable step comes from the stiffest, lightest ELEMENT (`2/w_max`), never from T1. The
    runners take their step as a fraction of a period, so T1 is fetched only to express the answer
    in the units they want. Sizing off T1 directly is the implicit habit and diverges here: the
    working explicit runs on this wall needed `steps_per_period` 750, not the runner's default 40.
    """
    from rclattice.builders import critical_time_step
    from rclattice.opensees import run_modal

    dt_crit, _w, _worst = critical_time_step(model, element_modulus)
    t1 = run_modal(model, 1)["periods"][0]
    return int(t1 / (safety * dt_crit)) + 1, dt_crit, t1
