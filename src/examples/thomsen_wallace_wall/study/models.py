"""Resolved parameters -> the RW2 model to analyse.

This module owns the MAPPING ONLY. Every number, grade and helper comes from the parent package
(`specimen.py`, `build.py`), so the study cannot drift away from the specimen it is about.
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import specimen                                                    # noqa: E402
from build import (calibrate, cracking_shear, describe, strut_groups,   # noqa: E402
                   strut_life, wall_lattice)

PAPER_TAIL = (specimen.AYDIN_A2, specimen.AYDIN_A3)      # Table 1: 80 / 350, fitted at 19 mm


def panel(params: dict) -> tuple[float, float]:
    return specimen.LW, specimen.HW


def drift_denominator(params: dict) -> float:
    return specimen.HW


def material_choice(params: dict):
    comp = params["comp"]
    if comp not in ("linear", "crushing", "capped"):
        raise SystemExit(f"unknown comp {comp!r}")
    material = "concrete02" if comp == "crushing" else "aydin"
    cap = float(params["fc"]) * float(params["fcx"]) if comp == "capped" else None
    tail = params["tail"]
    if tail == "paper":
        if comp == "crushing":
            raise NotImplementedError("tail='paper' means Aydin's printed a2/a3, which only exist "
                                      "on his trilinear backbone (comp linear or capped)")
        return material, PAPER_TAIL, cap
    if tail != "solved":
        raise SystemExit(f"unknown tail {tail!r}")
    return material, None, cap


def bond_kwargs(params: dict) -> dict:
    b = params["bond"]
    if b == "perfect":
        return {}
    if b.startswith("bond") and b[4:].isdigit():
        residual = int(b[4:]) / 100.0
        if not 0.0 < residual <= 1.0:
            raise SystemExit(f"bond residual out of range in {b!r}")
        return {"bond": True, "bond_residual": residual,
                "bond_damage": bool(params.get("bond_damage"))}
    raise SystemExit(f"unknown bond {b!r}; use 'perfect' or 'bondNN'")


def gf_factor(params: dict) -> float:
    return float(params["gf"]) / specimen.GF


def apply_material_overrides(params: dict) -> dict | None:
    """`fc`, `ft`, `nu` and `fy` as parameters, applied to the specimen's own grades.

    `epsc0` is RE-DERIVED as 2*fc/E (D56). Mutates `specimen.GRADES` / `specimen.STEEL` in place —
    the dict/object `build.py` reads — and returns the change so the run records it.
    """
    g = specimen.GRADES["wall"]
    fc, ft, nu = float(params["fc"]), float(params["ft"]), float(params["nu"])
    out = {}
    if fc != g.fc or ft != g.ft or nu != g.nu:
        epsc0 = 2.0 * fc / g.E
        specimen.GRADES["wall"] = dataclasses.replace(
            g, fc=fc, ft=ft, nu=nu, epsc0=epsc0, fcu=0.2 * fc, epsU=8.0 * epsc0)
        out.update({"fc": [g.fc, fc], "ft": [g.ft, ft], "nu": [g.nu, nu],
                    "epsc0": [g.epsc0, epsc0], "note": "epsc0 re-derived as 2fc/E (D56)"})
    fy = float(params["fy"])
    if fy != specimen.STEEL.fy:
        specimen.STEEL = dataclasses.replace(specimen.STEEL, name=f"Gr-fy{fy:g}", fy=fy)
        out["fy"] = [414.0, fy]
    return out or None


def steel_b(params: dict) -> float | None:
    b = float(params["steel_b"])
    return None if b == specimen.STEEL.b else b


def build(params: dict):
    """(model, calibration, meta) for one cell. No solver."""
    if params["fcx"] != 1.0 and params["comp"] != "capped":
        raise SystemExit("fcx only means anything under comp='capped'")
    length, height = panel(params)
    grid, mesh, horizon = params["grid"], float(params["mesh"]), float(params["horizon"])
    material, tail, fc_cap = material_choice(params)
    overrides = apply_material_overrides(params)
    cal = calibrate(mesh_size=mesh, horizon=horizon, field=params["field"])
    model, _edges = wall_lattice(
        cal.area, grid=grid, mesh_size=mesh, horizon=horizon,
        nonlinear=params["analysis"] != "elastic",
        gf_factor=gf_factor(params), material=material, tail_a2a3=tail, fc_cap=fc_cap,
        reinforced=bool(params["rebar"]), full_height_rebar=bool(params["rebar_top"]),
        steel_rupture=(float(params["steel_rupture"]) or None),
        concrete_residual=float(params["concrete_residual"]),
        steel_b=steel_b(params), length=length, height=height, **bond_kwargs(params))
    gff = gf_factor(params)
    meta = {
        "panel_mm": [length, height],
        "thickness_mm": specimen.TW,
        "drift_denominator_mm": height,
        "grid": grid,
        "field": cal.field,
        "nu": cal.nu,
        "material": material,
        "tail_a2a3": list(tail) if tail else None,
        "fc_cap_MPa": fc_cap,
        "steel_rupture": float(params["steel_rupture"]) or None,
        "concrete_residual": float(params["concrete_residual"]),
        "steel_b": float(params["steel_b"]),
        "gf_factor": gff,
        "material_overrides": overrides,
        "area_mm2": cal.area,
        "EA_N": cal.area * specimen.EC,
        "nodes": len(model.nodes),
        "elements": len(model.elements),
        "steel_nodes": len(model.steel_nodes),
        "describe": describe(model),
        "strut_life_orthogonal": strut_life(gff, mesh, grid=grid),
        "strut_life_diagonal": strut_life(gff, mesh, diagonal=True, grid=grid),
        "V_cr_N": cracking_shear(length=length, height=height, grid=grid, mesh_size=mesh)[0],
        "drift_cr": cracking_shear(length=length, height=height, grid=grid, mesh_size=mesh)[1],
        "axial_load_N": specimen.N_AXIAL,
    }
    return model, cal, meta


def load_path_groups(model, params: dict) -> dict | None:
    if not params["groups"]:
        return None
    mesh = float(params["mesh"])
    eg = {f"V:{k}": v for k, v in strut_groups(model, quantity="shear", mesh_size=mesh).items()}
    eg.update({f"M:{k}": v for k, v in strut_groups(model, quantity="moment", mesh_size=mesh).items()})
    return eg


def element_modulus(e) -> float:
    return specimen.STEEL.E0 if e.kind in ("longitudinal", "stirrup") else specimen.EC
