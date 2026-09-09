"""The single-strut cyclic-carrier check (PLAN.md §3, Stage 0 item 5).

**The question.** Both concrete laws in this study are `ElasticMultiLinear` and therefore
path-INDEPENDENT: a cracked strut recovers full stiffness on reload, so neither can be cycled as it
stands. A cyclic run needs a hysteretic carrier with the same envelope, and D79 measured that
`HystereticSM` holds **4.38 MPa at zero strain** after cracking — the wall does not pinch — while
`Concrete02` holds 0.29 and does. The pinching parameters are inert against this: they shape the
APPROACH to zero, not the value AT zero.

**What D79 left open**, and what this closes: Aydin's unloading is ORIGIN-ORIENTED (Fig. 1(c)), and
`HystereticSM`'s `beta` is exactly that parameter — beta = 1 is the secant through the origin. Both
of D79's runs used beta = 0. If beta = 1 brings the residual down, his envelope can be cycled; if
not, Stages 5-6 run on `Concrete02` and report that they are not running his envelope.

Two things are measured per candidate: the stress at zero strain after a tensile excursion that
cracks the strut, and whether the unloading path points at the origin.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
os.environ.setdefault("ALDEMIR_TW", "210.0")

import specimen                                                        # noqa: E402
from rclattice.materials import (concrete_lattice_aydin,               # noqa: E402
                                 concrete_lattice_aydin_cyclic,
                                 concrete_uniaxial_regularized)
from rclattice.model import Model                                      # noqa: E402
from rclattice.opensees import run_cyclic                              # noqa: E402

L = specimen.MESH          # one orthogonal strut at the study's mesh
AREA = 6324.7              # the calibrated strut area
GRADE = specimen.GRADES["wall"]
EPS_CR = GRADE.ft / GRADE.E


def one_strut(material) -> Model:
    """A 2-node truss: node 1 pinned, node 2 free along x and driven, with a small lumped mass.

    The mass is a numerical device — the runners size their step off a period, so a zero-mass DOF
    has none. It is small enough that the response is quasi-static at this speed.
    """
    m = Model(ndm=2, ndf=2)
    m.add_node(1, (0.0, 0.0))
    m.add_node(2, (L, 0.0))
    material.id = 1
    m.uniaxial_materials.append(material)
    m.add_element(1, "corotTruss", (1, 2), (AREA, 1), kind="concrete")
    from rclattice.model import Support
    m.supports.append(Support(1, (1, 1)))
    m.supports.append(Support(2, (0, 1)))
    return m


def probe(name: str, material, *, peak_strain_ratio: float = 12.0) -> dict:
    """Cycle one strut past cracking and back, and read what it holds at zero strain."""
    # STATIC displacement control, not dynamic relaxation. A single stiff strut has one DOF and a
    # period of ~1e-4 s, so relaxation marches at dt ~ 5e-7 and the solve falls over; there is also
    # nothing here for relaxation to buy, since a lone strut has no redistribution to trace through.
    peak = peak_strain_ratio * EPS_CR * L        # mm of stretch
    hist = [+peak, 0.0, -peak, 0.0]
    model = one_strut(material)
    from rclattice.model import Load
    res = run_cyclic(model, lateral_loads=[Load(2, (1.0, 0.0))], control_node=2, control_dof=1,
                     history=hist, dU=peak / 400.0, base_nodes=[1])
    disp, force = res["disp"], res["shear"]
    # The first return to zero AFTER the tensile peak: the value at zero is what decides pinching.
    i_peak = max(range(len(disp)), key=lambda i: disp[i])
    i_zero = min((i for i in range(i_peak, len(disp))), key=lambda i: abs(disp[i]))
    stress_at_zero = force[i_zero] / AREA
    # Origin-oriented unloading: compare the secant from the peak to the origin against the actual
    # unloading path, sampled midway down.
    # Origin-oriented unloading: the unloading slope against the secant to the origin. A ratio near
    # 1.0 means the path points at the origin, Aydin's own rule. It is MEANINGLESS for a
    # path-independent law, which has no separate unloading path — it retraces its own backbone, so
    # the "slope" it reports is just the local backbone slope, negative on the softening tail.
    i_mid = (i_peak + i_zero) // 2
    secant = force[i_peak] / disp[i_peak] if disp[i_peak] else float("nan")
    actual = ((force[i_peak] - force[i_mid]) / (disp[i_peak] - disp[i_mid])
              if disp[i_peak] != disp[i_mid] else float("nan"))
    return {"name": name, "stress_at_zero_MPa": stress_at_zero,
            "peak_stress_MPa": force[i_peak] / AREA,
            "unload_over_secant": actual / secant if secant else float("nan"),
            "converged": res["converged"], "steps": len(res["disp"])}


def main() -> None:
    gf = specimen.GF
    cands = [
        ("ElasticMultiLinear (lincomp/eppcomp as they stand)",
         concrete_lattice_aydin(GRADE, 1, L, Gf=gf, a3_over_a2=5.14)),
        ("HystereticSM, beta = 0   (D79's baseline)",
         concrete_lattice_aydin_cyclic(GRADE, 1, L, Gf=gf, a3_over_a2=5.14, beta=0.0)),
        ("HystereticSM, beta = 1   (ORIGIN-ORIENTED, the open item)",
         concrete_lattice_aydin_cyclic(GRADE, 1, L, Gf=gf, a3_over_a2=5.14, beta=1.0)),
        ("Concrete02, regularized  (the repo's control law)",
         concrete_uniaxial_regularized(GRADE, 1, L, Gf=gf)),
    ]
    print(f"one strut, L = {L:g} mm, A = {AREA:,.1f} mm2, cracking at eps = {EPS_CR:.3e} "
          f"(ft = {GRADE.ft} MPa)")
    print(f"driven to 12 x eps_cr in tension, back to zero, then the same in compression\n")
    print(f"  {'carrier':<52s}{'peak':>9s}{'AT ZERO':>10s}{'unload/secant':>15s}")
    rows = []
    for name, mat in cands:
        r = probe(name, mat)
        rows.append(r)
        print(f"  {name:<52s}{r['peak_stress_MPa']:>9.3f}"
              f"{r['stress_at_zero_MPa']:>10.3f}{r['unload_over_secant']:>15.3f}")
    print("\n  'AT ZERO' is the stress the strut still carries when its strain returns to zero.")
    print("  A wall pinches only if that is near zero; D79 measured 4.38 MPa for HystereticSM at")
    print("  beta = 0 and 0.29 for Concrete02. 'unload/secant' near 1.0 means the unloading path")
    print("  points at the origin, which is Aydin's own rule (Fig. 1c).")


if __name__ == "__main__":
    main()
