"""Stage 1 pushover machinery: gravity LoadControl + DisplacementControl pushover (D18/D19).

Elastic, so the pushover must be a straight line (constant lateral stiffness) reaching the
target drift, with base shear summed from the support reactions.
"""

from rclattice import viz
from rclattice.builders import build_lattice, select_nodes
from rclattice.model import Load
from rclattice.opensees import run_gravity, run_pushover
from rclattice.problem import BoxLoad, BoxSupport, ConcreteGrade, Problem, portal_frame

# small frame (kip-in style), coarse mesh for speed
H, SPAN, COL, BEAM, THK = 24.0, 48.0, 12.0, 12.0, 6.0
MESH, EPS = 6.0, 1e-6
DU, TARGET = 0.1, 1.0


def _frame() -> Problem:
    domain = portal_frame(height=H, span=SPAN, col_depth=COL, beam_depth=BEAM, thickness=THK)
    half = COL / 2.0
    return Problem(
        ndm=2, ndf=2, domain=domain,
        material=ConcreteGrade(name="RC", E=4030.0, nu=0.2, rho=2.25e-7, fc=6.0),
        supports=[BoxSupport(box=(-half - 1.0, SPAN + half + 1.0, -EPS, EPS), fix=(1, 1))],
        loads=[
            BoxLoad(box=(-half, half, H - EPS, H + EPS), total=(0.0, -50.0)),
            BoxLoad(box=(SPAN - half, SPAN + half, H - EPS, H + EPS), total=(0.0, -50.0)),
        ],
    )


def _lateral(model) -> list[Load]:
    ids = select_nodes(model, (-COL, SPAN + COL, H - EPS, H + BEAM + EPS))
    return [Load(nid, (10.0 / len(ids), 0.0)) for nid in ids]


def test_run_gravity_converges():
    model, _ = build_lattice(_frame(), MESH, strut_area=10.0)
    res = run_gravity(model)
    assert res["ok"] == 0
    # vertical gravity -> column tops settle (uy < 0); the portal splays outward antisymmetrically
    # (left top moves -x, right top +x), so drift is small relative to the settlement.
    left = select_nodes(model, (-EPS, EPS, H - EPS, H + EPS))[0]
    right = select_nodes(model, (SPAN - EPS, SPAN + EPS, H - EPS, H + EPS))[0]
    assert res["disps"][left][1] < 0.0 and res["disps"][right][1] < 0.0
    assert res["disps"][left][0] == -res["disps"][right][0] or \
        abs(res["disps"][left][0] + res["disps"][right][0]) < 1e-9  # antisymmetric splay


def test_elastic_pushover_is_a_straight_line_to_target():
    model, _ = build_lattice(_frame(), MESH, strut_area=10.0)
    control = select_nodes(model, (-EPS, EPS, H - EPS, H + EPS))[0]
    base = select_nodes(model, (-COL, SPAN + COL, -EPS, EPS))
    res = run_pushover(model, lateral_loads=_lateral(model), control_node=control,
                       control_dof=1, dU=DU, target=TARGET, base_nodes=base)

    assert res["converged"] and res["ok"] == 0
    disp, shear = res["disp"], res["shear"]
    assert disp[-1] >= TARGET - 1e-6                 # reached the target drift
    assert all(b > a for a, b in zip(disp, disp[1:]))  # monotonic push
    # elastic => constant *incremental* (tangent) stiffness: a straight line
    ks = [(s2 - s1) / (d2 - d1) for d1, d2, s1, s2 in zip(disp, disp[1:], shear, shear[1:])]
    assert max(ks) - min(ks) < 1e-6 * max(ks)
    assert min(ks) > 0.0                             # positive lateral stiffness


def test_base_nodes_default_to_supports():
    model, _ = build_lattice(_frame(), MESH, strut_area=10.0)
    control = select_nodes(model, (-EPS, EPS, H - EPS, H + EPS))[0]
    res = run_pushover(model, lateral_loads=_lateral(model), control_node=control,
                       control_dof=1, dU=DU, target=TARGET)
    assert set(res["base_nodes"]) == {s.node for s in model.supports}


def test_pushover_curve_renders(tmp_path):
    model, _ = build_lattice(_frame(), MESH, strut_area=10.0)
    control = select_nodes(model, (-EPS, EPS, H - EPS, H + EPS))[0]
    res = run_pushover(model, lateral_loads=_lateral(model), control_node=control,
                       control_dof=1, dU=DU, target=TARGET)
    out = tmp_path / "pushover.png"
    viz.figure_pushover([{"disp": res["disp"], "shear": res["shear"], "label": "lattice"}],
                        savepath=str(out))
    assert out.exists() and out.stat().st_size > 0
