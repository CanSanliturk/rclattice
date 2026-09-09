"""Plain-concrete cube in uniaxial COMPRESSION — the simplified specimen (D56).

A 200 mm plain-concrete cube as a 2D plane-stress square (200 x 200 mm face, 200 mm out-of-plane
thickness), squashed between smooth platens and traced by very slow dynamic relaxation. Deliberately
the simplest thing that still exercises Aydin's calibration in compression:

  * ONE material — Concrete02 at fc = 15.1 MPa, no elastic-compression variant;
  * standard small-displacement `Truss` struts, no corotational option;
  * rollers under the base with a single pinned node for stability;
  * strut area from Aydin's energy balance, unchanged.

Units: N, mm (stresses MPa). Backend-agnostic — nothing here imports openseespy (D8). Builders and
the calibration live in `build.py`, the entry script in `pushover.py`.
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

from rclattice.builders import select_nodes
from rclattice.model import Load
from rclattice.problem import BoxSupport, ConcreteGrade, Problem, RectangleDomain

# --- geometry ----------------------------------------------------------------------------------
L, THK = 200.0, 200.0        # cube edge and out-of-plane thickness (mm)
MESH, HORIZON = 20.0, 1.5    # grid spacing (10 divisions across the face); horizon (D9) -> 8 neighbours
EPS = 1e-6
A_FACE = L * THK             # section perpendicular to the push -> stress = N / A_FACE
OUT = Path(__file__).resolve().parent.parent / "output" / "compression_cube"
RUNS = OUT / "runs"


class _Tee:
    """Duplicate a stream into the run's `console.log`, writing through immediately.

    `vk3_wall` gets its console.log from a shell wrapper that logs to /tmp and moves the file in on
    exit — necessary there because the run directory is a timestamp nobody can redirect into ahead
    of time. The cube needs no wrapper for a 50 s run, but it has the same problem: the printed
    diagnostics ARE the result on this specimen (the modulus check, the load-path split, the
    capacity accounting), and a run whose numbers live only in a terminal scrollback cannot be
    quoted next to its own figures.
    """

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


def run_dir(kind: str, *, tag: str | None = None) -> Path:
    """A fresh TIMESTAMPED directory for one run. Nothing is ever overwritten.

    Ported from `vk3_wall` (D68). The cube's runs are short, which is exactly why a fixed output
    stem is tempting and exactly why it is wrong: every sweep — horizon, platen, solver, algorithm,
    Rmax/d — is a comparison ACROSS runs, and a fixed stem silently keeps only whichever variant ran
    last. Here each run keeps its own directory, so a sweep is a listing rather than a memory.

    Also drops a `command.txt` next to the results, because a figure whose command line is lost is
    a figure that cannot be reproduced — and on this specimen the command line IS the experiment.
    """
    name = f"{datetime.now():%Y-%m-%d_%H%M%S}_{kind}" + (f"_{tag}" if tag else "")
    d = RUNS / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "command.txt").write_text(" ".join(sys.argv) + "\n")
    link = RUNS / f"latest_{kind}"
    try:                                   # convenience pointer for the next run's baseline overlay
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(d.name)
    except OSError:
        pass                               # a filesystem without symlinks is not a reason to fail

    # Tee stdout AND stderr into the run directory. stderr matters most on the runs worth keeping
    # around: a transient that blows up leaves its traceback beside the partial results rather than
    # only in the terminal. Line-buffered, so `tail -f` follows a run in flight.
    if not isinstance(sys.stdout, _Tee):
        log = open(d / "console.log", "a", buffering=1, encoding="utf-8")
        sys.stdout, sys.stderr = _Tee(sys.stdout, log), _Tee(sys.stderr, log)
    return d


def find_run(kind: str | None = None, *, name: str | None = None, complete: bool = True) -> Path:
    """Locate a run directory: an explicit `name`, else the NEWEST one of `kind`.

    `complete=True` skips directories with no `data.json`. A run that is still in flight — or one
    that was cancelled — has already created its directory and written its `command.txt`, so a bare
    "newest of this kind" would hand the caller a run with nothing in it. Results are written only
    when the analysis finishes.
    """
    if name:
        d = RUNS / name
        if not d.is_dir():
            raise SystemExit(f"no run directory {d}")
        return d
    got = sorted((p for p in RUNS.glob("*")
                  if p.is_dir() and not p.name.startswith("latest_")
                  and (kind is None or kind in p.name)), key=lambda p: p.name)
    done = [p for p in got if (p / "data.json").exists()]
    if complete and done:
        return done[-1]
    if not got:
        raise SystemExit(f"no runs under {RUNS}"
                         + (f" matching {kind!r}" if kind else "")
                         + " — run an analysis first")
    if complete and not done:
        raise SystemExit(f"no COMPLETED {kind or 'run'} under {RUNS} — "
                         f"{len(got)} directory(ies) exist but none has a data.json "
                         f"(still running, or cancelled). Newest: {got[-1].name}")
    return got[-1]

# --- material ----------------------------------------------------------------------------------
# Low-strength concrete, fc = 15.1 MPa, with the measured modulus and tensile strength of the
# SW-NC-FF test unit (Sahinkaya et al. 2025) — the same grade the wall study uses, so the two
# examples share a material.
#
# epsc0 is DERIVED, not quoted. Concrete02's initial compressive tangent is 2*fc/epsc0 regardless of
# what the grade's `E` field says, so pinning epsc0 = 2*fc/E is what makes the material consistent
# with its own modulus. Quoting the paper's 0.0022 instead would give the struts a tangent of
# 13.7 GPa while Aydin's balance was handed 16.1 GPa — a silent 15% mismatch in the one number the
# calibration exists to reproduce.
FC = 15.1
E_C = 16100.0
NU = 0.20
EPSC0 = 2.0 * FC / E_C                      # 1.876e-3
CONCRETE = ConcreteGrade("fc15", E=E_C, nu=NU, rho=2.4e-9,
                         fc=FC, epsc0=EPSC0, fcu=4.0, epsU=0.012, ft=1.5)

# Tensile fracture energy, CEB-FIP MC90: Gf = 0.030*(fcm/10)^0.7 with fcm = fc + 8 = 23.1 MPa.
GF = 0.053
GFC_FACTOR = 250.0           # compressive fracture energy Gfc = GFC_FACTOR * Gf (D20 default)
RESIDUAL_RATIO = 0.2         # floor the crushing strength at this fraction of fc (D22)

EPS_CRACK = CONCRETE.ft / CONCRETE.E        # strut cracking strain, 9.3e-5
EPS_CRUSH = CONCRETE.epsc0                  # strut peak-compression strain

# --- loading -----------------------------------------------------------------------------------
TARGET_STRAIN = 0.01                        # ~5x epsc0: well past peak, into the residual plateau
TARGET = -TARGET_STRAIN * L                 # control displacement, mm (negative = compression)

# Very slow drive, mm/s. The response is rate-independent well above this (D54 swept 137 -> 0.8 mm/s
# with identical peak, strain at peak and load split); what a slower drive buys is a cleaner
# measurement, since the inertial and damping terms riding along in the recorded reaction shrink
# with speed. At this value they are negligible.
QUASI_STATIC_RATE = 0.5


def check_mesh_alignment(mesh_size: float) -> None:
    """Fail loudly if `mesh_size` puts no node at the face centre.

    The support layout pins ux at the bottom-edge CENTRE node, so a grid straddling x = L/2 would
    leave the block with no horizontal restraint at all — singular, and confusingly so.
    """
    half = L / 2.0
    if abs(half / mesh_size - round(half / mesh_size)) > 1e-9:
        raise ValueError(f"mesh_size={mesh_size:g} mm puts no node at the face centre x={half:g} mm "
                         f"— use a divisor of {half:g} (e.g. 20, 10, 5)")


def cube_problem(mesh_size: float = MESH) -> Problem:
    """The cube on ROLLERS, with one node pinned for stability.

    Bottom edge: uy fixed, ux free — a frictionless lower platen, so the block may spread laterally
    instead of being confined by the support. That is the low-friction bound of the van Vliet &
    van Mier platen study Aydin cites (Fig. 1.10); a bonded base would build confinement cones and
    read above fc, conflating platen friction with material strength.

    One node — the bottom centre — also has ux fixed. Rollers alone leave the cube free to slide
    sideways as a rigid body, so the stiffness matrix is singular; pinning a single node removes
    exactly that mode and adds no restraint to the deformation, and putting it on the centreline
    keeps the restraint symmetric.

    The three boxes PARTITION the bottom edge rather than overlapping: `build()` emits one `ops.fix`
    per Support with no de-duplication, so a node covered by two specs would be fixed twice.

    No loads: self-weight here is 0.005 MPa, 0.03% of fc.
    """
    check_mesh_alignment(mesh_size)
    half = L / 2.0
    domain = RectangleDomain(length=L, height=L, thickness=THK, origin=(0.0, 0.0))
    supports = [
        BoxSupport(box=(-EPS, half - EPS, -EPS, EPS), fix=(0, 1)),        # rollers, left of centre
        BoxSupport(box=(half - EPS, half + EPS, -EPS, EPS), fix=(1, 1)),  # centre node: the pin
        BoxSupport(box=(half + EPS, L + EPS, -EPS, EPS), fix=(0, 1)),     # rollers, right of centre
    ]
    return Problem(ndm=2, ndf=2, domain=domain, material=CONCRETE, supports=supports, loads=[])


def top_nodes(model) -> list[int]:
    """The loaded top face — the upper platen."""
    return select_nodes(model, (-EPS, L + EPS, L - EPS, L + EPS))


def axial_loads(model) -> list[Load]:
    """Reference load pattern for the STATIC pushover: a unit downward force spread evenly over the
    top face.

    Only its SHAPE matters — `DisplacementControl` scales the magnitude — so the total is 1 N split
    across the face. Spreading it uniformly keeps the applied traction uniform, which is what a
    smooth platen delivers; the dynamic solver reaches the same condition by prescribing the same
    downward ramp at every top node instead.
    """
    ids = top_nodes(model)
    return [Load(nid, (0.0, -1.0 / len(ids))) for nid in ids]


def control_node(model) -> int:
    """Displacement control node: top face centre (axial strain = |u| / L)."""
    return select_nodes(model, (L / 2.0 - EPS, L / 2.0 + EPS, L - EPS, L + EPS))[0]


def base_nodes(model) -> list[int]:
    """The supported bottom face — where the axial reaction is summed."""
    return select_nodes(model, (-EPS, L + EPS, -EPS, EPS))
