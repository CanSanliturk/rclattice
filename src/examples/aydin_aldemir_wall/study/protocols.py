"""Named cyclic protocols (PLAN.md §6).

**THE PROTOCOL IS INVENTED.** The 2019 paper never prints the loading history — not the number of
cycles, the amplitudes, or the order — and the 2017 test paper is not in the repo. Fig. 10(b)
constrains it without giving it back: the frame CLIPS at +16 mm while the test reached ~20 mm, dot
density shows no clean tip spikes, and 2-4 distinct dot bands between 6 and 14 mm suggest on the
order of 4-6 full loops rather than 16. So ONE CYCLE PER LEVEL is the defensible shape, and every
report says the protocol is an assumption rather than a recovery.

Amplitudes are DRIFTS, converted to millimetres against the panel's own height, because the drift
denominator moves with the panel (PLAN §2).

Costs below are from measured step rates: ~0.032 h per mm of total drive path without bond, ~0.110
with. The path is dominated by the largest amplitude, so a level ladder is priced by its tail.
"""
from __future__ import annotations

# name -> drift levels (one full reversed cycle each, push first)
PRESETS: dict[str, tuple[float, ...]] = {
    "none": (),
    "probe": (0.0005, 0.0010, 0.0015),
    "to0p3": (0.0005, 0.0010, 0.0015, 0.0020, 0.0030),
    "to0p5": (0.0005, 0.0010, 0.0020, 0.0035, 0.0050),
    "to1p0": (0.0005, 0.0010, 0.0015, 0.0020, 0.0030, 0.0050, 0.0075, 0.0100),
    # EXTENSIONS of `to1p0`: the same eight levels, then deeper ones. They exist so a deeper run
    # REPRODUCES the completed 1.0% ladder before going past it — anything that differs is then the
    # new levels and not a different loading history, which a rescaled ladder could not tell apart.
    # Past 1.0% there is no test record to compare against (Fig. 10(b)'s frame clips at +16 mm and
    # the specimen reached ~20 mm), so these levels characterise the model rather than validate it.
    "ext1p5": (0.0005, 0.0010, 0.0015, 0.0020, 0.0030, 0.0050, 0.0075, 0.0100, 0.0125, 0.0150),
    "ext2p0": (0.0005, 0.0010, 0.0015, 0.0020, 0.0030, 0.0050, 0.0075, 0.0100, 0.0125, 0.0150,
               0.0175, 0.0200),
}

# The SHAPE of the eight-level ladder above, as fractions of its own target. `--proto ladder` scales
# these to whatever `--drift` asks for, so the amplitude is an ARGUMENT and not a preset: at
# --drift 0.01 it reproduces `to1p0` exactly, and any other target keeps the same spacing.
LADDER = (0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.00)


def levels(spec: str, *, target_drift: float | None = None) -> tuple[float, ...]:
    """Resolve `--proto` to drift amplitudes. THREE FORMS, in the plan's own order (§6, §7):

      * a NAMED PRESET      `--proto to0p5`
      * a SCALED LADDER     `--proto ladder --drift 0.012`  — the eight-level shape, any target
      * AN EXPLICIT LIST    `--proto 0.0005,0.001,0.002`    — "any level list" (§7)

    The scaled ladder is what makes the cyclic target a parameter rather than a fixed 4 mm.
    """
    spec = str(spec).strip()
    if spec in PRESETS:
        return PRESETS[spec]
    if spec == "ladder":
        if not target_drift:
            raise SystemExit("--proto ladder needs a target: pass --drift (e.g. --drift 0.01)")
        return tuple(f * float(target_drift) for f in LADDER)
    if "," in spec or spec.replace(".", "").replace("e", "").replace("-", "").isdigit():
        try:
            got = tuple(float(v) for v in spec.split(",") if v.strip())
        except ValueError:
            raise SystemExit(f"could not read {spec!r} as a list of drift amplitudes")
        if not got:
            raise SystemExit(f"empty level list {spec!r}")
        return got
    raise SystemExit(f"unknown protocol {spec!r}; use a preset {sorted(PRESETS)}, "
                     f"'ladder' with --drift, or a comma-separated list of drifts")


def history(name: str, *, height: float, target_drift: float | None = None,
            cycles_per_level: int = 1) -> list[float]:
    """The drive history in mm: 0 -> +a -> -a -> 0 per cycle, per level, in order.

    Amplitudes are DRIFTS scaled by this panel's own height, because the drift denominator moves
    with the panel (PLAN §2).
    """
    out: list[float] = []
    for lvl in levels(name, target_drift=target_drift):
        a = lvl * height
        for _ in range(max(1, int(cycles_per_level))):
            out += [+a, -a, 0.0]
    return out


def drive_path_mm(name: str, *, height: float, target_drift: float | None = None,
                  cycles_per_level: int = 1) -> float:
    """Total distance the actuator travels — what the run time is proportional to."""
    hist = history(name, height=height, target_drift=target_drift,
                   cycles_per_level=cycles_per_level)
    cur, total = 0.0, 0.0
    for goal in hist:
        total += abs(goal - cur)
        cur = goal
    return total


# Hours per mm of drive path. NOBOND is MEASURED (111,017 steps / 750 s over a 6.75 mm path at
# mesh 50, 2026-09-05). BONDED is provisional and revised UP from the plan's original 0.110: the
# step COUNT is exact (248,088 for the same path — bond links sit on light steel nodes, which drops
# dt_crit from 10.0 to 4.5 us), and the per-step cost is scaled by the element ratio 34,775/13,531
# until a bonded run has actually finished. 5.7x overall, not the 3.4x the plan assumed.
# nobond REVISED 2026-09-06 from the finished 1.0% ladder itself: 4,514,700 steps in 7.8 h over a
# 274.5 mm path = 0.0284 h/mm, against 0.032 measured on a monotonic push. A cyclic run is slightly
# cheaper per mm than a push because its step is sized once for the whole history.
# bond REVISED 2026-09-08 from the running 1.5% ladder, the first bonded CYCLIC run: 15,133,338
# steps at a measured 45.0 steps/s alone on the machine = 93.4 h over a 412 mm path = 0.227 h/mm.
# The old 0.158 was inferred from bonded PUSHOVERS and is 31% low, so the 1.5% ladder was quoted at
# 65 h and will take ~93. Keyed by the CURRENT bond names, with the pre-D94 spellings accepted.
H_PER_MM = {"perfect": 0.0284, "bond": 0.227,
            "nobond": 0.0284}      # pre-D94 spelling of "perfect"


def cost_hours(name: str, *, height: float, bond: bool, target_drift: float | None = None,
               cycles_per_level: int = 1) -> float:
    """Estimated wall-clock hours for a protocol — see `H_PER_MM` for what is measured."""
    return drive_path_mm(name, height=height, target_drift=target_drift,
                         cycles_per_level=cycles_per_level) * H_PER_MM["bond" if bond else "perfect"]
