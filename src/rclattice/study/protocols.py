"""Cyclic protocols and their cost (Aldemir PLAN.md §6, lifted in D103).

A `ProtocolSet` is a specimen's named ladders plus the MEASURED cost per millimetre of drive path.
Amplitudes are DRIFTS, converted to millimetres against the run's own drift denominator, because
that denominator can move with the panel (Aldemir PLAN §2).

Whether the protocol is the test's own or an invention is a property of the specimen, and the
specimen says so through `invented` + `invented_note`; every report repeats it.
"""
from __future__ import annotations

from dataclasses import dataclass, field

# The SHAPE of the eight-level ladder, as fractions of its own target. `--proto ladder` scales
# these to whatever `--drift` asks for, so the amplitude is an ARGUMENT and not a preset.
LADDER = (0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.00)


@dataclass
class ProtocolSet:
    presets: dict[str, tuple[float, ...]]
    h_per_mm: dict[str, float]           # "perfect" / "bond" (+ any legacy spellings) -> hours/mm
    ladder: tuple[float, ...] = LADDER
    invented: bool = True
    invented_note: str = ""              # printed by the runner and quoted by every report

    def levels(self, spec: str, *, target_drift: float | None = None) -> tuple[float, ...]:
        """Resolve `--proto` to drift amplitudes. THREE FORMS:

          * a NAMED PRESET      `--proto to0p5`
          * a SCALED LADDER     `--proto ladder --drift 0.012`  — the eight-level shape, any target
          * AN EXPLICIT LIST    `--proto 0.0005,0.001,0.002`
        """
        spec = str(spec).strip()
        if spec in self.presets:
            return self.presets[spec]
        if spec == "ladder":
            if not target_drift:
                raise SystemExit("--proto ladder needs a target: pass --drift (e.g. --drift 0.01)")
            return tuple(f * float(target_drift) for f in self.ladder)
        if "," in spec or spec.replace(".", "").replace("e", "").replace("-", "").isdigit():
            try:
                got = tuple(float(v) for v in spec.split(",") if v.strip())
            except ValueError:
                raise SystemExit(f"could not read {spec!r} as a list of drift amplitudes")
            if not got:
                raise SystemExit(f"empty level list {spec!r}")
            return got
        raise SystemExit(f"unknown protocol {spec!r}; use a preset {sorted(self.presets)}, "
                         f"'ladder' with --drift, or a comma-separated list of drifts")

    def history(self, name: str, *, height: float, target_drift: float | None = None,
                cycles_per_level: int = 1) -> list[float]:
        """The drive history in mm: 0 -> +a -> -a -> 0 per cycle, per level, in order."""
        out: list[float] = []
        for lvl in self.levels(name, target_drift=target_drift):
            a = lvl * height
            for _ in range(max(1, int(cycles_per_level))):
                out += [+a, -a, 0.0]
        return out

    def drive_path_mm(self, name: str, *, height: float, target_drift: float | None = None,
                      cycles_per_level: int = 1) -> float:
        """Total distance the actuator travels — what the run time is proportional to."""
        hist = self.history(name, height=height, target_drift=target_drift,
                            cycles_per_level=cycles_per_level)
        cur, total = 0.0, 0.0
        for goal in hist:
            total += abs(goal - cur)
            cur = goal
        return total

    def cost_hours(self, name: str, *, height: float, bond: bool,
                   target_drift: float | None = None, cycles_per_level: int = 1) -> float:
        """Estimated wall-clock hours for a protocol — see `h_per_mm` for what is measured."""
        return self.drive_path_mm(name, height=height, target_drift=target_drift,
                                  cycles_per_level=cycles_per_level) \
            * self.h_per_mm["bond" if bond else "perfect"]
