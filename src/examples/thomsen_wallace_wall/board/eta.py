"""Print two data-driven ETAs for the running RW2 cyclic, from its own console.log:

    tangent = pace of the LAST 50k-step interval (di/dt of the last two samples)
    secant  = average pace since the run began (current step / current elapsed)

    uv run python examples/thomsen_wallace_wall/board/eta.py
"""
import re, datetime, pathlib

STUDY = pathlib.Path(__file__).parents[2] / "output" / "thomsen_wallace_wall" / "study"
GLOB = "2026-09-24_235853_cyclic_crushing*"
SAMPLE = re.compile(r"step\s+([\d,]+)/([\d,]+)\s+drift\s+[+-][\d.]+%\s+shear\s+[+-][\d.]+ kN\s+\[\s*(\d+)s\]")

d = sorted(STUDY.glob(GLOB))[-1]
m = SAMPLE.findall((d / "console.log").read_text())
s = [(int(a.replace(",", "")), int(n.replace(",", "")), int(t)) for a, n, t in m]
i, N, t = s[-1]
rem = N - i
now = datetime.datetime.now()

tan = (s[-1][0] - s[-2][0]) / max(1, s[-1][2] - s[-2][2]) if len(s) > 1 else i / max(1, t)
sec = i / max(1, t)   # from the beginning (t=0 ~ launch)
for name, rate in (("tangent (last interval)", tan), ("secant (avg since start)", sec)):
    h = rem / rate / 3600
    eta = now + datetime.timedelta(hours=h)
    print(f"{name:26s} {rate:5.1f} st/s  ->  {h:5.1f} h left  ETA {eta:%a %d %b %H:%M}")
