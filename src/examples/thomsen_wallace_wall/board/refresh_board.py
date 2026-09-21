"""Pull every progress sample from the running runs' logs into make_rw2_board.py and regenerate.

    python board/refresh_board.py "<stamp note>"
"""
import sys, pathlib, re, glob, subprocess
SC = pathlib.Path(__file__).parent; p = SC / "make_rw2_board.py"; s = p.read_text()
root = "/Users/cansanliturk/Desktop/Master/Thesis/rclattice/src/examples/output/thomsen_wallace_wall/study/"
def samples(pat):
    log = glob.glob(root + pat)[0] + "/console.log"
    return [(int(m.group(1).replace(",", "")), float(m.group(2)), float(m.group(3)), int(m.group(4)))
            for m in re.finditer(r"step\s+([\d,]+)/[\d,]+\s+drift\s+\+([\d.]+)%\s+shear\s+([+-][\d.]+) kN\s+\[\s*(\d+)s\]", open(log).read())]
uni = samples("2026-09-20_174437*"); gra = samples("2026-09-20_174435*")
fmt = lambda lst: "[" + ",".join(f"({a},{b},{c},{d})" for a, b, c, d in lst) + "]"
s = re.sub(r'("key": "uniform".*?"samples": )\[.*?\]\}', lambda m: m.group(1) + fmt(uni) + "}", s, flags=re.S)
s = re.sub(r'("key": "graded".*?"samples": )\[.*?\]\}', lambda m: m.group(1) + fmt(gra) + "}", s, flags=re.S)
note = sys.argv[1] if len(sys.argv) > 1 else None
if note: s = re.sub(r'"stamp_note": "[^"]*"', '"stamp_note": ' + repr(note), s)
p.write_text(s); subprocess.run([sys.executable, str(p)], check=True)
print("uniform", uni[-1], "graded", gra[-1])
