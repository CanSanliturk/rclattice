"""Score a finished run into board/out/rw2_finished.json (figures embedded) for the board.

    python board/finished.py <key> <run-dir-name>     e.g. finished.py graded 2026-09-20_174435_pushover_...
"""
import sys, json, base64, pathlib, numpy as np
B = pathlib.Path(__file__).parent; OUT = B / "out"; OUT.mkdir(exist_ok=True)
key, name = sys.argv[1], sys.argv[2]
run = B.parents[1] / "output" / "thomsen_wallace_wall" / "study" / name
d = json.load(open(run / "data.json")); H = 3660.0
b64 = lambda f: "data:image/png;base64," + base64.b64encode((run / "figures" / f).read_bytes()).decode()
g = d.get("groups") or {}; disp = np.array(d["disp"]); lp = []
for tgt in (0.3, 0.6, 1.0, 1.5, 2.0, 2.5):
    i = int(np.argmin(abs(disp / H * 100 - tgt)))
    tot = sum(abs(g[k][i]) for k in g if k.startswith("M:")) or 1
    sh = {k[2:]: abs(g[k][i]) / tot * 100 for k in g if k.startswith("M:")}
    lp.append((tgt, sh.get("diagonal", 0), sh.get("rebar", 0), sh.get("vertical", 0)))
fin = json.load(open(OUT / "rw2_finished.json")) if (OUT / "rw2_finished.json").exists() else {}
fin[key] = {"peak": d["peak_shear_smooth"] / 1e3, "peak_drift": d["drift_at_peak_smooth"] * 100,
            "plateau": [x * 100 for x in d["peak_plateau"]], "raw": d["peak_shear_raw"] / 1e3,
            "ring": d["peak_ringing_ratio"], "window_ms": d["peak_smoothing_window_s"] * 1e3,
            "end_shear": d["shear"][-1] / 1e3, "end_drift": d["end_drift"] * 100,
            "residual_pct": d["residual_ascending_p95"] / d["peak_shear_smooth"] * 100,
            "hours": d["elapsed_s"] / 3600, "steps": d["steps"], "capacity": d.get("drift_capacity"),
            "comparison": d["comparison"], "loadpath": lp,
            "fig_comparison": b64("comparison.png"), "fig_damage": b64("damage.png"), "fig_loadpath": b64("load_path.png")}
json.dump(fin, open(OUT / "rw2_finished.json", "w")); print("scored", key, round(fin[key]["peak"], 1), "kN")
