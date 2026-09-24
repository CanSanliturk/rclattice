"""One RW2 status tick: regenerate the board, write the live-doc widget code and running paragraph.

    uv run python examples/thomsen_wallace_wall/board/tick.py

Prints NEW or NOOP (a 50k-step sample landed since the last tick, or not) plus a one-line summary.
Writes out/doc_cyc_widget.js (whole widget source for the live force-deformation doc's cyclic chart,
node e5a0df93-74e5; the test cloud is blob/089201af-394c in that doc) and out/doc_running.md.
"""
import json, pathlib, subprocess, sys, datetime
B = pathlib.Path(__file__).parent; OUT = B / "out"
subprocess.run([sys.executable, str(B / "make_rw2_board.py")], check=True, stdout=subprocess.DEVNULL)
T = json.loads((OUT / "rw2_tick.json").read_text()); TEST = 163.284; H = 3660.0
state_f = OUT / "tick_state.json"
prev = json.loads(state_f.read_text()) if state_f.exists() else {}
step = T.get("step", 0)
new = step != prev.get("step") or T["done"] != prev.get("done")
state_f.write_text(json.dumps({"step": step, "done": T["done"]}))

pts = ",".join("{drift:%.4f,shear:%.1f}" % (s[1], s[2]) for s in T["samples"]) or "{drift:0,shear:0}"
lvl = f"{T.get('reached', 0):.2f}% reached" if step else "no sample yet"
code = ("export default () => <claude.Visualize data-claude-component='rw2cyc' sources={{"
        "cloud:{kind:'ref', ref:'blob/089201af-394c', adapter:{kind:'csv', hasHeader:true, columns:["
        "{type:{kind:'literal',type:'number'}},{type:{kind:'literal',type:'number'}}]}}, "
        f"model:{{kind:'data', data:[{pts}]}}}}}} "
        "definition={{kind:'chart', title:'Cyclic: model loops vs test', "
        f"note:'crushing cyclic, {lvl}; model = raw console samples every 50k steps; grey = digitized test loops', "
        "source:'model', axes:{x:[{column:'drift', type:'linear', title:'drift (%)', format:{kind:'number', suffix:'%'}}], "
        "y:[{title:'base shear', min:-200, max:200, format:{kind:'number', suffix:' kN'}}]}, "
        "marks:[{type:'point', y:'shear', source:'cloud', title:'test loops (digitized)', tone:{column:'shear', cutoffs:[-1000, 1000]}}, "
        "{type:'point', y:'shear', title:'model (crushing, b 0.015, iso 0.02)'}]}} />;")
(OUT / "doc_cyc_widget.js").write_text(code)

if step:
    eta = datetime.datetime.fromisoformat(T["eta"]); p, n = T["pos"], T["neg"]
    run = (f"**{T['pct']:.1f}%** — step {step:,} of {T['n']:,}, {T['rate']:.0f} steps/s, "
           f"{T['elapsed_h']:.1f} h in, ETA **{eta:%a %d %b %H:%M}** at that rate. Largest drift so far "
           f"{T['reached']:.3f}% ({T['reached']*H/100:.1f} of 73.2 mm). Peak raw samples so far "
           f"**{p[2]:+.1f} kN** ({p[2]/TEST:.3f} × test) at {p[1]:+.3f}%"
           + (f" and **{n[2]:+.1f} kN** ({abs(n[2])/TEST:.3f} ×) at {n[1]:+.3f}%" if n[2] < 0 else "")
           + " — unsmoothed 50k-step samples; the 5 ms peak comes at the end.")
else:
    run = "Built and stepping; no 50k-step sample yet."
(OUT / "doc_running.md").write_text(run)
print(("NEW" if new else "NOOP"), f"step {step:,}", f"{T.get('pct', 0):.2f}%", "alive" if T["alive"] else "SILENT", "done" if T["done"] else "")
