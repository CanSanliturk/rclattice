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
lvl = f"{T.get('reached', 0):.3f}% of 2.32% max" if step else "no sample yet"
# Hand-rolled SVG (small dots) so the loops read clearly: cloud from the blob (positional
# columns), model inline as small dots joined by the loop-path line. drift x-range +/-2.5%.
code = (
    "export default () => <claude.Visualize data-claude-component='rw2cyc' sources={{"
    "cloud:{kind:'ref', ref:'blob/089201af-394c', adapter:{kind:'csv', hasHeader:true, "
    "columns:[{type:{kind:'literal',type:'number'}},{type:{kind:'literal',type:'number'}}]}}, "
    "model:{kind:'data', data:[" + pts + "]}}}>"
    "{({cloud, model, datum}) => {"
    "const G='#8A94A6', M='#2E86FF';"
    "const X=(d)=>54+(d+2.5)/5.0*348, Y=(v)=>152-v/180*104;"
    "const cx=(r)=>Object.values(r)[0], cy=(r)=>Object.values(r)[1], ky=(r)=>Object.keys(r)[1];"
    "const path=model.map((p)=>X(p.drift).toFixed(1)+','+Y(p.shear).toFixed(1)).join(' ');"
    "return <svg viewBox='0 0 420 300' style={{width:'100%'}}>"
    "<text x='230' y='13' textAnchor='middle' fontSize='9' fill='currentColor' opacity='0.9'>RW2 cyclic — base shear (kN) vs drift (%)</text>"
    "<g fontSize='7.6' fill='currentColor'>"
    "<circle cx='150' cy='26' r='1' fill={G}/><text x='156' y='29' opacity='0.85'>test loops (digitized)</text>"
    "<line x1='250' y1='26.5' x2='268' y2='26.5' stroke={M} strokeWidth='1.4'/><circle cx='259' cy='26.5' r='1.2' fill={M}/><text x='273' y='29' opacity='0.85'>model (crushing)</text>"
    "</g>"
    "{[-150,-100,-50,0,50,100,150].map((v)=><g key={'g'+v}><line x1='54' y1={Y(v)} x2='406' y2={Y(v)} stroke='currentColor' opacity={v===0?0.45:0.1} strokeWidth='0.6'/><text x='49' y={Y(v)+3} textAnchor='end' fontSize='7.5' fill='currentColor' opacity='0.7'>{v}</text></g>)}"
    "{[-2,-1,0,1,2].map((d)=><g key={'x'+d}><line x1={X(d)} y1='42' x2={X(d)} y2='262' stroke='currentColor' opacity={d===0?0.45:0.1} strokeWidth='0.6'/><text x={X(d)} y='274' textAnchor='middle' fontSize='7.5' fill='currentColor' opacity='0.7'>{d}%</text></g>)}"
    "{cloud.map((p,i)=><circle key={'c'+i} cx={X(cx(p))} cy={Y(cy(p))} r='0.5' fill={G} opacity='0.5' {...datum(p, ky(p))}/>)}"
    "<polyline points={path} fill='none' stroke={M} strokeWidth='1.2' strokeLinejoin='round'/>"
    "{model.map((p,i)=><circle key={'m'+i} cx={X(p.drift)} cy={Y(p.shear)} r='1' fill={M} {...datum(p,'shear')}/>)}"
    "<text x='230' y='289' textAnchor='middle' fontSize='7.3' fill='currentColor' opacity='0.6'>crushing cyclic, " + lvl + "; small dots = raw 50k-step samples, line = loop path; grey = digitized test loops</text>"
    "</svg>;}}</claude.Visualize>;"
)
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
