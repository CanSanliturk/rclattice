"""Generate the personal RW2 run board (an HTML page republished to ONE fixed artifact URL).

    python board/refresh_board.py "<stamp note>"     # pulls samples from the run logs, regenerates
    python board/finished.py <run-dir-name>          # scores a finished run into rw2_finished.json

Then republish `board/out/rw2_board.html` to https://claude.ai/code/artifact/168c1d28-6f3a-4caa-8da9-1674c5cf6c8e
(pass that URL to the Artifact tool; a publish without it creates a second board, which the user
rejected). `rw2_ref.json` is the digitized record (from data/fig9b.npz) baked in for the charts.
STATE holds what is running; FIN (rw2_finished.json) holds finished runs with their figures.
"""
import json, pathlib, datetime
SC = pathlib.Path(__file__).parent
OUT = SC / "out"; OUT.mkdir(exist_ok=True)
REF = json.load(open(SC / "rw2_ref.json"))
NOW = datetime.datetime.now()
T0 = datetime.datetime(2026, 9, 20, 17, 44, 35)
TEST = 163.284
STATE = {
    "stamp_note": 'uniform run flat at ~124 kN through 1.28% drift; graded run at its peak region, 137.8 kN at 0.455%',
    "runs": [
        {"key": "uniform", "name": "uniform 30.5 mm grid", "sub": "bars snapped to the grid · grid-mode control",
         "elements": "4,961 nodes · 22,002 elements", "n": 2726277, "color": "var(--uni)",
         "samples": [(100000,0.0917,75.5,1904),(200000,0.1834,100.9,4057),(300000,0.2751,116.6,6207),(400000,0.3668,130.0,8335),(500000,0.4585,138.3,10372),(600000,0.5502,137.9,12108),(700000,0.6419,128.8,13829),(800000,0.7336,129.6,15613),(900000,0.8253,127.0,17156),(1000000,0.917,125.4,18792),(1100000,1.0087,126.0,20476),(1200000,1.1004,124.4,22169),(1300000,1.1921,124.1,23827),(1400000,1.2838,124.1,25469),(1500000,1.3755,122.1,27190),(1600000,1.4672,123.7,28899),(1700000,1.5589,123.4,30612),(1800000,1.6506,123.7,32310),(1900000,1.7423,122.0,33988),(2000000,1.834,123.6,35662),(2100000,1.9257,123.4,37393),(2200000,2.0174,121.6,39117),(2300000,2.1091,117.9,40792),(2400000,2.2008,115.4,42463),(2500000,2.2925,116.9,44143),(2600000,2.3842,115.7,45842),(2700000,2.4759,116.7,47514)]},
        {"key": "graded", "name": "graded 25 mm grid", "sub": "bars on nodes exactly · the deliverable run",
         "elements": "7,446 nodes · 32,423 elements", "n": 4395483, "color": "var(--gra)",
         "samples": [(100000,0.0569,53.6,2830),(200000,0.1138,85.0,6235),(300000,0.1706,98.1,9602),(400000,0.2275,108.4,13121),(500000,0.2844,118.4,16748),(600000,0.3413,126.6,20006),(700000,0.3981,134.8,23436),(800000,0.455,137.8,26563),(900000,0.5119,138.2,29194),(1000000,0.5688,135.8,31776),(1100000,0.6256,135.5,34267),(1200000,0.6825,135.7,36780),(1300000,0.7394,135.5,39340),(1400000,0.7963,134.7,41846),(1500000,0.8531,129.3,44364),(1600000,0.91,129.1,46885),(1700000,0.9669,128.9,49130),(1800000,1.0238,128.9,51161),(1900000,1.0807,129.2,53138),(2000000,1.1375,127.4,55424),(2100000,1.1944,127.9,58510),(2200000,1.2513,127.9,61658),(2300000,1.3082,127.3,64848),(2400000,1.365,127.7,68049),(2500000,1.4219,127.0,71304),(2600000,1.4788,125.9,74545),(2700000,1.5357,126.5,77806),(2800000,1.5925,126.5,81180),(2900000,1.6494,125.7,84629),(3000000,1.7063,121.6,88064),(3100000,1.7632,121.4,91467),(3200000,1.82,121.0,94336),(3300000,1.8769,118.8,96382),(3400000,1.9338,119.9,98427),(3500000,1.9907,119.7,100476),(3600000,2.0476,120.1,102526),(3700000,2.1044,120.2,104579),(3800000,2.1613,120.1,106630),(3900000,2.2182,119.8,108681),(4000000,2.2751,119.1,110733),(4100000,2.3319,118.7,112771),(4200000,2.3888,117.3,114783),(4300000,2.4457,117.2,117120)]},
    ],
}
el = (NOW - T0).total_seconds()

FIN = json.load(open(OUT / "rw2_finished.json")) if (OUT / "rw2_finished.json").exists() else {}

def finished_block(r, f):
    cmp_rows = "".join(f"<tr><td>{c['at_mm']:.1f} mm · {c['at_mm']/36.6:.2f}%</td><td>{c['model_kN']:.1f}</td>"
                       f"<td>{(c.get('test_envelope_kN') or float('nan')):.1f}</td><td>{(c.get('model_over_test') or float('nan')):.3f}</td></tr>"
                       for c in f["comparison"] if c.get("test_envelope_kN"))
    lp_rows = "".join(f"<tr><td>{a:.1f}%</td><td>{v:.0f}%</td><td>{rb:.0f}%</td><td>{dg:.0f}%</td></tr>" for a, dg, rb, v in f["loadpath"])
    cap = (f"{f['capacity']*100:.2f}% on the 5 ms curve — a ringing dip; windows of 10 ms and longer never cross 80%, and the run ends at {f['end_shear']/f['peak']:.3f} × peak"
           if f.get("capacity") else f"not reached — ends at {f['end_shear']/f['peak']:.3f} × peak")
    return f'''
    <article class="run" style="--c:{r['color']}">
      <header class="run-h">
        <div><div class="run-name">{r['name']}</div><div class="run-sub">{r['sub']} · {r['elements']}</div></div>
        <span class="chip ok"><i></i>finished · {f['hours']:.1f} h</span>
      </header>
      <div class="bar"><div class="fill" style="width:100%"></div></div>
      <dl class="kv">
        <div><dt>peak, {f['window_ms']:g} ms smoothed</dt><dd>{f['peak']:.1f} kN · <b>{f['peak']/TEST:.3f} × test</b></dd></div>
        <div><dt>at drift</dt><dd>{f['peak_drift']:.2f}% · plateau {f['plateau'][0]:.2f}–{f['plateau'][1]:.2f}%</dd></div>
        <div><dt>raw maximum</dt><dd>{f['raw']:.1f} kN · {f['ring']:.2f} × — ringing, not resistance</dd></div>
        <div><dt>end of run</dt><dd>{f['end_shear']:.1f} kN at {f['end_drift']:.2f}% · {f['end_shear']/TEST:.3f} × test</dd></div>
        <div><dt>residual, ascending p95</dt><dd>{f['residual_pct']:.1f}% of peak — the peak is real resistance</dd></div>
        <div><dt>drift capacity (80% drop)</dt><dd>{cap}</dd></div>
      </dl>
      <h3>at matched displacement</h3>
      <div class="scroll"><table><thead><tr><th>displacement</th><th>model kN</th><th>test envelope</th><th>model / test</th></tr></thead><tbody>{cmp_rows}</tbody></table></div>
      <h3>share of the overturning moment across the base cut</h3>
      <div class="scroll"><table><thead><tr><th>drift</th><th>vertical struts</th><th>rebar</th><th>diagonals</th></tr></thead><tbody>{lp_rows}</tbody></table></div>
      <p>The concrete verticals shed from 44% to 32% while the bars pick up 28% → 39%: the compressed boundary element passes its peak strain and hands its share to the steel, which at nominal f<sub>y</sub> without hardening cannot grow. That is the 138 → 125 kN step at 0.6%. Only 66 of 22,002 struts are past ε<sub>c0</sub> at the end, all in the compression toe.</p>
      <img src="{f['fig_comparison']}" alt="this run against the digitized record and Aydin's curve" style="border-radius:6px;border:1px solid var(--line)">
      <img src="{f['fig_damage']}" alt="damage pattern at peak and at the end" style="border-radius:6px;border:1px solid var(--line)">
    </article>'''

def run_block(r):
    if r["key"] in FIN:
        return finished_block(r, FIN[r["key"]])
    i, d, v, t = r["samples"][-1]; rate = i / t; total_h = r["n"] / rate / 3600
    pct = min(100, (i + rate * (el - t)) / r["n"] * 100); eta = T0 + datetime.timedelta(hours=total_h)
    peak = max(s[2] for s in r["samples"]); pk_d = [s[1] for s in r["samples"] if s[2] == peak][0]
    rows = "".join(f"<tr><td>{s[1]:.3f}%</td><td>{s[2]:.1f}</td><td>{s[2]/TEST:.3f}</td>"
                   f"<td class='dim'>{s[0]//1000:,}k</td></tr>" for s in r["samples"])
    return f'''
    <article class="run" style="--c:{r['color']}">
      <header class="run-h">
        <div><div class="run-name">{r['name']}</div><div class="run-sub">{r['sub']} · {r['elements']}</div></div>
        <span class="chip live"><i></i>running</span>
      </header>
      <div class="bar"><div class="fill" style="width:{pct:.1f}%"></div></div>
      <dl class="kv">
        <div><dt>progress</dt><dd>{pct:.0f}% · step {i:,} of {r['n']:,}</dd></div>
        <div><dt>rate</dt><dd>{rate:.0f} steps/s</dd></div>
        <div><dt>eta</dt><dd>{eta:%a %H:%M}</dd></div>
        <div><dt>peak so far</dt><dd>{peak:.1f} kN · {peak/TEST:.3f} × test · at {pk_d:.2f}%</dd></div>
      </dl>
      <div class="scroll"><table><thead><tr><th>drift</th><th>base shear kN</th><th>/ test</th><th>step</th></tr></thead><tbody>{rows}</tbody></table></div>
    </article>'''

def sw(kind, color):
    if kind == "run":
        return f'<svg width="34" height="12" viewBox="0 0 34 12"><line x1="1" y1="6" x2="33" y2="6" stroke="{color}" stroke-width="2.2"/><circle cx="8" cy="6" r="3" fill="{color}"/><circle cx="26" cy="6" r="3" fill="{color}"/></svg>'
    if kind == "dash":
        return f'<svg width="34" height="12" viewBox="0 0 34 12"><line x1="1" y1="6" x2="33" y2="6" stroke="{color}" stroke-width="1.8" stroke-dasharray="5 4"/></svg>'
    if kind == "dot":
        return f'<svg width="34" height="12" viewBox="0 0 34 12"><line x1="1" y1="6" x2="33" y2="6" stroke="{color}" stroke-width="1.6" stroke-dasharray="2 3"/></svg>'
    if kind == "fine":
        return f'<svg width="34" height="12" viewBox="0 0 34 12"><line x1="1" y1="6" x2="33" y2="6" stroke="{color}" stroke-width="1.2" stroke-dasharray="1 3"/></svg>'
    if kind == "cloud":
        return f'<svg width="34" height="12" viewBox="0 0 34 12">' + "".join(f'<rect x="{x}" y="{y}" width="1.6" height="1.6" fill="{color}"/>' for x, y in ((3,4),(7,8),(11,3),(15,7),(19,5),(23,9),(27,4),(31,7))) + '</svg>'

LEGEND = (
    '<div class="legend">'
    + f'<span>{sw("run","var(--uni)")}<b>uniform 30.5 grid</b> — this study, console samples</span>'
    + f'<span>{sw("run","var(--gra)")}<b>graded 25 grid</b> — this study, console samples</span>'
    + f'<span>{sw("dash","var(--test)")}<b>test envelope</b> — running maximum of the digitized Fig. 9(b) loops</span>'
    + f'<span>{sw("fine","var(--test)")}<b>163.3 kN</b> — measured peak, Table 4</span>'
    + f'<span>{sw("dot","var(--aydin)")}<b>Aydin 2019</b> — his lattice, horizon 1.5d, 169.8 kN</span>'
    + f'<span>{sw("cloud","var(--faint)")}<b>test loops</b> — the digitized cyclic record</span>'
    + '</div>')

uni, gra = STATE["runs"]
u_peak = max(s[2] for s in uni["samples"]); u_last = uni["samples"][-1]
html = f'''<title>RW2 Run Board</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Archivo:wght@500;600;700&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root{{
  --ground:#EEF1F3; --surface:#FFFFFF; --sunken:#E3E8EC; --ink:#171B20; --muted:#5A6470; --faint:#8791.9C;
  --faint:#87919C; --line:#D5DBE1; --line-strong:#B9C2CB;
  --test:#A63D2F; --test-soft:#A63D2F1A; --gra:#1E5A86; --uni:#2E7D5B; --aydin:#7A6AA6;
  --run:#B27A19; --run-soft:#B27A191C; --ok:#1F6E4A; --ok-soft:#1F6E4A18; --warn:#A63D2F;
  --display:"Archivo",-apple-system,BlinkMacSystemFont,sans-serif; --sans:"IBM Plex Sans",-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; --mono:"IBM Plex Mono",ui-monospace,Menlo,monospace;
}}
@media (prefers-color-scheme:dark){{ :root:not([data-theme="light"]){{
  --ground:#101418; --surface:#181E24; --sunken:#1F262E; --ink:#E6EAEE; --muted:#98A2AD; --faint:#6C7681;
  --line:#2A333C; --line-strong:#3D4854; --test:#D9705C; --test-soft:#D9705C22; --gra:#6FA8D6; --uni:#63B58E; --aydin:#B0A3D6;
  --run:#D6A64A; --run-soft:#D6A64A22; --ok:#5DB58A; --ok-soft:#5DB58A1E; --warn:#D9705C; }} }}
:root[data-theme="dark"]{{
  --ground:#101418; --surface:#181E24; --sunken:#1F262E; --ink:#E6EAEE; --muted:#98A2AD; --faint:#6C7681;
  --line:#2A333C; --line-strong:#3D4854; --test:#D9705C; --test-soft:#D9705C22; --gra:#6FA8D6; --uni:#63B58E; --aydin:#B0A3D6;
  --run:#D6A64A; --run-soft:#D6A64A22; --ok:#5DB58A; --ok-soft:#5DB58A1E; --warn:#D9705C; }}
*{{box-sizing:border-box}}
body{{background:var(--ground); color:var(--ink); font:15px/1.55 var(--sans); -webkit-font-smoothing:antialiased; margin:0}}
.wrap{{max-width:700px; margin:0 auto; padding:26px 18px 56px; display:flex; flex-direction:column; gap:28px}}
.eyebrow{{font:500 11px/1 var(--mono); letter-spacing:.14em; text-transform:uppercase; color:var(--faint)}}
h1{{font:700 30px/1.1 var(--display); letter-spacing:-.02em; margin:6px 0 0; text-wrap:balance}}
h2{{font:600 18px/1.2 var(--display); letter-spacing:-.01em; margin:0}}
h3{{font:500 11px/1 var(--mono); letter-spacing:.09em; text-transform:uppercase; color:var(--faint); margin:0; padding-top:4px}}
.stamp{{font:12px/1.5 var(--mono); color:var(--faint); border-top:1px solid var(--line); padding-top:8px; margin-top:10px}}
.sect{{display:flex; flex-direction:column; gap:12px}}
p{{margin:0; color:var(--muted); font-size:14px; max-width:64ch}} p b{{color:var(--ink); font-weight:600}}
.tiles{{display:grid; grid-template-columns:repeat(3,1fr); gap:10px}} @media (max-width:560px){{.tiles{{grid-template-columns:1fr 1fr}}}}
.tile{{background:var(--surface); border:1px solid var(--line); border-radius:6px; padding:12px 13px; display:flex; flex-direction:column; gap:4px}}
.tile .k{{font:500 10.5px/1 var(--mono); letter-spacing:.08em; text-transform:uppercase; color:var(--faint)}}
.tile .v{{font:600 22px/1.1 var(--display); font-variant-numeric:tabular-nums; letter-spacing:-.01em}}
.tile .s{{font:12px/1.4 var(--sans); color:var(--muted)}}
.tile.test{{border-color:var(--test); background:var(--test-soft)}}
.chip{{display:inline-flex; align-items:center; gap:6px; font:500 10.5px/1 var(--mono); letter-spacing:.06em; text-transform:uppercase; padding:4px 8px; border-radius:3px; white-space:nowrap}}
.chip i{{width:6px; height:6px; border-radius:50%; background:currentColor}}
.chip.live{{background:var(--run-soft); color:var(--run)}} .chip.ok{{background:var(--ok-soft); color:var(--ok)}} .chip.todo{{background:var(--sunken); color:var(--faint)}}
.run{{background:var(--surface); border:1px solid var(--line); border-left:4px solid var(--c); border-radius:6px; padding:14px 16px; display:flex; flex-direction:column; gap:11px}}
.run-h{{display:flex; justify-content:space-between; align-items:flex-start; gap:10px}}
.run-name{{font:600 15px/1.2 var(--display)}} .run-sub{{font:12px/1.4 var(--sans); color:var(--muted); margin-top:2px}}
.bar{{height:6px; background:var(--sunken); border-radius:3px; overflow:hidden}} .fill{{height:100%; background:var(--c); border-radius:3px}}
.kv{{display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:8px 14px; margin:0}}
.kv div{{display:flex; flex-direction:column; gap:1px}} .kv dt{{font:500 10px/1 var(--mono); letter-spacing:.08em; text-transform:uppercase; color:var(--faint)}}
.kv dd{{margin:0; font:13.5px/1.3 var(--mono); font-variant-numeric:tabular-nums}}
.scroll{{overflow-x:auto; border:1px solid var(--line); border-radius:4px; background:var(--surface)}}
table{{border-collapse:collapse; width:100%; min-width:360px}}
th,td{{text-align:right; padding:7px 10px; font:12px/1.3 var(--mono); font-variant-numeric:tabular-nums; white-space:nowrap; border-bottom:1px solid var(--line)}}
th{{font-size:10px; letter-spacing:.07em; text-transform:uppercase; color:var(--faint); font-weight:500; background:var(--sunken)}}
th:first-child,td:first-child{{text-align:left}} tbody tr:last-child td{{border-bottom:0}} td.dim{{color:var(--faint)}} td.hi{{font-weight:600}}
canvas{{width:100%; display:block; background:var(--surface); border:1px solid var(--line); border-radius:6px}}
.legend{{display:grid; grid-template-columns:1fr 1fr; gap:6px 16px; font:12px/1.35 var(--sans); color:var(--muted); padding:2px 2px 0}} .legend span{{display:flex; align-items:center; gap:8px}} .legend svg{{flex:none}} .legend b{{color:var(--ink); font-weight:600}} @media (max-width:520px){{.legend{{grid-template-columns:1fr}}}}
.ln{{width:18px; height:0; border-top:2px solid var(--c)}} .ln.dash{{border-top-style:dashed}} .ln.dot{{border-top-style:dotted}}
.gates{{border:1px solid var(--line); border-radius:6px; background:var(--surface)}}
.gate{{display:flex; justify-content:space-between; align-items:center; gap:12px; padding:10px 14px; border-bottom:1px solid var(--line)}} .gate:last-child{{border-bottom:0}}
.gate .t{{font-size:13.5px; font-weight:500}} .gate .d{{font-size:12px; color:var(--muted)}} .gate .v{{font:500 14px/1 var(--mono); font-variant-numeric:tabular-nums; text-align:right; flex:none}}
.v.ok{{color:var(--ok)}} .v.warn{{color:var(--warn)}}
.plan{{display:flex; flex-direction:column; border:1px solid var(--line); border-radius:6px; background:var(--surface)}}
.step{{display:grid; grid-template-columns:96px 1fr auto; gap:12px; align-items:center; padding:10px 14px; border-bottom:1px solid var(--line)}} .step:last-child{{border-bottom:0}}
.step .n{{font:500 11px/1 var(--mono); letter-spacing:.08em; text-transform:uppercase; color:var(--faint)}} .step .w{{font-size:13.5px}} .step .w small{{display:block; font-size:12px; color:var(--muted)}}
.find{{display:flex; flex-direction:column; gap:12px}} .f{{padding-left:12px; border-left:2px solid var(--line-strong); display:flex; flex-direction:column; gap:2px}} .f.alert{{border-left-color:var(--test)}}
.f .t{{font-size:14px; font-weight:600}} .f .b{{font-size:13.5px; color:var(--muted)}}
footer{{font:12px/1.6 var(--mono); color:var(--faint); border-top:1px solid var(--line); padding-top:12px}}
</style>
<div class="wrap">
  <header>
    <div class="eyebrow">Thomsen &amp; Wallace RW2 · lattice study</div>
    <h1>RW2 run board</h1>
    <div class="stamp">{NOW:%H:%M · %Y-%m-%d} — {STATE['stamp_note']}</div>
  </header>

  <section class="sect">
    <h2>Stage 1 · complete</h2>
    <div class="tiles">
      <div class="tile"><span class="k">uniform 30.5 · finished</span><span class="v">{FIN['uniform']['peak']/TEST:.3f} ×</span><span class="s">peak {FIN['uniform']['peak']:.1f} kN at {FIN['uniform']['peak_drift']:.2f}% · to 2.5%, converged</span></div>
      <div class="tile"><span class="k">graded 25 · finished</span><span class="v">{FIN['graded']['peak']/TEST:.3f} ×</span><span class="s">peak {FIN['graded']['peak']:.1f} kN at {FIN['graded']['peak_drift']:.2f}% · to 2.5%, converged</span></div>
      <div class="tile test"><span class="k">measured · Aydin</span><span class="v">163.3</span><span class="s">kN, Table 4 · Aydin's own lattice 1.040 ×</span></div>
    </div>
    <p>Both runs are the same cell — <b>crushing / solved / perfect bond</b>, Concrete02, b = 0.01, f<sub>y</sub> 414 (nominal), no rupture switch, ζ = 0.5, 7.6 mm/s, explicit CentralDifference, target <b>2.5% drift</b> — and differ only in the grid. Both runs are finished and scored from their stored series (5 ms smoothing — T1 is 20 ms here, D106). Nothing is running.</p>
  </section>

  <section class="sect">
    <h2>Against the record</h2>
    <canvas id="cFull" style="height:300px" aria-label="Both runs over the digitized Fig. 9(b) record, full range"></canvas>
    {LEGEND}
    <h3>the range reached so far</h3>
    <canvas id="cZoom" style="height:300px" aria-label="Both runs over the record, zoomed to the drift reached so far"></canvas>
    {LEGEND}
    <p><b>The two grids give one answer.</b> Peaks 139.4 (uniform) and 138.3 kN (graded), 0.854 and 0.847 × the measured; at matched drift the graded/uniform ratio runs 0.98–1.05 from 0.3% to 2.5%. So the peak, the 0.6% step and the post-peak slide are properties of the model, not of the discretization — the same conclusion Aldemir's mesh 25/50 pair gave (D98). Both are the flexural-yield plateau the section check predicted (134 kN at nominal f<sub>y</sub>, no hardening), and both fall behind a test that keeps rising to 163 kN at ~1.5%.</p>
  </section>

  <section class="sect">
    <h2>Runs</h2>
    {run_block(uni)}
    {run_block(gra)}
  </section>

  <section class="sect">
    <h2>Stage 0 · gates</h2>
    <div class="gates">
      <div class="gate"><div><div class="t">Harness regression after the lift (D103)</div><div class="d">master report, payload, 2 smoke runs, 38 run sheets — 45 Aldemir runs before vs after</div></div><span class="v ok">4 / 4 identical</span></div>
      <div class="gate"><div><div class="t">Fig. 9(b) digitized</div><div class="d">envelope peak vs Table 4; the two simulated plateaus 0.994 / 1.018; record unclipped</div></div><span class="v ok">1.003 ×</span></div>
      <div class="gate"><div><div class="t">K<sub>lattice</sub> / K<sub>continuum</sub> — graded grid</div><div class="d">same grid, rebar, BCs. Aldemir gave 0.999; here the D53 flexural under-read (ν<sub>eff</sub> 0.41)</div></div><span class="v">0.951</span></div>
      <div class="gate"><div><div class="t">K<sub>lattice</sub> / K<sub>continuum</sub> — uniform 30.5</div><div class="d">graded is no worse, and puts every bar on its axis</div></div><span class="v">0.941</span></div>
      <div class="gate"><div><div class="t">Published equibiaxial field</div><div class="d">worse again, as on Aldemir — the uniaxial field stays</div></div><span class="v">0.821</span></div>
      <div class="gate"><div><div class="t">Table 4 “initial stiffness” 35.19 kN/mm</div><div class="d">vs the uncracked transformed cantilever 29.9 — impossible for a top-displacement secant (cf. D73); Aydin's 32.26 also sits above it</div></div><span class="v warn">1.18 × section</span></div>
      <div class="gate"><div><div class="t">Strut life at 25 mm</div><div class="d">orthogonal / diagonal — no G<sub>f</sub> factor needed</div></div><span class="v ok">45.3 / 32.4</span></div>
      <div class="gate"><div><div class="t">V<sub>flex</sub> at nominal f<sub>y</sub>, no hardening, N included</div><div class="d">the prediction on record: the baseline lands below the test</div></div><span class="v">0.82 × test</span></div>
    </div>
  </section>

  <section class="sect">
    <h2>The specimen, as the 2019 paper reports it</h2>
    <div class="scroll"><table>
      <thead><tr><th>quantity</th><th>value</th><th>source</th></tr></thead>
      <tbody>
        <tr><td>peak force</td><td class="hi">163.284 kN</td><td class="dim">Table 4</td></tr>
        <tr><td>initial stiffness</td><td>35.19 kN/mm</td><td class="dim">Table 4 — see gate</td></tr>
        <tr><td>max displacement</td><td>72 mm ≈ 1.97% drift</td><td class="dim">Table 2</td></tr>
        <tr><td>axial load</td><td>378 kN = 0.071 A<sub>g</sub>f<sub>c</sub></td><td class="dim">Table 1</td></tr>
        <tr><td>geometry</td><td>1220 × 3660 × 102 · aspect 3.0</td><td class="dim">Fig. 9(a)</td></tr>
        <tr><td>boundary elements</td><td>8-#3 each end · hoops 4.76 @ 76</td><td class="dim">Fig. 9(a)</td></tr>
        <tr><td>web</td><td>8-#2 vertical · #2 @ 191 horizontal</td><td class="dim">Fig. 9(a)</td></tr>
        <tr><td>f<sub>c</sub> / f<sub>t</sub> / E<sub>c</sub></td><td>42.8 / 2.03 / 31,030 MPa</td><td class="dim">Table 1</td></tr>
        <tr><td>f<sub>y</sub></td><td>414 MPa — nominal Grade 60</td><td class="dim">Table 1</td></tr>
        <tr><td class="dim">Aydin, horizon 1.5</td><td class="dim">169.834 kN · 32.26 kN/mm · 1.040 ×</td><td class="dim">Table 4</td></tr>
        <tr><td class="dim">Aydin, horizon 3.01</td><td class="dim">174.544 kN · 38.38 kN/mm</td><td class="dim">Table 4</td></tr>
      </tbody></table></div>
    <p>Neither primary source is in the repo: every number is second-hand through Aydin (2019). Not printed there: the loading levels, coupon strengths, hardening, rupture strain, failure mode.</p>
  </section>

  <section class="sect">
    <h2>Plan</h2>
    <div class="plan">
      <div class="step"><span class="n">stage 0</span><span class="w">Elastic gates<small>continuum, cantilever, both calibration fields, both grids</small></span><span class="chip ok"><i></i>done</span></div>
      <div class="step"><span class="n">stage 1</span><span class="w">Baseline pushover to 2.5%<small>crushing / solved / perfect — uniform 0.854 ×, graded 0.847 ×; grid-objective to 5% at every drift</small></span><span class="chip ok"><i></i>done</span></div>
      <div class="step"><span class="n">stage 2</span><span class="w">The matrix<small>compression law × tension tail; priced, parallel</small></span><span class="chip todo">next</span></div>
      <div class="step"><span class="n">stage 3</span><span class="w">Failure model<small>--steel-rupture × --concrete-residual 0, then b × ε<sub>su</sub>; likely also f<sub>y</sub> and a confined boundary grade</small></span><span class="chip todo">queued</span></div>
      <div class="step"><span class="n">stage 4</span><span class="w">Cyclic<small>on the cell that reproduces the monotonic peak; levels invented, said so</small></span><span class="chip todo">queued</span></div>
    </div>
  </section>

  <section class="sect">
    <h2>Findings so far</h2>
    <div class="find">
      <div class="f alert"><span class="t">Peak here is a yield quantity, not a cracking one</span><span class="b">V<sub>flex</sub> at nominal f<sub>y</sub> is 0.82 × the test and the run plateaus at 0.85 ×. The f<sub>y</sub> / b / ε<sub>su</sub> axes will matter on this wall the way the compression law did not on Aldemir.</span></div>
      <div class="f"><span class="t">Grid-objective (D107)</span><span class="b">Graded 25 mm (bars on nodes, 32,423 elements, 33.2 h) against uniform 30.5 (bars snapped, 22,002 elements, 13.3 h): peak 0.992 ×, matched-drift ratio 0.98–1.05 throughout, residual 0.6% vs 0.9%. The graded grid buys exact bar placement at 2.5 × the cost and changes no answer; the matrix can run on the uniform grid.</span></div>
      <div class="f alert"><span class="t">Both grids, final: 0.85 × the test, and the model is stiff-then-weak</span><span class="b">139.4 kN at 0.50% against 163.3 measured; at 0.5% drift the model is 1.09 × the test envelope, at 1.5% it is 0.81 ×. No hardening and a toe that sheds to the bars at 0.6% against a test that keeps gaining to 1.5% on hardening and a confined boundary.</span></div>
      <div class="f alert"><span class="t">Early softening under axial load</span><span class="b">7% off the peak by 0.64% drift. First candidate: the compressed boundary element passes ε<sub>c0</sub> = 0.00276 early under 378 kN on a 102 mm wall, and Concrete02 softens to 0.2 f<sub>c</sub> while the real boundary is confined by hoops at 76 mm — no confined grade here (the D66 reasoning was written for a 150 mm wall). Sampled, not yet measured.</span></div>
      <div class="f"><span class="t">5% soft against the continuum, by construction</span><span class="b">The uniaxial-field balance pins E/(1−ν²); bending reads C<sub>11</sub>(1−ν<sub>eff</sub>²) with the lattice's own ν<sub>eff</sub> = 0.41. Recorded as a property of the method, not fitted — a property Aldemir's squat wall hid.</span></div>
      <div class="f"><span class="t">Bars on nodes without a common divisor (D104)</span><span class="b"><code>--grid rebar</code>: bar axes become grid lines, gaps fill at the target, horizon in index space, strut areas scale with tributary width. Costs 2.4 × the uniform grid per run; buys exact bar placement, not a different answer so far.</span></div>
    </div>
  </section>

  <section class="sect">
    <h2>Placed against the repo</h2>
    <div class="scroll"><table>
      <thead><tr><th>specimen</th><th>aspect</th><th>K / K<sub>cont</sub></th><th>peak / test</th><th>bond</th></tr></thead>
      <tbody>
        <tr><td>SW-NC-FF</td><td>3.00</td><td class="dim">0.93 vs transformed</td><td>1.08 / 1.05 cyclic</td><td class="dim">plain bars — unfair</td></tr>
        <tr><td>WSH3</td><td>2.28</td><td class="dim">0.93 vs transformed</td><td>0.950 cyclic</td><td class="dim">deformed — fair</td></tr>
        <tr><td>Aldemir</td><td>0.75</td><td>0.999</td><td>0.997 push · 0.966 cyclic</td><td class="dim">not reported</td></tr>
        <tr><td class="hi">RW2</td><td>3.00</td><td class="hi">0.951</td><td class="hi">0.854 / 0.847 push (uniform / graded)</td><td class="dim">not reported</td></tr>
      </tbody></table></div>
  </section>

  <footer>Owed when the runs land: smoothed peak, ascending-branch residual, load-path split, damage figure; then master report, advisor page, run sheets, DECISIONS. Source: <code>examples/thomsen_wallace_wall/</code>, runs under <code>examples/output/thomsen_wallace_wall/study/</code>.</footer>
</div>
<script>
(function(){{
  var R={json.dumps(REF)}, TEST={TEST};
  var RUNS={json.dumps([{ "c": r["color"], "pts": [[s[1], s[2]] for s in r["samples"]] } for r in STATE["runs"]])};
  function css(v){{return getComputedStyle(document.documentElement).getPropertyValue(v.slice(4,-1)).trim();}}
  function draw(id,xlim,ylim){{
    var c=document.getElementById(id); if(!c) return; var W=c.clientWidth||660,Hh=c.clientHeight||300;
    c.width=W*2; c.height=Hh*2; var g=c.getContext('2d'); g.scale(2,2);
    var ink=css('var(--ink)'),faint=css('var(--faint)'),line=css('var(--line)'),test=css('var(--test)'),ayd=css('var(--aydin)');
    var L=46,Rr=12,T=14,B=32, X=function(x){{return L+(x-xlim[0])/(xlim[1]-xlim[0])*(W-L-Rr)}}, Y=function(y){{return T+(ylim[1]-y)/(ylim[1]-ylim[0])*(Hh-T-B)}};
    g.font='10px IBM Plex Mono, monospace'; g.fillStyle=faint; g.strokeStyle=line; g.lineWidth=1;
    for(var y=ylim[0]; y<=ylim[1]+1e-9; y+=50){{ g.beginPath(); g.moveTo(L,Y(y)); g.lineTo(W-Rr,Y(y)); g.stroke(); g.fillText(y,4,Y(y)+3); }}
    var xs=(xlim[1]-xlim[0])>1?0.5:0.1; for(var x=Math.ceil(xlim[0]/xs-1e-9)*xs; x<=xlim[1]+1e-9; x+=xs){{ g.beginPath(); g.moveTo(X(x),T); g.lineTo(X(x),Hh-B); g.stroke(); g.fillText(x.toFixed(1)+'%',X(x)-11,Hh-B+14); }}
    g.fillStyle=faint; g.globalAlpha=.55; R.cloud.forEach(function(p){{ if(p[0]>=xlim[0]&&p[0]<=xlim[1]&&p[1]>=ylim[0]&&p[1]<=ylim[1]) g.fillRect(X(p[0])-.7,Y(p[1])-.7,1.4,1.4); }}); g.globalAlpha=1;
    function poly(pts,col,w,dash){{ g.strokeStyle=col; g.lineWidth=w; g.setLineDash(dash||[]); g.beginPath(); var f=true; pts.forEach(function(p){{ if(p[0]<xlim[0]||p[0]>xlim[1]) return; var px=X(p[0]),py=Y(Math.max(ylim[0],Math.min(ylim[1],p[1]))); if(f){{g.moveTo(px,py);f=false;}} else g.lineTo(px,py); }}); g.stroke(); g.setLineDash([]); }}
    poly(R.env,test,1.6,[5,4]); poly(R.aydin,ayd,1.4,[2,3]);
    g.strokeStyle=test; g.setLineDash([1,3]); g.beginPath(); g.moveTo(L,Y(TEST)); g.lineTo(W-Rr,Y(TEST)); g.stroke(); g.setLineDash([]);
    RUNS.forEach(function(r){{ var col=css(r.c); poly(r.pts,col,2.2); g.fillStyle=col; r.pts.forEach(function(p){{ if(p[0]>=xlim[0]&&p[0]<=xlim[1]){{ g.beginPath(); g.arc(X(p[0]),Y(p[1]),3.2,0,6.3); g.fill(); }} }}); }});
    g.fillStyle=ink; g.font='11px IBM Plex Sans, sans-serif'; g.fillText('base shear (kN)  vs  drift (%) = top displacement / 3660 mm', L+6, T+12);
  }}
  function all(){{ draw('cFull',[-2.1,2.1],[-200,200]); draw('cZoom',[0,0.8],[0,200]); }}
  all(); window.addEventListener('resize',all);
  if(window.matchMedia) window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change',all);
}})();
</script>
'''
(OUT / "rw2_board.html").write_text(html)
print(len(html))
