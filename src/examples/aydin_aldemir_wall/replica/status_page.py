"""Write a self-contained status page for the running cyclic analysis (figure embedded).

    uv run python examples/aydin_aldemir_wall/replica/status_page.py

Regenerates report.png, then writes an HTML file with the figure inlined as a data URI, so the
page needs no network and can be published as an artifact and opened on a phone.
"""
from __future__ import annotations

import base64, datetime, json, re, subprocess, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent.parent / "output" / "aydin_aldemir_wall" / "replica"
LOG = OUT / "_live_cyclic.log"
HW = 2680.0


def main() -> None:
    subprocess.run([sys.executable, str(HERE / "report.py")], cwd=HERE.parent.parent.parent,
                   capture_output=True)
    run = max((d for d in OUT.glob("*cyclic*") if d.is_dir()), key=lambda d: d.stat().st_mtime)
    subprocess.run([sys.executable, str(HERE / "damage_extra.py")],
                   cwd=HERE.parent.parent.parent, capture_output=True)
    b64 = base64.b64encode((run / "report.png").read_bytes()).decode()

    def embed(name):
        f = run / name
        return base64.b64encode(f.read_bytes()).decode() if f.exists() else None

    b64_rev = embed("damage_reversals.png")
    b64_evo = embed("damage_evolution.png")

    r = re.compile(r"step\s+(\d+)/(\d+).*\[\s*(\d+)s\]")
    p = [(int(m.group(1)), int(m.group(2)), int(m.group(3)))
         for l in LOG.read_text().splitlines() if (m := r.search(l))]
    (i0, n, t0), (i1, _n, t1) = p[1], p[-1]
    sps = (t1 - t0) / (i1 - i0)
    rem = (n - i1) * sps
    eta = datetime.datetime.now() + datetime.timedelta(seconds=rem)
    running = subprocess.run(["pgrep", "-f", "replica/run.py --cyclic"],
                             capture_output=True).returncode == 0

    par = json.loads((run / "params.json").read_text())
    z = np.load(run / "hysteresis.npz")
    u, s = z["disp"], z["shear"] / 1e3
    d = np.sign(np.diff(u))
    turns = np.nonzero(np.diff(d))[0] + 1
    rev = [(u[k], s[k]) for k in turns]
    levels = [("0.050%", 1.34), ("0.100%", 2.68), ("0.150%", 4.02)]
    rows = ""
    for li, (name, mm) in enumerate(levels):
        pos = next((f"{v:+,.0f}" for x, v in rev if abs(x - mm) < .08), "—")
        neg = next((f"{v:+,.0f}" for x, v in rev if abs(x + mm) < .08), "—")
        ratio = "—"
        if pos != "—" and neg != "—":
            ratio = f"{abs(float(neg.replace(',','')) / float(pos.replace(',',''))):.3f}"
        rows += (f"<tr><td>{name}</td><td>{mm:.2f}</td><td>{pos}</td>"
                 f"<td>{neg}</td><td>{ratio}</td></tr>")

    # --- parameter tables -------------------------------------------------------------------
    import os
    os.environ.setdefault("ALDEMIR_TW", "210")
    sys.path.insert(0, str(HERE))
    import specimen as sp
    from build import calibrate as _cal
    at = {h: _cal(horizon=h).area for h in (1.5, 3.01)}
    tw = par.get("thickness", sp.TW)

    def rows_of(pairs):
        return "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in pairs)

    shared = rows_of([
        ("Panel (from his Table 2 counts)", f"{sp.LW:.0f} &times; {sp.HW:.0f} mm"),
        ("Thickness", f"{tw:.0f} mm <span class='q'>never published; set by K_sim</span>"),
        ("Grid / horizon", f"{sp.MESH:.0f} mm &middot; 1.5 and 3.01"),
        ("Concrete", f"f<sub>c</sub> {sp.FC:.0f}, f<sub>t</sub> {sp.FT:.2f}, "
                     f"E {sp.EC:,.0f} MPa, G<sub>f</sub> {sp.GF*1e3:.0f} N/m"),
        ("Poisson", f"{sp.NU:.3f} <span class='q'>his Appendix</span>"),
        ("Steel", f"f<sub>y</sub> {sp.FY:.0f} MPa, E {sp.ES/1e3:,.0f} GPa, "
                  f"elastic&ndash;perfectly plastic"),
        ("Reinforcement", f"3&#8211;&#216;8 @ {sp.S_BAR:.0f} mm both ways "
                          f"({sp.BAR_AREA:.0f} mm&sup2;/line), run to the top"),
        ("Calibration", f"equibiaxial energy balance &middot; A<sub>t</sub> "
                        f"{at[1.5]:,.0f} mm&sup2; (h 1.5), {at[3.01]:,.0f} mm&sup2; (h 3.01)"),
        ("Tension tail", f"a1 = horizon, a2 {sp.A2:.0f}, a3 {sp.A3:.0f}, "
                         f"b1 {sp.B1}, b2 {sp.B2} <span class='q'>b1/b2 not published</span>"),
        ("Axial load", "none"),
    ])

    def run_row(name, hz, law, drive, nodes, els, res):
        return (f"<tr><td>{name}</td><td>{hz}</td><td>{law}</td><td>{drive}</td>"
                f"<td>{nodes}</td><td>{els}</td><td>{res}</td></tr>")

    lvl_txt = ", ".join(f"{float(v):.3%}" for v in str(par.get("cyclic", "")).split(",") if v)
    law_name = {"concrete02": "Concrete02", "aydin": "HystereticSM"}
    ours = run_row(f"cyclic <span class='q'>(this run)</span>", "1.5",
                   law_name.get(par.get("cyclic_law", "aydin"), "?"),
                   f"&plusmn;{lvl_txt}, 1 cyc", "20,385", "88,474",
                   f"+{max(s):,.0f} / {min(s):,.0f} kN")
    for d in sorted(OUT.glob("*cyclic*")):
        if d == run or not (d / "hysteresis.npz").exists():
            continue
        zz = np.load(d / "hysteresis.npz"); ss = zz["shear"] / 1e3
        pj = json.loads((d / "params.json").read_text())
        ours += run_row("cyclic <span class='q'>(stopped &mdash; no pinching)</span>", "1.5",
                        law_name.get(pj.get("cyclic_law", "aydin"), "?"),
                        f"&plusmn;{lvl_txt}, 1 cyc", "20,385", "88,474",
                        f"+{max(ss):,.0f} / {min(ss):,.0f} kN")
    ours += (
             run_row("pushover", "1.5", "ElasticMultiLinear", "monotonic to 0.30%",
                      "20,385", "88,474", "897 kN peak")
            + run_row("pushover", "3.01", "ElasticMultiLinear", "monotonic to 0.30%",
                      "20,385", "288,050", "&ge;861 kN, no peak"))

    refs = rows_of([
        ("Aydin lattice, h 1.5", "K 943.2 kN/mm &middot; F 1,164.4 kN &middot; 80,684 struts"),
        ("Aydin lattice, h 3.01", "K 1,052.4 kN/mm &middot; F 1,325.7 kN &middot; 280,260 struts"),
        ("Aydin drive", "PID force control, dt 5&times;10<sup>-8</sup> s "
                        "<span class='q'>we prescribe displacement</span>"),
        ("Measured test", "K 1,038.4 kN/mm &middot; F 963.6 kN at ~20 mm, no degradation"),
        ("Test protocol", "<span class='q'>never published, and not recoverable from Fig. 10(b)</span>"),
    ])

    html = f"""<title>Aldemir Cyclic Watch</title>
<style>
:root {{
  --ground:#f6f7f8; --surface:#ffffff; --ink:#16202b; --muted:#5d6b78;
  --line:#dde3e8; --accent:#b03a2e; --support:#16786b; --bar:#e6eaee;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --ground:#121820; --surface:#1a222c; --ink:#e8eef4; --muted:#93a3b2;
    --line:#2a3640; --accent:#e2705f; --support:#4cbfae; --bar:#26313c;
  }}
}}
:root[data-theme="dark"] {{
  --ground:#121820; --surface:#1a222c; --ink:#e8eef4; --muted:#93a3b2;
  --line:#2a3640; --accent:#e2705f; --support:#4cbfae; --bar:#26313c;
}}
* {{ box-sizing:border-box; }}
body {{
  margin:0; background:var(--ground); color:var(--ink);
  font:16px/1.6 ui-sans-serif,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
}}
.wrap {{ max-width:760px; margin:0 auto; padding:28px 18px 56px; }}
header {{ display:flex; flex-direction:column; gap:6px; margin-bottom:22px; }}
.eyebrow {{
  font:600 11px/1.4 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
  letter-spacing:.14em; text-transform:uppercase; color:var(--muted);
}}
h1 {{ margin:0; font-size:26px; line-height:1.22; text-wrap:balance; letter-spacing:-.01em; }}
.sub {{ color:var(--muted); font-size:14.5px; }}
.live {{ display:inline-flex; align-items:center; gap:7px; color:var(--support); font-weight:600; }}
.dot {{ width:8px; height:8px; border-radius:50%; background:var(--support); }}
@media (prefers-reduced-motion:no-preference) {{
  .dot {{ animation:pulse 2.2s ease-in-out infinite; }}
  @keyframes pulse {{ 0%,100%{{opacity:1}} 50%{{opacity:.35}} }}
}}
.metrics {{ display:grid; grid-template-columns:repeat(3,1fr); gap:10px; margin:0 0 14px; }}
.metric {{
  background:var(--surface); border:1px solid var(--line); border-radius:10px; padding:12px 13px;
}}
.metric .k {{
  font:600 10.5px/1.3 ui-monospace,SFMono-Regular,Menlo,monospace;
  letter-spacing:.12em; text-transform:uppercase; color:var(--muted);
}}
.metric .v {{
  font:600 21px/1.25 ui-monospace,SFMono-Regular,Menlo,monospace;
  font-variant-numeric:tabular-nums; margin-top:3px;
}}
.track {{ height:7px; background:var(--bar); border-radius:4px; overflow:hidden; margin-bottom:26px; }}
.fill {{ height:100%; width:{i1/n*100:.1f}%; background:var(--accent); }}
figure {{ margin:0 0 26px; }}
figure img {{
  width:100%; height:auto; display:block; border:1px solid var(--line);
  border-radius:10px; background:var(--surface);
}}
figcaption {{ color:var(--muted); font-size:13.5px; margin-top:9px; }}
h2 {{ font-size:14px; letter-spacing:.02em; margin:0 0 10px; }}
table {{ width:100%; border-collapse:collapse; font-variant-numeric:tabular-nums; font-size:14.5px; }}
th, td {{ text-align:right; padding:8px 6px; border-bottom:1px solid var(--line); }}
th:first-child, td:first-child {{ text-align:left; }}
th {{
  font:600 10.5px/1.3 ui-monospace,SFMono-Regular,Menlo,monospace;
  letter-spacing:.1em; text-transform:uppercase; color:var(--muted);
}}
.note {{ color:var(--muted); font-size:13.5px; margin-top:14px; }}
.note code {{ font-size:12.5px; background:var(--bar); padding:1px 5px; border-radius:4px; }}
h2 {{ margin-top:30px; }}
table.kv td:first-child {{ width:38%; color:var(--muted); }}
table.kv td {{ text-align:left; }}
.q {{ color:var(--muted); font-style:italic; font-size:.88em; }}
.scroll {{ overflow-x:auto; }}
</style>
<div class="wrap">
  <header>
    <div class="eyebrow">Aldemir wall · replica · horizon 1.5</div>
    <h1>Reversed-cyclic analysis</h1>
    <div class="sub">
      {"<span class='live'><span class='dot'></span>running</span>" if running else "stopped"}
      · step {i1:,} of {n:,} · {sps*1e3:.0f} ms/step
    </div>
  </header>

  <div class="metrics">
    <div class="metric"><div class="k">Elapsed</div><div class="v">{t1/3600:.1f} h</div></div>
    <div class="metric"><div class="k">Remaining</div><div class="v">{rem/3600:.1f} h</div></div>
    <div class="metric"><div class="k">Finishes</div><div class="v">{eta:%a %H:%M}</div></div>
  </div>
  <div class="track"><div class="fill"></div></div>

  <figure>
    <img src="data:image/png;base64,{b64}" alt="Hysteresis loops and damage panels">
    <figcaption>Top: measured test hysteresis with our cyclic response overlaid, full range then
      zoomed. Bottom: crack pattern at four instants, sharing one strain scale so the panels are
      directly comparable. Pinch to zoom.</figcaption>
  </figure>

  <h2>Reversal peaks</h2>
  <div class="scroll">
    <table>
      <thead><tr><th>Level</th><th>mm</th><th>Push (kN)</th><th>Pull (kN)</th><th>Pull/Push</th></tr></thead>
      <tbody>{rows}</tbody>
    </table>
  </div>
  {"" if b64_rev is None else f'''
  <h2>Damage at each loop tip</h2>
  <figure>
    <img src="data:image/png;base64,{b64_rev}" alt="Crack pattern at each reversal">
    <figcaption>The crack field at every displacement reversal &mdash; the worst instant of each
      half-cycle. One shared strain scale, so the six panels are directly comparable. Tension-only
      classification: these struts have no compressive strength to exceed.</figcaption>
  </figure>'''}

  {"" if b64_evo is None else f'''
  <h2>How the damage accumulated</h2>
  <figure>
    <img src="data:image/png;base64,{b64_evo}" alt="Cracked fraction against drift and step">
    <figcaption>Cracked fraction against drift (upper &mdash; the loops trace back over themselves,
      dotted lines mark the reversals) and against step (lower). Damage is monotonic in time: it
      rises through reversals and on unloading, and never recovers.</figcaption>
  </figure>'''}

  <h2>Specimen and model</h2>
  <div class="scroll"><table class="kv"><tbody>{shared}</tbody></table></div>

  <h2>Runs in this figure</h2>
  <div class="scroll">
    <table>
      <thead><tr><th>Run</th><th>h</th><th>Concrete law</th><th>Drive</th>
        <th>Nodes</th><th>Elements</th><th>Result</th></tr></thead>
      <tbody>{ours}</tbody>
    </table>
  </div>
  <p class="note">All three: thickness {tw:.0f} mm, grid {sp.MESH:.0f} mm, explicit
    CentralDifference, {par.get('rate')} mm/s drive, {par.get('damping'):.0%} damping, longitudinal
    bars to the top. The cyclic run swaps the concrete law because
    <code>ElasticMultiLinear</code> is path-independent &mdash; a cracked strut would recover full
    stiffness on reload.</p>

  <h2>References</h2>
  <div class="scroll"><table class="kv"><tbody>{refs}</tbody></table></div>

  <p class="note">Currently {u[-1]:+.2f} mm at {s[-1]:+,.0f} kN, {max(abs(u))/HW:.4%} drift reached
    against the test's ~0.71%. The pull half-cycle running consistently stronger than the push is
    crack closure: each reversal drives material cracked on the previous half-cycle into
    compression, where closed cracks carry load.</p>
</div>
"""
    dst = run / "status.html"
    dst.write_text(html)
    print(f"saved {dst}  ({len(html)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
