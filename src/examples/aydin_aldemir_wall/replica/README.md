# Aydin (2019) replica — Aldemir wall

His model, not ours with his numbers. Built to separate two questions that six failed experiments
in the parent study kept conflating: does our model collapse because of **how we built it**, or
because of something the lattice method does to **this specimen**?

| | parent study | this replica | source |
|---|---|---|---|
| panel | 3000 x 2250 | **3000 x 2680** | his Table 2 counts, solved exactly — but see the caution below |
| mesh | 50 mm | **20 mm** | his Tables 1 and 2 |
| nu | 0.20 | **1/3** | his Appendix (Hrennikoff) |
| calibration | uniaxial (thesis) | **equibiaxial (his Appendix)** | C = 0.621 |
| concrete | Concrete02, bilinear tail | **tension-only, trilinear** | his Fig. 1(c) |
| tension tail | solved from Gf | **a2 = 70, a3 = 360 as published** | his Table 1 |
| steel | Steel02, b = 0.01 | **elastic-perfectly plastic, b = 0** | his Fig. 1(d) |

**CAUTION — the panel may be on its side (D82).** The Table 2 inversion is unique only up to
TRANSPOSITION: 3000 x 2680 and 2680 x 3000 both give his 20,385 nodes and 80,684 struts exactly.
This replica takes the first. Fig. 4(f) dimensions 2680 along the HORIZONTAL, and that figure is a
detail view of the bottom 1500 mm of his model rather than a panel, so its width label would be the
model's true full width — which points at the second. If that is right, every ratio below moves.
Unresolved; it needs the author.

**Targets** (his Table 4, horizon 1.5d): K_sim = 943.16 kN/mm, F_sim = 1164.413 kN.
Measured: K = 1038.44 kN/mm, F = 963.592 kN.

## Bond (`--bond`)

His steel nodes are separate from the concrete and tied to it by elastic-brittle links (a = 0.7).
Shared nodes were the replica's largest known departure; `--bond` removes it, using the same
force-slip law the parent study needed (D75) — a truss couples `k = EA/L` and `F = ft*A`, so the
link area alone cannot set stiffness and strength, and sizing it for stiffness alone collapsed the
parent wall *earlier* than no bond at all.

Cost: **20,385 nodes / 88,324 elements → 28,081 / 149,802** (61,478 links, a ring of 8 per steel
node). His Table 2 counts are unchanged, because they count concrete struts only.

**The area ratio is inferred, not his.** He publishes `a = 0.7` and nothing else — not the link
area, peak slip or peak force. The default 0.01 was calibrated against perfect bond at the parent
study's **50 mm** mesh; the replica meshes at **20 mm**, and the failure slip scales with link
length, so the same ratio gives 0.149 mm here against 0.372 there. Both sit inside the Model Code
`s1` band for a deformed bar (0.1–1 mm), but this one sits on its floor, and 0.003 went singular at
mesh 50 — so the usable window is narrow and is not known at this mesh. Check it before trusting a
pushover, and read `K/K_perfect` as a ratio between the two elastic runs:

```
run.py --elastic            # perfect bond, the denominator
run.py --elastic --bond     # bond can only ADD flexibility; > 1 is the D72 artefact, and sizes it
```

## Still not replicated

`specimen.not_replicated(bond=…, explicit=…)`, printed by every run with only the entries that
still apply. With `--bond --explicit` what remains is his unpublished bond parameters, the
horizon-ring topology (which cannot reduce to perfect bond in the limit, D72), his PID force drive
against our prescribed displacement, `b1`/`b2`, and the thickness.

## Use

```
uv run python examples/aydin_aldemir_wall/replica/run.py --elastic              # cheap check first
uv run python examples/aydin_aldemir_wall/replica/preflight.py --bond --explicit  # price it
uv run python examples/aydin_aldemir_wall/replica/run.py --drift 0.0025 --bond --explicit
```

`preflight.py` takes the same flags as `run.py` and must be given them: an implicit probe does not
price an explicit run, and the no-bond model does not price the bond one.

## Where it stood before bond

Without bond the replica reached **0.42 of his published peak** and collapsed at 0.098% drift.
In the parent study bond moved peak ×1.14 and collapse drift ×1.21 — a real gain that did not
remove the failure. Budget past the no-bond collapse drift, not to it.
