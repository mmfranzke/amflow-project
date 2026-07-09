# TASK: finish the analytic pole check of the two-loop pentabox I[1] via Kira on Linux

You are an autonomous coding agent running on the ETH Euler HPC cluster (Linux x86_64). Complete the
sub-sector pole assembly ("Phase B") for the pentabox counterterm cross-check, install any software
you need, and report the outcome. Work in THIS directory (a subfolder of the git repo) and write all
outputs here so they can be committed.

## Where this lives (repository context)
This folder is `pentabox_kira_check/` inside the git repository **amflow-project**
(github.com/mmfranzke/amflow-project). That repo already contains the *numeric* half of this study —
a **working AMFlow setup** that evaluated `I[1]` directly and matched the counterterm poles to 12
digits. Two consequences:
- **Reuse the repo's AMFlow** for Step 4's numeric fallback: the pole-carrying sub-sector masters are
  ≤7-line integrals; evaluating them at the reference point with the existing AMFlow (see the repo's
  `amflow_pentabox.wl` / `amflow_logs/` / `cache/`) is the reliable, already-installed route — no need
  to reinstall AMFlow.
- When done, **commit your results** into this folder (`git add pentabox_kira_check && git commit`),
  do NOT touch other parts of the repo. Leave the working tree clean otherwise.

## Background (what is already established — do NOT redo)
An independent analytic cross-check of the poles of the bare two-loop pentabox
`I[1] = Pentabox(1,…,1)` is wanted. A basis-lever IBP reduction (LiteRed, on a Mac) already produced
the make-or-break number: with the scalar corner swapped out of the master basis for a dotted top
integral, the top-sector coefficient `c_top = 1/a = -6459375/(1224464 (d-5))` is **O(ε⁰)** (finite,
no 1/ε), and the pentabox **top sector is 4-dimensional** (masters: scalar corner + single dots on
propagators 1, 2, 5). All four top-master coefficients are O(ε⁰). This EXCLUDES the "blocked"
outcome. What remains is to confirm the poles of `I[1]` come entirely from the (known) sub-sector
masters — i.e. the full symmetrized reduction + pole assembly + magnitude match. That is your job.

## The integral and kinematics
```
I[1] = ∫∫ 1/(A1…A8),  ∫ ≡ ∫ d^dℓ/(iπ^{d/2}),  d = 4-2ε
A1=ℓ2²  A2=(ℓ2+p1)²  A3=(ℓ2+p1+p2)²  A4=(ℓ2+p1+p2+p3)²
A5=ℓ1²  A6=(ℓ1+p5)²  A7=(ℓ1+p5+p4)²  A8=(ℓ1+ℓ2)²
p1²=0 only; p2,p3,p4,p5 off-shell; Σp=0.   (k2≡ℓ2 pentagon, k1≡ℓ1 box)
```
Six invariants entering the poles: `m2²,s12,m23²,s123,m5²,s15` (plus pole-irrelevant `m3²,m4²`).
Numeric Euclidean reference point `{m2²,s12,m23²,s123,m5²,s15} = {-5,-7,-14,-13,-24,-26}`
(with derived m3²=-45, m4²=-3, s34=-33). Reference pole coefficients to match (MAGNITUDE; a global
sign flip is a KNOWN bookkeeping issue — report the sign, do NOT force agreement):
```
c_-2 = -J0/2                   |c_-2| = 1.914807547604452136e-5
c_-1 = γ_E J0 - J1/2 + ĝ0      |c_-1| = 2.197205645392600987e-4
```

## Step 1 — install Kira (Linux x86_64)
Try, in order, whatever works on this cluster (use `module avail`/`module load` if present):
1. Prebuilt static binary: download from https://kira.hepforge.org (kira / kira-2.x static). `chmod +x`.
2. conda/mamba: `mamba create -n kira -c conda-forge kira firefly` (if the feedstock is available on
   linux-64 — it is not on osx-arm64 but may be on linux-64).
3. Build from source (meson): deps GiNaC+CLN, yaml-cpp, zlib, FLINT, Fermat at runtime. FireFly is a
   meson subproject. **Build WITH FLINT** (`meson setup -Dflint=true …`) — the finite-field
   modular arithmetic needs it. Get Fermat 7.x for Linux from https://home.bway.net/lewis/ferm7.html
   (`Ferl7.tar.gz`), `chmod +x`, and set `FERMATPATH=/path/to/fer64`.
Verify the install by reducing Kira's shipped `examples/1-loop-box` and confirming it yields a
NONZERO master list (a "0 masters" result means the finite-field backend is broken — fix FLINT/Fermat
before proceeding). Record tool + version.

## Step 2 — reduce I[1] with symmetries ON and the basis lever
The Kira family is already written (`config/`, `jobs_lever.yaml`, `preferred`, `target`):
- `config/kinematics.yaml` keeps the 9 invariants SYMBOLIC (preferred: gives P2,P1 as functions).
  If full symbolic reconstruction is too slow, switch to the numeric point in
  `config/kinematics_numeric.yaml.bak` (only d symbolic).
- `preferred` lists FOUR dotted top masters (not the scalar corner) → forces `I[1]` to reduce onto
  a quasi-finite top basis + sub-sectors.  `jobs_lever.yaml` has `run_symmetries: true`.
Run: `FERMATPATH=… ./kira jobs_lever.yaml`.  Outputs: `results/pentabox/masters` (master list; report
the count and per-sector propagator content, compare to FAIR's 239 masters found WITHOUT symmetries),
and `results/pentabox/kira_target.m`  ( `I[1] = Σ c_j(ε,{s}) M_j` in Mathematica form).

## Step 3 — classify the masters and their coefficients
Laurent-expand every coefficient `c_j` in ε. Split the masters into: (a) the four TOP-sector masters
(8 lines) and (b) SUB-sector masters (≤7 lines: degenerate boxes, triangles, bubbles, factorizable
one-loop×one-loop products). Confirm the four top-master coefficients are O(ε⁰) (consistency with the
LiteRed result). List which SUB-sector coefficients carry 1/ε² or 1/ε.

## Step 4 — evaluate the pole-carrying sub-sector masters and assemble
For each SUB-sector master with a 1/ε² or 1/ε coefficient, obtain its Laurent expansion:
- Prefer literature: Ellis–Zanderighi (arXiv:0712.1851) — the 3-mass box is their Box 5 (already
  transcribed in `pentabox.tex`); standard one-/two-loop bubbles, triangles, boxes; factorizable
  two-loop masters = products of one-loop results. **Cite each.**
- If a master's closed pole part is not readily available, evaluate it numerically at the reference
  point with `pySecDec` (install via pip) or AMFlow — these are ≤7-line integrals and are cheap.
- If a genuinely-two-loop NON-factorizable sub-master with a pole is untabulated for this
  one-lightlike/four-off-shell kinematics: name it and stop (**outcome 3**).
Assemble `[1/ε²]I[1] = P2` and `[1/ε]I[1] = P1` from the sub-sector contributions (the four top-master
terms contribute nothing to the poles iff those masters are IR-finite — which the match below tests).

## Step 5 — compare and classify the outcome
Evaluate `P2, P1` at `{m2²,s12,m23²,s123,m5²,s15}={-5,-7,-14,-13,-24,-26}` to ≥10 digits. Compare
MAGNITUDES to `|c_-2|`, `|c_-1|` above; report the sign. If `P2` is compact, also attempt a symbolic
match `P2 = -J0/2`. Report exactly one outcome:
1. **Closes:** sub-sector poles reproduce `c_-2,c_-1` (magnitude; sign noted) ⇒ the four top masters
   are IR-finite by consistency ⇒ analytic pole check closes. (Strongest result.)
2. **Top master carries a pole after all:** if the sub-sector poles do NOT reproduce `c_-2,c_-1`
   (i.e. a top master must supply a pole) ⇒ route needs the untabulated top master. Report and stop.
3. **A sub-sector two-loop master is unknown for our kinematics:** name it; partial result.

## Deliverables
`REPORT_CLUSTER.md` with: tool+version; master count (symmetries ON) vs FAIR's 239; the explicit
`I[1] = Σ c_j M_j`; the four top-master coefficients (confirm O(ε⁰)); the pole-carrying sub-sector
masters + their references; `P2,P1` vs `c_-2,c_-1` (magnitude + sign); and the outcome (1/2/3). Keep
the Kira config/run files and any evaluation scripts. Do NOT edit `pentabox.tex` unless outcome 1.

## Rules
- Symmetries ON; symbolic kinematics if feasible (numeric point otherwise). Compare pole magnitudes,
  report signs, do NOT force agreement. Cap any single reduction at ~30 min wall and report progress;
  Kira+FireFly (finite-field) should not hang. Verify the install produces nonzero masters on the
  shipped example BEFORE trusting the pentabox run.
```
```
FILES PRESENT IN THIS DIRECTORY (uploaded):
  config/kinematics.yaml, config/integralfamilies.yaml, config/kinematics_numeric.yaml.bak
  jobs_lever.yaml, jobs_natural.yaml, preferred, target
  CTOP_RESULT.md         (the established c_top=O(ε⁰) result + 4-master finding)
  pentabox.tex           (propagator defs, invariants, Ellis–Zanderighi Box 5, reference values)
```
