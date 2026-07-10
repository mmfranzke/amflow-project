# TASK (cluster, Euler): the decisive sub-sector pole assembly for the pentabox I[1]

Autonomous agent on the ETH Euler cluster. Everything is Python (no Mathematica). Work in this dir.

## Goal
Decide whether the poles of the two-loop pentabox `I[1]` are isolated to its sub-sector masters, by
assembling the sub-sector contribution and comparing to the known double pole `c_-2 = -J0/2 =
-1.914807547604452136e-5`.

Background (established): the reduction `I[1] = Sum_j c_j(d) M_j` (in `kira_target.m`) has 4 top-sector
masters with O(eps^0) coefficients and 51 genuinely-new "class-D" sub-sector masters that are PROVEN
FINITE. The open question: do the 4 (untabulated) TOP masters contribute to the poles? They are dots
on jet lines, so likely divergent -> likely DO contribute. Test it directly:

  P2_sub := [1/eps^2] ( Sum_{j != top} c_j M_j )   vs   c_-2.
  * P2_sub ~ c_-2  -> top masters net zero -> poles isolated to sub-sectors (route closes).
  * P2_sub = O(1)  -> top masters contribute -> route blocked by them; the mismatch measures it.

## Steps
1. `pip install --user pySecDec sympy` (Euler has Python; use a venv or --user; needs a C++ compiler,
   `module load gcc` if required). Verify `python -c "import pySecDec, sympy"`.
2. Run the master batch:  `python eval_sub_masters.py`  (writes `master_laurents.json`, checkpointed;
   evaluates the 161 sub-sector masters via pySecDec QMC to ~15 digits; skips the 4 top masters and
   the 23 exact bubble-chains). This is the long step (~hours); it is restartable (skips done masters).
   Submit as an sbatch job or run in `tmux`/`nohup` on a login node with internet for the pip install.
3. When it prints `DONE SUB MASTERS`, assemble:  `python assemble_sub.py`  -> prints P2_sub, P1_sub,
   the `[1/eps^3]` cancellation check (expect ~0), and `P2_sub - c_-2`.
4. Report `REPORT_POLE2.md`: the P2_sub value, `[1/eps^3]`, whether P2_sub matches c_-2 (magnitude +
   sign), and the verdict (isolated / top masters contribute). If any masters FAILED in pySecDec,
   list them (7-line masters are the hard ones); a few failures may need higher `maxeval` or AMFlow.

## Files here
- `eval_sub_masters.py`  — pySecDec QMC batch (family + kinematics baked in; korobov3, minn=1e6).
- `assemble_sub.py`      — pure-Python assembler (sympy exact coeffs + Gamma bubble-chains + pySecDec).
- `kira_target.m`        — the reduction `I[1] = Sum c_j M_j` (Mathematica syntax; parsed by sympy).

## Notes
- Precision: the assembly has an ~8-order cancellation, so masters need ~14-15 digits (QMC minn=1e6
  delivers this; validated: pySecDec reproduced the sector-53 bubble-chain and sector-150 triangle to
  14-16 digits vs exact). Do NOT lower the precision.
- No new IBP runs; reuse `kira_target.m`.
