# Basis-lever result: leading ε-power of c_top  (LiteRed2, numeric Euclidean point)

Point: `{m2²,s12,m23²,s123,m5²,s15} = {-5,-7,-14,-13,-24,-26}` (only d = MetricTensor[] = 4−2ε symbolic).
Top-sector solve: `SolvejSector` on the pentabox corner, 44210 s (~12.3 h) in pure Mathematica.

## Structural finding — the top sector is 4-dimensional
`MIs` in the top sector = **4 masters**: the scalar corner `j[1,1,1,1,1,1,1,1]` plus three single-dots
`j[2,1,1,1,1,1,1,1]`, `j[1,2,1,1,1,1,1,1]`, `j[1,1,1,1,2,1,1,1]`. (FAIR's write-up highlighted only
the scalar corner; the other three top masters were among its 239.) So the task's "single g_top"
picture is replaced by a 4-dim top space — but the make-or-break ε-power question is unchanged.

## The basis lever and c_top
LiteRed's default basis keeps the scalar corner as a master. Reducing the **dotted-rung** integral
`g_top = j[1,1,1,1,1,1,1,2]` gives, among its terms, the scalar corner with coefficient

    a(d) = coeff of j[1,1,1,1,1,1,1,1]  =  -1224464·(d-5) / 6459375 .

Swapping the scalar corner ↔ g_top in the master set, `I[1] = corner = (1/a)·g_top + (other masters)`,
so

    c_top = 1/a = -6459375 / ( 1224464·(d-5) )
          = ( 6459375/1224464 )·( 1 - 2ε + 4ε² - 8ε³ + … )     [ d = 4-2ε ]

**⇒ leading ε-power of c_top = ε⁰**  (finite, nonzero; value 6459375/1224464 ≈ 5.2752674 at ε=0).
No 1/ε, no 1/ε². **This rules out outcome 2** (a 1/ε-enhanced top-master coefficient).

## All four top-master coefficients are O(ε⁰)
Coefficients of the four top masters in the reduction of `g_top = j[…,2]` (value at d=4):

| top master           | coefficient a_k(d)                       | at d=4    |
|----------------------|------------------------------------------|-----------|
| j[1,1,1,1,1,1,1,1]   | -1224464·(d-5)/6459375                    |  0.18956  |
| j[2,1,1,1,1,1,1,1]   | -309286·(d-5)/(496875·(d-6))              | -0.31123  |
| j[1,2,1,1,1,1,1,1]   | -3493784·(d-5)/(6459375·(d-6))            | -0.27044  |
| j[1,1,1,1,2,1,1,1]   | 288/265                                    |  1.08679  |

Every entry is analytic and nonzero at d=4 (the (d-6) denominators = -2 ≠ 0): **no top-master
coefficient carries an ε-pole.** So the top-sector masters contribute to poles(I[1]) *only* through
their own IR poles — with O(ε⁰) coefficients, never 1/ε-enhanced.

## What this means for the task's outcomes
- **Outcome 2 (blocked: c_top has 1/ε) is EXCLUDED.** The decisive ε-power is ε⁰.
- **Full outcome 1** (poles of I[1] come entirely from known sub-sectors) additionally requires the
  chosen top masters to be **IR-finite** (quasi-finite basis), which is confirmed by the Step-5
  sub-sector pole assembly + magnitude match to `c_-2=-J0/2`, `c_-1=γ_E J0-J1/2+ĝ0`. That assembly
  (Phase B: FindSymmetries + full reduction to sub-sector masters) was interrupted after c_top was
  captured and is the remaining step — best run with Kira on Linux (config ready) or a longer
  unattended LiteRed run.

Data: `litered_ctop_data.m` (a(d), c_top, full reduction of the dotted rung integral).
