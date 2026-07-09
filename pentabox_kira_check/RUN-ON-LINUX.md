# Finishing the pentabox basis-lever pole check on Linux (Kira)

The decisive number is already in hand (`CTOP_RESULT.md`: `c_top` is **O(ε⁰)**, outcome 2 excluded).
What remains is **Phase B** — the full symmetrized reduction of `I[1]` to sub-sector masters and the
pole-magnitude match to `c_-2=-J0/2`, `c_-1=γ_E J0-J1/2+ĝ0`. Kira does this in minutes on Linux
x86_64 (its prebuilt static binaries work there; the macOS-arm64 source build on this box is broken —
see REPORT.md). Everything needed is in this directory.

## 1. Get Kira (Linux x86_64) — no build
    wget https://kira.hepforge.org/downloads?f=kira-2.0-static -O kira && chmod +x kira
    # (or the latest static release listed on https://kira.hepforge.org)
Fermat is bundled/optional for the static build; if asked, point FERMATPATH at any Fermat 7.x
(the arm64 `Ferm7a/fer64` here is Mac-only — on Linux use Ferl7 from home.bway.net/lewis/ferm7.html).

## 2. Files already prepared (copy this whole `kira_pentabox/` dir)
- `config/kinematics.yaml`      — SYMBOLIC kinematics (9 invariants); for the numeric point instead,
                                   use the numeric SP table in `config/kinematics_numeric.yaml.bak`.
- `config/integralfamilies.yaml`— 8 props (Eq.props) + 3 ISPs, top sector 255, k2=pentagon,k1=box.
- `jobs_natural.yaml`           — symmetries ON, reduce the corner (master-count vs FAIR 239).
- `jobs_lever.yaml`             — **the basis lever**: `preferred_masters: preferred` + reduce target.
- `preferred`                   — dotted-rung top master `pentabox[1,1,1,1,1,1,1,2,0,0,0]`.
- `target`                      — `pentabox[1,1,1,1,1,1,1,1,0,0,0]` (= I[1]).

NOTE the top sector is **4-dimensional** (LiteRed found masters: corner + dots on lines 1,2,5). To
push *all four* off the scalar corner, list four dotted/numerator top integrals in `preferred`, e.g.
    pentabox[1,1,1,1,1,1,1,2,0,0,0]
    pentabox[2,1,1,1,1,1,1,1,0,0,0]
    pentabox[1,2,1,1,1,1,1,1,0,0,0]
    pentabox[1,1,1,1,2,1,1,1,0,0,0]
(or numerators via negative ISP indices) so `I[1]` reduces onto a quasi-finite top basis + subsectors.

## 3. Run
    FERMATPATH=/path/to/fermat ./kira jobs_lever.yaml
    # results/pentabox/masters              -> master list (per-sector props); compare count to FAIR 239
    # results/pentabox/kira_target.m        -> I[1] = Σ c_j(ε,{s}) M_j   (Mathematica form)

## 4. Pole assembly + comparison (Mathematica)
From `results/pentabox/kira_target.m`:
1. Identify each master M_j with a 1/ε² or 1/ε coefficient; these are ≤7-prop sub-sectors
   (degenerate boxes/triangles/bubbles, factorizable 1-loop×1-loop). The 3-mass box = Ellis–Zanderighi
   (0712.1851) Box 5, already transcribed in `pentabox.tex`.
2. Substitute their Laurent expansions, collect `[1/ε²]I[1]=P2`, `[1/ε]I[1]=P1`.
3. Evaluate at `{m2²,s12,m23²,s123,m5²,s15}={-5,-7,-14,-13,-24,-26}`; compare **magnitudes** to
   `|c_-2|=1.914807547604452136e-5`, `|c_-1|=2.197205645392600987e-4`; **report the sign** (a flip is
   the known bookkeeping issue, not a failure).
4. If P2,P1 match → the four top masters are quasi-finite by consistency → **outcome 1** (analytic
   pole check closes). If a genuinely-2-loop non-factorizable sub-master with a pole is untabulated
   for our kinematics → name it (**outcome 3**).

Since `c_top` (and all top-master coefficients) are already known to be **O(ε⁰)**, a magnitude match
in step 3 closes the analytic check without ever evaluating the untabulated top masters.
