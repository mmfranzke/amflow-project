(* Pure Mathematica diagnostic for the transformation from Eq. (29) to
   Eq. (43) in derivations/snippets/dckernel/main.tex.

   This script does not use AMFlow.  It only checks scalar-product algebra,
   Feynman-parameter prefactors, the shift, and the delta-collapse Jacobian. *)

ClearAll["Global`*"];

printValue[label_, value_] := Print[label, " = ", InputForm[value]];

Print["== Eq. (29) to Eq. (43): symbolic denominator check =="];

(* Scalar-product representation:
   l2 = l.l, lp = l.p, p2 = s. *)
A = l2 - ml2;
B = l2 - 2 lp + s - Ml2;
Cscale = muKlambda - a (s - 2 lp + l2) - b s;

Ddirect = u A + v B + w Cscale;

U = u + v - a w;
V = v - a w;
W = (
  (v - w (a + b)) s
  - u ml2
  - v Ml2
  + w muKlambda
);

Dquad = U l2 - 2 V lp + W;

denominatorDifference = FullSimplify[Ddirect - Dquad];

printValue["A", A];
printValue["B", B];
printValue["C[l]", Cscale];
printValue["Ddirect = u A + v B + w C[l]", Ddirect];
printValue["Dquad = U l2 - 2 V lp + W", Dquad];
printValue["FullSimplify[Ddirect - Dquad]", denominatorDifference];

If[denominatorDifference =!= 0,
  Print["WARNING: Ddirect - Dquad did not simplify to zero."];
];

Print["\n== Shifted quadratic form check =="];

Delta = W - V^2/U s;
Dshift = U q2 + Delta;

shiftRules = {
  l2 -> q2 + 2 (V/U) qp + (V/U)^2 s,
  lp -> qp + (V/U) s
};

shiftDifference = FullSimplify[
  (Dquad /. shiftRules) - Dshift,
  U != 0
];

printValue["Delta = W - V^2/U s", Delta];
printValue["Dshift = U q2 + Delta", Dshift];
printValue["FullSimplify[(Dquad /. shiftRules) - Dshift]", shiftDifference];

If[shiftDifference =!= 0,
  Print["WARNING: shifted form did not simplify to zero."];
];

Print["\n== Feynman-parameter prefactor check =="];

Print[
  "General formula: 1/(A^1 B^1 C^eps) = ",
  "Gamma[2+eps]/(Gamma[1] Gamma[1] Gamma[eps]) ",
  "Integral[du dv dw delta(1-u-v-w) w^(eps-1)/D^(2+eps)]."
];
prefactorAfterOutsideGamma = FullSimplify[
  Gamma[eps] * Gamma[2 + eps]/(Gamma[1] Gamma[1] Gamma[eps])
];
printValue[
  "Outside Gamma[eps] times Feynman coefficient",
  prefactorAfterOutsideGamma
];
printValue[
  "Gamma[eps] Gamma[2]/Gamma[2+eps]",
  FullSimplify[Gamma[eps] Gamma[2]/Gamma[2 + eps]]
];

Print["\n== Gaussian/Fourier step check =="];

Print[
  "After l = q + (V/U) p, delta[x - l.eta] = ",
  "delta[x - q.eta - V/U]."
];
Print[
  "The Fourier source is proportional to eta, and eta^2 = 0, ",
  "so the q-Gaussian is independent of omega and returns delta[x - V/U]."
];

d = 4 - 2 eps;
gaussianFactor = U^(-d/2) Gamma[2 eps]/Gamma[2 + eps] Delta^(-2 eps);
eq43KernelFactor = FullSimplify[Gamma[2 + eps] gaussianFactor];

printValue[
  "Gaussian factor before multiplying Feynman prefactor",
  gaussianFactor
];
printValue[
  "Gamma[2+eps] times Gaussian factor",
  eq43KernelFactor
];

Print["\n== Delta-collapse check from Eq. (43) to Eq. (44) =="];

Ucollapsed = FullSimplify[U /. u -> 1 - v - w];
Vcollapsed = FullSimplify[V /. u -> 1 - v - w];

vStar = a w + x Ucollapsed;
uStar = FullSimplify[(1 - vStar - w)];

deltaArgument = x - Vcollapsed/Ucollapsed;
deltaDerivativeCheck = FullSimplify[
  D[deltaArgument, v] + 1/Ucollapsed,
  Ucollapsed != 0
];
vStarCheck = FullSimplify[
  deltaArgument /. v -> vStar,
  Ucollapsed != 0
];
uStarCheck = FullSimplify[
  uStar - (1 - x) Ucollapsed,
  Ucollapsed != 0
];

uPowerCheck = FullSimplify[
  (Ucollapsed^(-d/2) Ucollapsed)/(Ucollapsed^(-1 + eps)) /. d -> 4 - 2 eps,
  Ucollapsed > 0
];

printValue["U /. u -> 1-v-w", Ucollapsed];
printValue["V /. u -> 1-v-w", Vcollapsed];
printValue["vStar = a w + x U", vStar];
printValue["uStar", uStar];
printValue["FullSimplify[(x - V/U) /. v -> vStar]", vStarCheck];
printValue["FullSimplify[uStar - (1-x) U]", uStarCheck];
printValue["FullSimplify[D[x - V/U, v] + 1/U]", deltaDerivativeCheck];
printValue[
  "FullSimplify[U^(-d/2) * U / U^(-1+eps)] with d=4-2 eps",
  uPowerCheck
];

If[deltaDerivativeCheck =!= 0,
  Print["WARNING: delta-collapse Jacobian check did not simplify to zero."];
];
If[uPowerCheck =!= 1,
  Print["WARNING: U power check did not simplify to one."];
];

Print["\n== Branch-safe numerical sanity point =="];

xVal = 3/10;
yVal = 1/4;
sVal = 0;
ml2Val = 49/100;
Ml2Val = 4;
mk2Val = 1/4;
Mk2Val = 81/100;

XVal = 1 - xVal;
lambdaVal = yVal/XVal;
aVal = lambdaVal (1 - lambdaVal);
bVal = lambdaVal (1 - lambdaVal xVal);
muLxVal = (1 - xVal) ml2Val + xVal Ml2Val;
muKlambdaVal = (1 - lambdaVal) mk2Val + lambdaVal Mk2Val;
Delta0Val = muLxVal - xVal (1 - xVal) sVal;
DeltaWVal = (
  ((1 + aVal) xVal (1 - xVal) + bVal) sVal
  - muKlambdaVal
  + aVal Ml2Val
  - (1 + aVal) muLxVal
);
Delta1Val = (1 + aVal) Delta0Val + DeltaWVal;

printValue["lambdaVal", lambdaVal];
printValue["aVal", aVal];
printValue["bVal", bVal];
printValue["muKlambdaVal", muKlambdaVal];
printValue["Delta0Val", Delta0Val];
printValue["DeltaWVal", DeltaWVal];
printValue["Delta1Val", Delta1Val];

printValue["lambdaVal == 5/14", FullSimplify[lambdaVal == 5/14]];
printValue["aVal == 45/196", FullSimplify[aVal == 45/196]];
printValue["bVal == 125/392", FullSimplify[bVal == 125/392]];
printValue["muKlambdaVal == 9/20", FullSimplify[muKlambdaVal == 9/20]];
printValue["Delta0Val == 1543/1000", FullSimplify[Delta0Val == 1543/1000]];
printValue["DeltaWVal == -40009/28000", FullSimplify[DeltaWVal == -40009/28000]];
printValue["Delta1Val == 459/980", FullSimplify[Delta1Val == 459/980]];

