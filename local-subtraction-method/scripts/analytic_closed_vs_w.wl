(* Pure analytic check of the pre-collapsed w representation against the
   compact closed form.  This script uses only exact rational arithmetic.

   Equation-label map to derivations/snippets/dckernel/main.tex:
     - support point and theta region: eq:support-theta
     - lambda: eq:lambda
     - a, b: eq:method1-a-b-def
     - mu_l^2 and mu_k^2 definitions: eq:mu-l-def, eq:lambda
     - Method-1 w representation: eq:method1-w-integral-before-euler
     - Delta_x(w): eq:method1-Delta-x-def
     - Delta0, DeltaW: eq:method1-Delta0-Deltaw-def
     - Euler variable and z: eq:method1-euler-map, eq:method1-theta-z-def
     - hypergeometric representation: eq:method1-hypergeometric-form
     - Delta1 from DeltaW: eq:method1-Delta1-from-Deltaw
     - compact kernel: eq:closed-form-kernel
     - final Delta0, Delta1: eq:Delta0-final, eq:Delta1-final *)

ClearAll["Global`*"];

x = 3/10;
y = 1/4;
s = 0; (* s is p^2 in eq:method1-Delta0-Deltaw-def and eq:Delta1-final. *)
ml2 = 49/100;
Ml2 = 4;
mk2 = 1/4;
Mk2 = 81/100;

(* eq:lambda, with X = 1 - x from the support variables. *)
X = 1 - x;
lambda = y/X;

(* eq:method1-a-b-def *)
a = lambda (1 - lambda);
b = lambda (1 - lambda x);

(* eq:mu-l-def and eq:lambda for mu_k^2(lambda). *)
muLx = (1 - x) ml2 + x Ml2;
muKlambda = (1 - lambda) mk2 + lambda Mk2;

(* eq:method1-Delta0-Deltaw-def, with p^2 -> s. *)
Delta0 = muLx - x (1 - x) s;
DeltaW = (
  ((1 + a) x (1 - x) + b) s
  - muKlambda
  + a Ml2
  - (1 + a) muLx
);

(* eq:method1-Delta1-from-Deltaw *)
Delta1 = (1 + a) Delta0 + DeltaW;

(* eq:Delta1-final, used as a cross-check of eq:method1-Delta1-from-Deltaw. *)
Delta1Alt = (
  -muKlambda
  + lambda (1 - lambda) Ml2
  + lambda (1 - lambda x) s
);

(* eq:method1-w-integral-before-euler and eq:method1-Delta-x-def *)
U[w_] := 1 - (1 + a) w;
DeltaX[w_] := Delta0 + w DeltaW;

(* eq:method1-theta-z-def *)
z = -DeltaW/((1 + a) Delta0);

(* Euler/hypergeometric form of the w integral.  The direct Integrate form is:
   Gamma[2 eps]/X Integrate[
     w^(eps - 1) U[w]^(-1 + eps) DeltaX[w]^(-2 eps),
     {w, 0, 1/(1 + a)}, Assumptions -> eps > 0]
   We use the Euler form to keep endpoint factors explicit and avoid a slow
   general-purpose integration step. *)
(* eq:method1-hypergeometric-form, before applying
   Hypergeometric2F1[A, B, A, z] = (1 - z)^(-B). *)
IwHyper[eps_] := Gamma[2 eps]/X *
  (1 + a)^(-eps) *
  Delta0^(-2 eps) *
  Beta[eps, eps] *
  Hypergeometric2F1[2 eps, eps, 2 eps, z];

(* eq:method1-hypergeometric-form after the hypergeometric identity and before
   replacing (1 + a) Delta0 (1 - z) by Delta1. *)
IwCollapsed[eps_] := Gamma[2 eps]/X *
  (1 + a)^(-eps) *
  Delta0^(-2 eps) *
  Beta[eps, eps] *
  (1 - z)^(-eps);

(* eq:closed-form-kernel with eq:Delta0-final and eq:Delta1-final. *)
Iclosed[eps_] := Gamma[eps]^2/X * Delta0^(-eps) * Delta1^(-eps);

coefficients[expr_] := FullSimplify[
  Table[
    SeriesCoefficient[FunctionExpand[expr], {eps, 0, power}],
    {power, -2, 0}
  ],
  eps > 0
];

printValue[label_, value_] := Print[label, " = ", InputForm[value]];

Print["== Kinematic point =="];
printValue["x (eq:support-theta)", x];
printValue["y (eq:support-theta)", y];
printValue["s = p^2", s];
printValue["ml2", ml2];
printValue["Ml2", Ml2];
printValue["mk2", mk2];
printValue["Mk2", Mk2];

Print["\n== Derived scales =="];
printValue["X", X];
printValue["lambda (eq:lambda)", lambda];
printValue["a (eq:method1-a-b-def)", a];
printValue["b (eq:method1-a-b-def)", b];
printValue["muLx = mu_l^2(x) (eq:mu-l-def)", muLx];
printValue["muKlambda = mu_k^2(lambda) (eq:lambda)", muKlambda];
printValue["Delta0 (eq:method1-Delta0-Deltaw-def)", Delta0];
printValue["DeltaW (eq:method1-Delta0-Deltaw-def)", DeltaW];
printValue["Delta1 (eq:method1-Delta1-from-Deltaw)", FullSimplify[Delta1]];
printValue["Delta1Alt (eq:Delta1-final)", FullSimplify[Delta1Alt]];
printValue[
  "Delta1 - Delta1Alt (eq:method1-Delta1-from-Deltaw vs eq:Delta1-final)",
  FullSimplify[Delta1 - Delta1Alt]
];
printValue["z (eq:method1-theta-z-def)", FullSimplify[z]];

Print["\n== Expected sanity checks =="];
printValue["lambda == 5/14", FullSimplify[lambda == 5/14]];
printValue["Delta0 == 1543/1000", FullSimplify[Delta0 == 1543/1000]];
printValue["Delta1 == 459/980", FullSimplify[Delta1 == 459/980]];

Print["\n== Hypergeometric identity check =="];
printValue[
  "FunctionExpand[IwHyper - IwCollapsed] (eq:method1-hypergeometric-form)",
  FullSimplify[FunctionExpand[IwHyper[eps] - IwCollapsed[eps]], eps > 0]
];
printValue[
  "FunctionExpand[IwCollapsed - Iclosed] (eq:method1-theta-z-def + eq:method1-Delta1-from-Deltaw -> eq:closed-form-kernel)",
  FullSimplify[FunctionExpand[IwCollapsed[eps] - Iclosed[eps]], eps > 0]
];
printValue[
  "FunctionExpand[IwHyper - Iclosed] (eq:method1-hypergeometric-form -> eq:closed-form-kernel)",
  FullSimplify[FunctionExpand[IwHyper[eps] - Iclosed[eps]], eps > 0]
];

Print["\n== Laurent coefficients {eps^-2, eps^-1, eps^0} =="];
printValue["IwHyper (eq:method1-hypergeometric-form)", coefficients[IwHyper[eps]]];
printValue["IwCollapsed (hypergeometric identity applied)", coefficients[IwCollapsed[eps]]];
printValue["Iclosed (eq:closed-form-kernel)", coefficients[Iclosed[eps]]];
printValue["IwHyper - Iclosed", coefficients[IwHyper[eps] - Iclosed[eps]]];
printValue[
  "IwCollapsed - Iclosed",
  coefficients[IwCollapsed[eps] - Iclosed[eps]]
];
printValue[
  "IwHyper - IwCollapsed",
  coefficients[IwHyper[eps] - IwCollapsed[eps]]
];
printValue["Iclosed/2 diagnostic", coefficients[Iclosed[eps]/2]];

Print["\n== Differences =="];
printValue[
  "IwHyper - Iclosed (eq:method1-hypergeometric-form -> eq:closed-form-kernel)",
  FullSimplify[FunctionExpand[IwHyper[eps] - Iclosed[eps]], eps > 0]
];
printValue[
  "IwCollapsed - Iclosed (eq:method1-theta-z-def + eq:method1-Delta1-from-Deltaw -> eq:closed-form-kernel)",
  FullSimplify[FunctionExpand[IwCollapsed[eps] - Iclosed[eps]], eps > 0]
];
printValue[
  "IwHyper - IwCollapsed (hypergeometric identity)",
  FullSimplify[FunctionExpand[IwHyper[eps] - IwCollapsed[eps]], eps > 0]
];

Print["\n== Cheap numerical check at epsTest = 1/1000 =="];
epsTest = 1/1000;
printValue["N[IwHyper[epsTest], 50]", N[IwHyper[epsTest], 50]];
printValue["N[IwCollapsed[epsTest], 50]", N[IwCollapsed[epsTest], 50]];
printValue["N[Iclosed[epsTest], 50]", N[Iclosed[epsTest], 50]];
printValue[
  "N[IwHyper[epsTest] - Iclosed[epsTest], 50]",
  N[IwHyper[epsTest] - Iclosed[epsTest], 50]
];
printValue[
  "N[IwCollapsed[epsTest] - Iclosed[epsTest], 50]",
  N[IwCollapsed[epsTest] - Iclosed[epsTest], 50]
];
