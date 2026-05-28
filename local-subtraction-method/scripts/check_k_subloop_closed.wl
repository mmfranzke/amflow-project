(* Pure analytic Laurent check for the closed k-subloop kernel K(q;y).
   This script does not use AMFlow. *)

ClearAll["Global`*"];

Xval = 7/10;
yval = 1/4;
q2val = 1;
mk2val = 1/4;
Mk2val = 81/100;

lambdaVal = yval/Xval;
muKlambdaVal = (1 - lambdaVal) mk2val + lambdaVal Mk2val;
CKval = muKlambdaVal - lambdaVal (1 - lambdaVal) q2val;

Kclosed[eps_] := Gamma[eps]/Xval * CKval^(-eps);

coefficients[expr_] := FullSimplify[
  Table[
    SeriesCoefficient[expr, {eps, 0, power}],
    {power, -1, 1}
  ],
  eps > 0
];

printValue[label_, value_] := Print[label, " = ", InputForm[value]];

Print["== k-subloop closed-form check =="];
printValue["Xval", Xval];
printValue["yval", yval];
printValue["q2val", q2val];
printValue["mk2val", mk2val];
printValue["Mk2val", Mk2val];

Print["\n== Derived scales =="];
printValue["lambdaVal", lambdaVal];
printValue["muKlambdaVal", muKlambdaVal];
printValue["CKval", CKval];
printValue["Sign[CKval]", Sign[CKval]];
printValue["CKval == 54/245", FullSimplify[CKval == 54/245]];

Print["\n== Laurent coefficients {eps^-1, eps^0, eps^1} =="];
printValue["Kclosed", coefficients[Kclosed[eps]]];
printValue[
  "leading coefficient == 10/7",
  FullSimplify[SeriesCoefficient[Kclosed[eps], {eps, 0, -1}] == 10/7]
];

