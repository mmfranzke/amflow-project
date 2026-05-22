(* targets/CompareOneLoop.wl
   Fresh one-loop check: computes F+, J-, and compares directly.
   This avoids stale exported result files. *)

projectDir = DirectoryName[DirectoryName[$InputFileName]];

Get[projectDir, "config", "LoadAMFlow.wl"];
Get[projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "OneLoopKernelUncutPlusFamily.wl"}]];
Get[FileNameJoin[{projectDir, "families", "OneLoopKernelUncutMinusFamily.wl"}]];

SetBasicAMFlowOptions[4];

(* Keep this in sync with coeffPowers and the GaugeLink output convention. *)
precisionGoal = 20;
epsOrder = 5;

(* Coefficients compared in the final Laurent expansion. *)
coeffPowers = {-1, 0, 1, 2, 3};

(* Must match the Numeric rules in both uncut family files. *)
xval = 1/3;
m2val = 1/10;
M2val = 1/5;
sval = 0;

Aval = (1 - xval) m2val + xval M2val - xval (1 - xval) sval;

analyticExpr =
  N[Normal[Series[Gamma[1 + eps]/eps * Aval^(-eps), {eps, 0, 3}]], 30];

Print["Full one-loop discontinuity check"];
Print["Precision goal: ", precisionGoal];
Print["Epsilon order passed to SolveIntegralsGaugeLink: ", epsOrder];
Print["Compared coefficient powers: ", coeffPowers];
Print["Kinematic point: x=", xval, ", m2=", m2val, ", M2=", M2val, ", s=p^2=", sval];
Print["Analytic A = (1-x)m2 + x M2 - x(1-x)s = ", Aval];
Print["Analytic expression:"];
Print[analyticExpr];

DefineOneLoopKernelUncutPlusFamily[];
plusTarget = j[oneloopkerneluncutplus, 1, 1, 1];

(* F+ contains 1/(L+i0). It is not physical by itself. *)
Print["Running +L GaugeLink integral: ", plusTarget];
plusResult = SolveIntegralsGaugeLink[{plusTarget}, precisionGoal, epsOrder];
Fplus = plusTarget /. plusResult;

Print["+L result:"];
Print[plusResult];

DefineOneLoopKernelUncutMinusFamily[];
minusTarget = j[oneloopkerneluncutminus, 1, 1, 1];

(* J- contains 1/(-L+i0). Together with F+ it gives the delta cut. *)
Print["Running -L GaugeLink integral: ", minusTarget];
minusResult = SolveIntegralsGaugeLink[{minusTarget}, precisionGoal, epsOrder];
Jminus = minusTarget /. minusResult;

Print["-L result:"];
Print[minusResult];

(* If Jminus is the integral with denominator (-L+i0), then
   1/(L-i0) = -1/(-L+i0). *)
discExpr = (-Jminus - Fplus)/(2 Pi I);

Print["Discontinuity expression (-Jminus - Fplus)/(2 Pi I):"];
Print[discExpr];

comparison =
  Table[
    With[
      {
        cAM = N[Coefficient[discExpr, eps, p], 30],
        cAN = N[Coefficient[analyticExpr, eps, p], 30]
      },
      <|
        "power" -> p,
        "amflow" -> cAM,
        "analytic" -> cAN,
        "difference" -> N[cAM - cAN, 30],
        "ratio" -> N[cAM/cAN, 30],
        ""
      |>
    ],
    {p, coeffPowers}
  ];

Print["Coefficient comparison:"];
Print[comparison];

(* Store the raw ingredients too, so the comparison can be audited later. *)
fullResult =
  <|
    "metadata" -> <|
      "description" -> "Fresh full one-loop discontinuity check; no pre-existing result files used.",
      "precisionGoal" -> precisionGoal,
      "epsOrder" -> epsOrder,
      "coeffPowers" -> coeffPowers,
      "x" -> xval,
      "m2" -> m2val,
      "M2" -> M2val,
      "s" -> sval,
      "A" -> Aval
    |>,
    "plusResult" -> plusResult,
    "minusResult" -> minusResult,
    "discontinuityExpression" -> discExpr,
    "analyticExpression" -> analyticExpr,
    "comparison" -> comparison
  |>;

Export[
  FileNameJoin[{$ResultsDirectory, "compare_oneloop_result.wl"}],
  fullResult
];

Export[
  FileNameJoin[{$ResultsDirectory, "compare_oneloop_result.txt"}],
  ToString[fullResult, InputForm],
  "Text"
];
