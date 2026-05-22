(* targets/CompareOneLoopAuxDen.wl
   Fresh one-loop check with an auxiliary denominator at exponent 0.
   Computes F+, J-, forms the discontinuity, and compares to the analytic result. *)

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "OneLoopKernelAuxDenFamilies.wl"}]];

SetBasicAMFlowOptions[4];

precisionGoal = 20;
epsOrder = 5;
coeffPowers = {-1, 0, 1, 2, 3};

xval = 1/3;
m2val = 1/10;
M2val = 1/5;
sval = 0;

Aval = (1 - xval) m2val + xval M2val - xval (1 - xval) sval;

analyticExpr =
  N[Normal[Series[Gamma[1 + eps]/eps * Aval^(-eps), {eps, 0, 3}]], 30];

DefineOneLoopKernelAuxDenPlusFamily[];
plusTarget = j[oneloopkernelauxplus, 1, 1, 1, 0];

Print["Running auxiliary-denominator +L integral: ", plusTarget];
plusResult = SolveIntegralsGaugeLink[{plusTarget}, precisionGoal, epsOrder];
Fplus = plusTarget /. plusResult;

DefineOneLoopKernelAuxDenMinusFamily[];
minusTarget = j[oneloopkernelauxminus, 1, 1, 1, 0];

Print["Running auxiliary-denominator -L integral: ", minusTarget];
minusResult = SolveIntegralsGaugeLink[{minusTarget}, precisionGoal, epsOrder];
Jminus = minusTarget /. minusResult;

discExpr = (-Jminus - Fplus)/(2 Pi I);

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
        "ratio" -> N[cAM/cAN, 30]
      |>
    ],
    {p, coeffPowers}
  ];

Print["Auxiliary-denominator one-loop comparison:"];
Print[comparison];

fullResult =
  <|
    "metadata" -> <|
      "description" -> "Fresh one-loop auxiliary-denominator check.",
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
  FileNameJoin[{$ResultsDirectory, "compare_oneloop_aux_result.wl"}],
  fullResult
];

Export[
  FileNameJoin[{$ResultsDirectory, "compare_oneloop_aux_result.txt"}],
  ToString[fullResult, InputForm],
  "Text"
];
