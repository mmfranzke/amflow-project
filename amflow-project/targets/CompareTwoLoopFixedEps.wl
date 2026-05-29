(* targets/CompareTwoLoopFixedEps.wl
   Fresh two-loop check evaluated at one fixed positive epsilon. *)

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "config", "TwoLoopKinematicPoints.wl"}]];
$TwoLoopDefaultPointName = "equal_mass_onshell_branch";
Get[FileNameJoin[{projectDir, "families", "TwoLoopKernelUncutFamilies.wl"}]];
Get[FileNameJoin[{projectDir, "targets", "TwoLoopKernelTargetTools.wl"}]];

SetBasicAMFlowOptions[4];

precisionGoal = TwoLoopPrecisionGoal[10];
epsOrder = TwoLoopEpsOrder[4];
epsEnv = Environment["TWO_LOOP_EPS_VALUE"];
epsListEnv = Environment["TWO_LOOP_EPS_LIST"];
epsValue =
  If[StringQ[epsEnv] && StringLength[epsEnv] > 0,
    ToExpression[epsEnv],
    1/10
  ];
If[! NumericQ[N[epsValue]],
  epsValue = 1/10;
];
epsList =
  If[StringQ[epsListEnv] && StringLength[epsListEnv] > 0,
    ToExpression /@ StringSplit[epsListEnv, ","],
    {epsValue}
  ];

point = TwoLoopSelectedPoint[];
xval = point["x"];
yval = point["y"];
sval = point["p2"];
ml2val = point["ml2"];
Ml2val = point["Ml2"];
mk2val = point["mk2"];
Mk2val = point["Mk2"];

Xval = 1 - xval;
lambdaVal = yval/Xval;

muLxVal = (1 - xval) ml2val + xval Ml2val;
muKlambdaVal = (1 - lambdaVal) mk2val + lambdaVal Mk2val;

delta0 = muLxVal - xval (1 - xval) sval;
delta1 = (
  -muKlambdaVal
  + lambdaVal (1 - lambdaVal) Ml2val
  + lambdaVal (1 - xval) sval
);

etaReg = 10^-30;
method12Value =
  N[
    Gamma[epsValue]^2/Xval *
      (delta0 - I etaReg)^(-epsValue) *
      (delta1 - I etaReg)^(-epsValue),
    40
  ];

pieces = TwoLoopKernelPieces[];
pieceResults =
  Association[
    Table[
      piece -> SolveTwoLoopKernelPiece[piece, precisionGoal, epsOrder],
      {piece, pieces}
    ]
  ];

pieceExprs =
  Association[
    Table[
      piece -> (TwoLoopKernelTarget[piece] /. pieceResults[piece]),
      {piece, pieces}
    ]
  ];

rawPieceCombination =
  pieceExprs["PP"] + pieceExprs["PM"] + pieceExprs["MP"] + pieceExprs["MM"];

discExpr = rawPieceCombination/(2 Pi I)^2;
amflowValue = N[discExpr /. eps -> epsValue, 40];
amflowEpsListValues = N[discExpr /. eps -> #, 40] & /@ epsList;

Print["Fixed-eps two-loop comparison"];
TwoLoopPrintPoint[point];
Print["eps = ", epsValue];
Print["Delta0 = ", delta0];
Print["Delta1 = ", delta1];
Print["METHOD12_FIXED_EPS_VALUE=", InputForm[method12Value]];
Print["AMFLOW_FIXED_EPS_VALUE=", InputForm[amflowValue]];
Do[
  Print[
    "AMFLOW_EPS_LIST_VALUE_INDEX=", i,
    " VALUE=", InputForm[amflowEpsListValues[[i]]]
  ],
  {i, Length[amflowEpsListValues]}
];

Export[
  FileNameJoin[{$ResultsDirectory, "compare_twoloop_fixed_eps_result.wl"}],
  <|
    "pointName" -> point["name"],
    "eps" -> epsValue,
    "epsList" -> epsList,
    "delta0" -> delta0,
    "delta1" -> delta1,
    "method12Closed" -> method12Value,
    "amflowOriginal" -> amflowValue,
    "amflowEpsListValues" -> amflowEpsListValues,
    "pieceExpressions" -> pieceExprs,
    "doubleDiscontinuityExpression" -> discExpr
  |>
];
