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
coeffPowers = TwoLoopCoeffPowers[epsOrder];
epsEnv = Environment["TWO_LOOP_EPS_VALUE"];
epsListEnv = Environment["TWO_LOOP_EPS_LIST"];
reuseEnv = Environment["AMFLOW_REUSE_RESULT"];
reuseAMFlowResult = StringQ[reuseEnv] && MemberQ[{"1", "true", "True", "yes", "YES"}, reuseEnv];
legacyCacheFile = FileNameJoin[{$ResultsDirectory, "compare_twoloop_fixed_eps_result.wl"}];
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
cacheFile = FileNameJoin[
  {$ResultsDirectory, "compare_twoloop_fixed_eps_result_" <> point["name"] <> ".wl"}
];
If[reuseAMFlowResult && ! FileExistsQ[cacheFile] && FileExistsQ[legacyCacheFile],
  cacheFile = legacyCacheFile;
];
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

If[reuseAMFlowResult,
  If[! FileExistsQ[cacheFile],
    Print["AMFLOW_CACHE_MISS=", cacheFile];
    Abort[];
  ];
  cached = Import[cacheFile];
  If[! AssociationQ[cached] || ! KeyExistsQ[cached, "doubleDiscontinuityExpression"],
    Print["AMFLOW_CACHE_INVALID=", cacheFile];
    Abort[];
  ];
  If[KeyExistsQ[cached, "pointName"] && cached["pointName"] =!= point["name"],
    Print[
      "AMFLOW_CACHE_POINT_MISMATCH cached=",
      cached["pointName"],
      " requested=",
      point["name"]
    ];
    Abort[];
  ];
  pieceExprs = If[KeyExistsQ[cached, "pieceExpressions"], cached["pieceExpressions"], <||>];
  discExpr = cached["doubleDiscontinuityExpression"];
  rawPieceCombination = (2 Pi I)^2 discExpr;
  Print["AMFLOW_REUSED_CACHE=", cacheFile],

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
  Print["AMFLOW_FRESH_SOLVE_DONE=1"];
];
amflowValue = N[discExpr /. eps -> epsValue, 40];
amflowEpsListValues = N[discExpr /. eps -> #, 40] & /@ epsList;
rawCoeffAssoc =
  Association[
    Table[
      power -> N[SeriesCoefficient[rawPieceCombination, {eps, 0, power}], 40],
      {power, coeffPowers}
    ]
  ];
normalizedCoeffAssoc =
  Association[
    Table[
      power -> N[SeriesCoefficient[discExpr, {eps, 0, power}], 40],
      {power, coeffPowers}
    ]
  ];

Print["Fixed-eps two-loop comparison"];
TwoLoopPrintPoint[point];
Print["eps = ", epsValue];
Print["Delta0 = ", delta0];
Print["Delta1 = ", delta1];
Print["METHOD12_FIXED_EPS_VALUE=", InputForm[method12Value]];
Print["AMFLOW_FIXED_EPS_VALUE=", InputForm[amflowValue]];
Print["AMFLOW_OBJECT_POINT=", point["name"]];
Print["AMFLOW_OBJECT_TARGET=compare-twoloop-fixed-eps"];
Print["AMFLOW_OBJECT_FAMILIES=", InputForm[TwoLoopKernelFamilySymbol /@ TwoLoopKernelPieces[]]];
Print["AMFLOW_OBJECT_PIECES=", InputForm[TwoLoopKernelPieces[]]];
Print["AMFLOW_OBJECT_TARGET_INTEGRALS=", InputForm[TwoLoopKernelTarget /@ TwoLoopKernelPieces[]]];
Print["AMFLOW_OBJECT_PROPAGATOR_EXPONENTS=", InputForm[{1, 1, 1, 1, 1, 1, 0}]];
Print["AMFLOW_OBJECT_LINEAR_DENOMINATORS={sx (x - n.l), sy (y - n.k)}"];
Print["AMFLOW_OBJECT_PRESCRIPTION=PP+PM+MP+MM uncut GaugeLink double discontinuity"];
Print["AMFLOW_OBJECT_RAW_COMBINATION=pieceExprs[PP]+pieceExprs[PM]+pieceExprs[MP]+pieceExprs[MM]"];
Print["AMFLOW_OBJECT_NORMALIZATION_FACTOR=1/(2 Pi I)^2"];
Print["AMFLOW_OBJECT_REAL_PART_TAKEN=False"];
Print["AMFLOW_OBJECT_SUMMED_PIECES=True"];
Print["AMFLOW_OBJECT_DOUBLE_DISCONTINUITY=raw combination/(2 Pi I)^2"];
Print["AMFLOW_COEFF_POWERS=", InputForm[coeffPowers]];
Do[
  Print[
    "AMFLOW_COEFF_POWER=", power,
    " RAW=", InputForm[rawCoeffAssoc[power]],
    " NORMALIZED=", InputForm[normalizedCoeffAssoc[power]]
  ],
  {power, coeffPowers}
];
Do[
  Print[
    "AMFLOW_EPS_LIST_VALUE_INDEX=", i,
    " VALUE=", InputForm[amflowEpsListValues[[i]]]
  ],
  {i, Length[amflowEpsListValues]}
];

Export[
  cacheFile,
  <|
    "pointName" -> point["name"],
    "eps" -> epsValue,
    "epsList" -> epsList,
    "delta0" -> delta0,
    "delta1" -> delta1,
    "method12Closed" -> method12Value,
    "amflowOriginal" -> amflowValue,
    "amflowEpsListValues" -> amflowEpsListValues,
    "rawPieceCombinationCoefficients" -> rawCoeffAssoc,
    "doubleDiscontinuityCoefficients" -> normalizedCoeffAssoc,
    "pieceExpressions" -> pieceExprs,
    "doubleDiscontinuityExpression" -> discExpr
  |>
];

Export[
  legacyCacheFile,
  <|
    "pointName" -> point["name"],
    "eps" -> epsValue,
    "epsList" -> epsList,
    "delta0" -> delta0,
    "delta1" -> delta1,
    "method12Closed" -> method12Value,
    "amflowOriginal" -> amflowValue,
    "amflowEpsListValues" -> amflowEpsListValues,
    "rawPieceCombinationCoefficients" -> rawCoeffAssoc,
    "doubleDiscontinuityCoefficients" -> normalizedCoeffAssoc,
    "pieceExpressions" -> pieceExprs,
    "doubleDiscontinuityExpression" -> discExpr
  |>
];
