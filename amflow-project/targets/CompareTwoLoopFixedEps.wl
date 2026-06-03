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
residueScaleEnv = Environment["AMFLOW_RESIDUE_SCALE"];
residueScaleComparison =
  If[StringQ[residueScaleEnv] && StringLength[residueScaleEnv] > 0,
    residueScaleEnv,
    "not specified"
  ];
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
pointMetadata = KeyTake[
  point,
  {"name", "pPlus", "pMinus", "p2", "pPerp2", "x", "y", "ml2", "Ml2", "mk2", "Mk2"}
];
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
delta1OldCompare = (
  -muKlambdaVal
  + lambdaVal (1 - lambdaVal) Ml2val
  + lambdaVal (1 - xval) sval
);
delta1 = (
  -muKlambdaVal
  + lambdaVal (1 - lambdaVal) Ml2val
  + lambdaVal (1 - lambdaVal*xval) sval
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
  cachedPointMetadata =
    If[KeyExistsQ[cached, "pointMetadata"],
      cached["pointMetadata"],
      If[KeyExistsQ[cached, "pointName"],
        <|"name" -> cached["pointName"]|>,
        <||>
      ]
    ];
  If[KeyExistsQ[cachedPointMetadata, "name"] && cachedPointMetadata["name"] =!= point["name"],
    Print[
      "AMFLOW_CACHE_POINT_MISMATCH cached=",
      cachedPointMetadata["name"],
      " requested=",
      point["name"]
    ];
    Abort[];
  ];
  If[KeyExistsQ[cached, "pointMetadata"] && cached["pointMetadata"] =!= pointMetadata,
    Print["AMFLOW_CACHE_METADATA_MISMATCH"];
    Print["cached point metadata: ", InputForm[cached["pointMetadata"]]];
    Print["requested point metadata: ", InputForm[pointMetadata]];
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
  cachedCoeffPowers =
    If[KeyExistsQ[cached, "coefficientPowers"],
      cached["coefficientPowers"],
      If[KeyExistsQ[cached, "doubleDiscontinuityCoefficients"],
        Keys[cached["doubleDiscontinuityCoefficients"]],
        Range[-2, 0]
      ]
    ];
  missingCoeffPowers = Complement[coeffPowers, cachedCoeffPowers];
  If[Length[missingCoeffPowers] > 0,
    Print[
      "AMFLOW_CACHE_COEFFICIENTS_UNAVAILABLE cached powers=",
      InputForm[cachedCoeffPowers],
      " requested powers=",
      InputForm[coeffPowers]
    ];
    Print[
      "Cached AMFlow result has coefficients ",
      InputForm[cachedCoeffPowers],
      ". Rerun with --amflow fresh --amflow-eps-order ",
      epsOrder,
      " or larger."
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
method12CoeffAssoc =
  Association[
    Table[
      power -> N[
        SeriesCoefficient[
          Gamma[eps]^2/Xval *
            (delta0 - I etaReg)^(-eps) *
            (delta1 - I etaReg)^(-eps),
          {eps, 0, power}
        ],
        40
      ],
      {power, coeffPowers}
    ]
  ];
comboExprs =
  If[AssociationQ[pieceExprs] && And @@ (KeyExistsQ[pieceExprs, #] & /@ {"PP", "PM", "MP", "MM"}),
    <|
      "combo_all_plus" -> pieceExprs["PP"] + pieceExprs["PM"] + pieceExprs["MP"] + pieceExprs["MM"],
      "combo_delta_standard" -> pieceExprs["PP"] - pieceExprs["PM"] - pieceExprs["MP"] + pieceExprs["MM"],
      "combo_delta_alt1" -> pieceExprs["PP"] - pieceExprs["PM"] + pieceExprs["MP"] - pieceExprs["MM"],
      "combo_delta_alt2" -> pieceExprs["PP"] + pieceExprs["PM"] - pieceExprs["MP"] - pieceExprs["MM"],
      "combo_single_x_disc" -> pieceExprs["PP"] - pieceExprs["MP"] + pieceExprs["PM"] - pieceExprs["MM"],
      "combo_single_y_disc" -> pieceExprs["PP"] - pieceExprs["PM"] + pieceExprs["MP"] - pieceExprs["MM"]
    |>,
    <||>
  ];
comboNormalizedCoeffAssoc =
  Association[
    KeyValueMap[
      Function[
        {name, expr},
        name -> Association[
          Table[
            power -> N[SeriesCoefficient[expr/(2 Pi I)^2, {eps, 0, power}], 40],
            {power, coeffPowers}
          ]
        ]
      ],
      comboExprs
    ]
  ];
comboScores =
  Association[
    KeyValueMap[
      Function[
        {name, coeffs},
        name -> N[
          Sqrt[
            Total[
              Table[
                Abs[coeffs[power] - method12CoeffAssoc[power]]^2,
                {power, coeffPowers}
              ]
            ]
          ],
          40
        ]
      ],
      comboNormalizedCoeffAssoc
    ]
  ];

Print["Fixed-eps two-loop comparison"];
TwoLoopPrintPoint[point];
Print["eps = ", epsValue];
Print["Delta0 = ", delta0];
Print["Delta1 current PDF = ", delta1];
Print["Delta1 old compare diagnostic = ", delta1OldCompare];
If[sval =!= 0,
  Print["WARNING: off-shell point is sensitive to the Delta1 p2-term convention."];
  Print["Current-PDF default uses lambda*(1-lambda*x)*p2."];
  Print["Old compare convention used lambda*(1-x)*p2."];
];
Print["METHOD12_FIXED_EPS_VALUE=", InputForm[method12Value]];
Print["AMFLOW_FIXED_EPS_VALUE=", InputForm[amflowValue]];
Print["BEGIN_AMFLOW_METADATA"];
Print["point = ", point["name"]];
Print["pPlus = ", InputForm[point["pPlus"]]];
Print["pMinus = ", InputForm[point["pMinus"]]];
Print["pPerp2 = ", InputForm[point["pPerp2"]]];
Print["p2 = ", InputForm[point["p2"]]];
Print["x = ", InputForm[point["x"]]];
Print["y = ", InputForm[point["y"]]];
Print["ml2 = ", InputForm[point["ml2"]]];
Print["Ml2 = ", InputForm[point["Ml2"]]];
Print["mk2 = ", InputForm[point["mk2"]]];
Print["Mk2 = ", InputForm[point["Mk2"]]];
Print["AMFLOW_EPS_ORDER = ", InputForm[epsOrder]];
Print["AMFLOW_PRECISION_GOAL = ", InputForm[precisionGoal]];
Print["AMFLOW_TARGET_INTEGRALS = ", InputForm[TwoLoopKernelTarget /@ TwoLoopKernelPieces[]]];
Print["AMFLOW_NORMALIZATION = 1/(2 Pi I)^2"];
Print["residue_scale_comparison = ", residueScaleComparison];
Print["END_AMFLOW_METADATA"];
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
Print["AMFLOW_SIGN_PATTERN_DIAGNOSTIC=default object remains combo_all_plus; signed delta patterns are diagnostic only"];
Print["AMFLOW_SIGN_PATTERN_combo_all_plus=PP+PM+MP+MM"];
Print["AMFLOW_SIGN_PATTERN_combo_delta_standard=PP-PM-MP+MM, corresponding to [P-M]_x [P-M]_y if P/M are opposite i0 sides"];
Print["AMFLOW_SIGN_PATTERN_combo_delta_alt1=PP-PM+MP-MM, diagnostic alternate single-axis convention"];
Print["AMFLOW_SIGN_PATTERN_combo_delta_alt2=PP+PM-MP-MM, diagnostic alternate single-axis convention"];
Print["AMFLOW_SIGN_PATTERN_combo_single_x_disc=PP-MP+PM-MM, x-linear denominator discontinuity diagnostic"];
Print["AMFLOW_SIGN_PATTERN_combo_single_y_disc=PP-PM+MP-MM, y-linear denominator discontinuity diagnostic"];
Print["AMFLOW_COEFF_POWERS=", InputForm[coeffPowers]];
Do[
  Print[
    "AMFLOW_COEFF_POWER=", power,
    " RAW=", InputForm[rawCoeffAssoc[power]],
    " NORMALIZED=", InputForm[normalizedCoeffAssoc[power]]
  ],
  {power, coeffPowers}
];
KeyValueMap[
  Function[
    {comboName, coeffs},
    Do[
      Print[
        "AMFLOW_COMBO_COEFF NAME=", comboName,
        " POWER=", power,
        " VALUE=", InputForm[coeffs[power]],
        " METHOD12=", InputForm[method12CoeffAssoc[power]],
        " DIFF=", InputForm[N[coeffs[power] - method12CoeffAssoc[power], 40]],
        " RATIO=", InputForm[N[coeffs[power]/method12CoeffAssoc[power], 40]],
        " SCORE=", InputForm[comboScores[comboName]]
      ],
      {power, coeffPowers}
    ]
  ],
  comboNormalizedCoeffAssoc
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
    "pointMetadata" -> pointMetadata,
    "eps" -> epsValue,
    "epsList" -> epsList,
    "delta0" -> delta0,
    "delta1" -> delta1,
    "delta1OldCompare" -> delta1OldCompare,
    "method12Closed" -> method12Value,
    "amflowOriginal" -> amflowValue,
    "amflowEpsListValues" -> amflowEpsListValues,
    "AMFLOW_EPS_ORDER" -> epsOrder,
    "coefficientPowers" -> coeffPowers,
    "selectedFamilies" -> (TwoLoopKernelFamilySymbol /@ TwoLoopKernelPieces[]),
    "selectedIntegrals" -> (TwoLoopKernelTarget /@ TwoLoopKernelPieces[]),
    "normalizationFactor" -> HoldForm[1/(2 Pi I)^2],
    "rawPieceCombinationCoefficients" -> rawCoeffAssoc,
    "doubleDiscontinuityCoefficients" -> normalizedCoeffAssoc,
    "method12Coefficients" -> method12CoeffAssoc,
    "doubleDiscontinuitySignPatternDiagnostics" -> comboNormalizedCoeffAssoc,
    "doubleDiscontinuitySignPatternScores" -> comboScores,
    "pieceExpressions" -> pieceExprs,
    "doubleDiscontinuityExpression" -> discExpr
  |>
];

Export[
  legacyCacheFile,
  <|
    "pointName" -> point["name"],
    "pointMetadata" -> pointMetadata,
    "eps" -> epsValue,
    "epsList" -> epsList,
    "delta0" -> delta0,
    "delta1" -> delta1,
    "delta1OldCompare" -> delta1OldCompare,
    "method12Closed" -> method12Value,
    "amflowOriginal" -> amflowValue,
    "amflowEpsListValues" -> amflowEpsListValues,
    "AMFLOW_EPS_ORDER" -> epsOrder,
    "coefficientPowers" -> coeffPowers,
    "selectedFamilies" -> (TwoLoopKernelFamilySymbol /@ TwoLoopKernelPieces[]),
    "selectedIntegrals" -> (TwoLoopKernelTarget /@ TwoLoopKernelPieces[]),
    "normalizationFactor" -> HoldForm[1/(2 Pi I)^2],
    "rawPieceCombinationCoefficients" -> rawCoeffAssoc,
    "doubleDiscontinuityCoefficients" -> normalizedCoeffAssoc,
    "method12Coefficients" -> method12CoeffAssoc,
    "doubleDiscontinuitySignPatternDiagnostics" -> comboNormalizedCoeffAssoc,
    "doubleDiscontinuitySignPatternScores" -> comboScores,
    "pieceExpressions" -> pieceExprs,
    "doubleDiscontinuityExpression" -> discExpr
  |>
];
