(* targets/CompareTwoLoop.wl
   Fresh two-loop check: computes PP, PM, MP, MM and compares the double
   discontinuity to the analytic closed form at one light-like point. *)

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "TwoLoopKernelUncutFamilies.wl"}]];
Get[FileNameJoin[{projectDir, "targets", "TwoLoopKernelTargetTools.wl"}]];

SetBasicAMFlowOptions[4];

(* Override with AMFLOW_PRECISION_GOAL and AMFLOW_EPS_ORDER for quick tests. *)
precisionGoal = TwoLoopPrecisionGoal[10];
epsOrder = TwoLoopEpsOrder[4];
coeffPowers = TwoLoopCoeffPowers[epsOrder];

(* Analytic comparison point. It must match TwoLoopKernelUncutFamilies.wl. *)
xval = 3/10;
yval = 1/4;
sval = 0;
ml2val = 49/100;
Ml2val = 4;
mk2val = 1/4;
Mk2val = 81/100;

Xval = 1 - xval;
lambdaVal = yval/Xval;
aVal = lambdaVal (1 - lambdaVal);
bVal = lambdaVal (1 - lambdaVal xval);

muLxVal = (1 - xval) ml2val + xval Ml2val;
muKlambdaVal = (1 - lambdaVal) mk2val + lambdaVal Mk2val;

delta0 = muLxVal - xval (1 - xval) sval;

deltaW = (
  ((1 + aVal) xval (1 - xval) + bVal) sval
  - muKlambdaVal
  + aVal Ml2val
  - (1 + aVal) muLxVal
);

delta1FromDeltaW = (1 + aVal) delta0 + deltaW;

delta1Final = (
  -muKlambdaVal
  + lambdaVal (1 - lambdaVal) Ml2val
  + lambdaVal (1 - lambdaVal xval) sval
);

delta1ConsistencyDifference = FullSimplify[delta1FromDeltaW - delta1Final];

deltaSign[name_, value_] := Which[
  N[value] > 0, name <> " > 0",
  N[value] < 0, name <> " < 0",
  True, name <> " = 0"
];

(* Tiny -i0 regulator fixes the branch for negative scales. *)
etaReg = 10^-30;

analyticExpr =
  N[
    Normal[
      Series[
        Gamma[eps]^2/Xval *
          (delta0 - I etaReg)^(-eps) *
          (delta1Final - I etaReg)^(-eps),
        {eps, 0, 0}
      ]
    ],
    30
  ];

Print["Fresh two-loop double-discontinuity check"];
Print["Precision goal: ", precisionGoal];
Print["Epsilon order: ", epsOrder];
Print["Compared coefficient powers: ", coeffPowers];
Print["Kinematic point: x=", xval, ", y=", yval, ", s=p^2=", sval];
Print["Mass squares: ml2=", ml2val, ", Ml2=", Ml2val, ", mk2=", mk2val, ", Mk2=", Mk2val];
Print["lambda = y/(1-x): ", lambdaVal];
Print["a = lambda (1-lambda): ", aVal];
Print["b = lambda (1-lambda x): ", bVal];
Print["muLx = (1-x) ml2 + x Ml2: ", muLxVal];
Print["muKlambda = (1-lambda) mk2 + lambda Mk2: ", muKlambdaVal];
Print["delta0: ", delta0, " (", deltaSign["delta0", delta0], ")"];
Print["delta0 numeric: ", N[delta0, 30]];
Print["deltaW: ", deltaW];
Print["delta1FromDeltaW = (1+a) delta0 + deltaW: ", delta1FromDeltaW];
Print["delta1Final = -muKlambda + lambda (1-lambda) Ml2 + lambda (1-lambda x) s: ", delta1Final, " (", deltaSign["delta1Final", delta1Final], ")"];
Print["delta1FromDeltaW - delta1Final: ", delta1ConsistencyDifference];
If[delta1ConsistencyDifference =!= 0,
  Print["WARNING: delta1FromDeltaW != delta1Final after FullSimplify"];
];
Print["Check that the family Numeric point printed below matches this analytic point."];
Print["Analytic expression:"];
Print[analyticExpr];

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

Print["Two-loop piece expressions:"];
Do[
  Print[piece, " -> ", pieceExprs[piece]],
  {piece, pieces}
];

(* Stored M pieces use denominators -L+i0. Therefore
   integral[1/(L-i0)] = - stored M. The physical bracket
   (I_x^- - I_x^+) (I_y^- - I_y^+) is represented by the stored sum. *)
rawPieceCombination =
  pieceExprs["PP"] + pieceExprs["PM"] + pieceExprs["MP"] + pieceExprs["MM"];

(* Diagnostic only: this would be the sign pattern if PM and MP were already
   physical L-i0 pieces rather than stored -L+i0 pieces. *)
altSignedPieceCombination =
  pieceExprs["PP"] - pieceExprs["PM"] - pieceExprs["MP"] + pieceExprs["MM"];

discontinuityPrefactor = 1/(2 Pi I)^2;

Print["Raw stored-piece combination PP + PM + MP + MM:"];
Print[rawPieceCombination];
Print["Alternate diagnostic combination PP - PM - MP + MM:"];
Print[altSignedPieceCombination];
Print["Double-discontinuity prefactor:"];
Print[discontinuityPrefactor];

rawDiscExpr =
  discontinuityPrefactor rawPieceCombination;

(* Keep this explicit. The old s=-3 run showed a leading-pole ratio 1/2,
   so rerun compare-twoloop-from-files with doubleCutNormalization = 1/2
   if the new light-like point shows the same constant factor. *)
doubleCutNormalization = 1;
discExpr = rawDiscExpr/doubleCutNormalization;

leadingPower = First[coeffPowers];
inferredLeadingNormalization =
  N[
    Coefficient[rawDiscExpr, eps, leadingPower]/
      Coefficient[analyticExpr, eps, leadingPower],
    30
  ];

Print["Double-discontinuity expression:"];
Print[discExpr];
Print["doubleCutNormalization used: ", doubleCutNormalization];
Print["Inferred leading-pole normalization raw/analytic: ", inferredLeadingNormalization];
Print["Analytic expression:"];
Print[analyticExpr];
Print["Half analytic expression (diagnostic only, not physics):"];
Print[analyticExpr/2];

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

Print["Coefficient comparison:"];
Print[comparison];

halfAnalyticComparison =
  Table[
    With[
      {
        cAM = N[Coefficient[discExpr, eps, p], 30],
        cHalf = N[Coefficient[analyticExpr/2, eps, p], 30]
      },
      <|
        "power" -> p,
        "amflow" -> cAM,
        "halfAnalytic" -> cHalf,
        "differenceFromHalfAnalytic" -> N[cAM - cHalf, 30],
        "ratioToHalfAnalytic" -> N[cAM/cHalf, 30]
      |>
    ],
    {p, coeffPowers}
  ];

Print["Comparison against half of the analytic expression (diagnostic only, not physics):"];
Print[halfAnalyticComparison];

fullResult =
  <|
    "metadata" -> <|
      "description" -> "Fresh two-loop double-discontinuity check; no pre-existing result files used.",
      "precisionGoal" -> precisionGoal,
      "epsOrder" -> epsOrder,
      "coeffPowers" -> coeffPowers,
      "x" -> xval,
      "y" -> yval,
      "s" -> sval,
      "ml2" -> ml2val,
      "Ml2" -> Ml2val,
      "mk2" -> mk2val,
      "Mk2" -> Mk2val,
      "lambda" -> lambdaVal,
      "a" -> aVal,
      "b" -> bVal,
      "muLx" -> muLxVal,
      "muKlambda" -> muKlambdaVal,
      "delta0" -> delta0,
      "deltaW" -> deltaW,
      "delta1FromDeltaW" -> delta1FromDeltaW,
      "delta1Final" -> delta1Final,
      "delta1ConsistencyDifference" -> delta1ConsistencyDifference,
      "delta0Sign" -> deltaSign["delta0", delta0],
      "delta1FinalSign" -> deltaSign["delta1Final", delta1Final],
      "doubleCutNormalization" -> doubleCutNormalization,
      "inferredLeadingNormalization" -> inferredLeadingNormalization
    |>,
    "pieceResults" -> pieceResults,
    "pieceExpressions" -> pieceExprs,
    "rawStoredPieceCombination" -> rawPieceCombination,
    "alternateSignedPieceCombination" -> altSignedPieceCombination,
    "doubleDiscontinuityPrefactor" -> discontinuityPrefactor,
    "rawDoubleDiscontinuityExpression" -> rawDiscExpr,
    "doubleDiscontinuityExpression" -> discExpr,
    "analyticExpression" -> analyticExpr,
    "comparison" -> comparison,
    "halfAnalyticComparison" -> halfAnalyticComparison
  |>;

Export[
  FileNameJoin[{$ResultsDirectory, "compare_twoloop_result.wl"}],
  fullResult
];

Export[
  FileNameJoin[{$ResultsDirectory, "compare_twoloop_result.txt"}],
  ToString[fullResult, InputForm],
  "Text"
];
