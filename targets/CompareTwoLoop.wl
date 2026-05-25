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
Ml2val = 1;
mk2val = 1/4;
Mk2val = 81/100;

Xval = 1 - xval;
lamval = yval/Xval;

delta0 =
  (1 - xval) ml2val + xval Ml2val - xval (1 - xval) sval;

(* Current analytic comparison scale. *)
delta1Current =
  ((1 - lamval) mk2val + lamval Mk2val)
  - lamval (1 - lamval) Ml2val
  - lamval (1 - lamval xval) sval;

(* Effective scale from the AMFlow routing:
   D3 = k^2 - mk2, D4 = (l+k-p)^2 - Mk2, q = l-p,
   n.q = x-1, and lambda = y/(1-x). *)
delta1FromRouting =
  ((1 - lamval) mk2val + lamval Mk2val)
  - lamval (1 - lamval) Ml2val
  - lamval (1 - lamval xval) sval;

delta1 = delta1Current;

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
        Gamma[eps]^2/Xval
          (delta0 - I etaReg)^(-eps)
          (delta1 - I etaReg)^(-eps),
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
Print["lambda = y/(1-x): ", lamval];
Print["Delta0 analytic scale: ", delta0, " (", deltaSign["Delta0", delta0], ")"];
Print["Delta0 numeric: ", N[delta0, 30]];
Print["Delta1 currently used: ", delta1Current, " (", deltaSign["Delta1", delta1Current], ")"];
Print["Delta1 currently used numeric: ", N[delta1Current, 30]];
Print["Delta1 from routing formula: ", delta1FromRouting, " (", deltaSign["Delta1 routing", delta1FromRouting], ")"];
Print["Delta1 from routing numeric: ", N[delta1FromRouting, 30]];
Print["Delta1 current - routing: ", N[delta1Current - delta1FromRouting, 30]];
Print["Routing formula: ((1-lambda) mk2 + lambda Mk2) - lambda (1-lambda) Ml2 - lambda (1-lambda x) s"];
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

(* Product of the two one-dimensional discontinuity brackets:
   [-I_Mx - I_Px] [-I_My - I_Py] gives the sum of all four sign pieces. *)
rawDiscExpr =
  Total[Values[pieceExprs]]/(2 Pi I)^2;

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
      "Delta0" -> delta0,
      "Delta1" -> delta1,
      "Delta1Current" -> delta1Current,
      "Delta1FromRouting" -> delta1FromRouting,
      "Delta0Sign" -> deltaSign["Delta0", delta0],
      "Delta1Sign" -> deltaSign["Delta1", delta1],
      "doubleCutNormalization" -> doubleCutNormalization,
      "inferredLeadingNormalization" -> inferredLeadingNormalization
    |>,
    "pieceResults" -> pieceResults,
    "pieceExpressions" -> pieceExprs,
    "rawDoubleDiscontinuityExpression" -> rawDiscExpr,
    "doubleDiscontinuityExpression" -> discExpr,
    "analyticExpression" -> analyticExpr,
    "comparison" -> comparison
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
