(* checks/CompareTwoLoopFromFiles.wl
   Post-processes previously exported PP, PM, MP, MM two-loop results.
   Use compare-twoloop for a fresh calculation and comparison. *)

ClearAll["Global`*"];

Print["Warning: this comparison reads previously exported two-loop result files."];
Print["Regenerate PP, PM, MP, MM after changing p2, masses, or epsOrder."];

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

resultsDir = FileNameJoin[{projectDir, "results"}];

Get[FileNameJoin[{projectDir, "families", "TwoLoopKernelUncutFamilies.wl"}]];
Get[FileNameJoin[{projectDir, "targets", "TwoLoopKernelTargetTools.wl"}]];

pieces = TwoLoopKernelPieces[];

resultFile[piece_] :=
  FileNameJoin[{resultsDir, TwoLoopKernelResultBase[piece] <> ".wl"}];

Do[
  If[! FileExistsQ[resultFile[piece]],
    Print["Missing two-loop result file for piece ", piece, ": ", resultFile[piece]];
    Abort[];
  ],
  {piece, pieces}
];

pieceRules =
  Association[
    Table[
      piece -> Get[resultFile[piece]],
      {piece, pieces}
    ]
  ];

pieceExprs =
  Association[
    Table[
      piece -> (TwoLoopKernelTarget[piece] /. pieceRules[piece]),
      {piece, pieces}
    ]
  ];

Print["Two-loop piece expressions from files:"];
Do[
  Print[piece, " -> ", pieceExprs[piece]],
  {piece, pieces}
];

(* Stored M pieces use denominators -L+i0. Therefore
   integral[1/(L-i0)] = - stored M. *)
rawPieceCombination =
  pieceExprs["PP"] + pieceExprs["PM"] + pieceExprs["MP"] + pieceExprs["MM"];

(* Diagnostic only: this is the naive signed pattern for physical L-i0 pieces. *)
altSignedPieceCombination =
  pieceExprs["PP"] - pieceExprs["PM"] - pieceExprs["MP"] + pieceExprs["MM"];

discontinuityPrefactor = 1/(2 Pi I)^2;

Print["Raw stored-piece combination PP + PM + MP + MM:"];
Print[rawPieceCombination];
Print["Alternate diagnostic combination PP - PM - MP + MM:"];
Print[altSignedPieceCombination];
Print["Double-discontinuity prefactor:"];
Print[discontinuityPrefactor];

rawDiscExpr = discontinuityPrefactor rawPieceCombination;

(* Same numerical point as TwoLoopKernelUncutFamilies.wl. *)
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

delta1 =
  ((1 - lamval) mk2val + lamval Mk2val)
  - lamval (1 - lamval) Ml2val
  - lamval (1 - lamval xval) sval;

deltaSign[name_, value_] := Which[
  N[value] > 0, name <> " > 0",
  N[value] < 0, name <> " < 0",
  True, name <> " = 0"
];

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

(* With epsOrder = 4, compare poles and the finite coefficient. *)
coeffPowers = {-2, -1, 0};

(* Change to 1/2 if the leading-pole diagnostic remains raw/analytic = 1/2. *)
doubleCutNormalization = 1;
discExpr = rawDiscExpr/doubleCutNormalization;

leadingPower = First[coeffPowers];
inferredLeadingNormalization =
  N[
    Coefficient[rawDiscExpr, eps, leadingPower]/
      Coefficient[analyticExpr, eps, leadingPower],
    30
  ];

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

Print["Delta0: ", delta0, " (", deltaSign["Delta0", delta0], ")"];
Print["Delta1: ", delta1, " (", deltaSign["Delta1", delta1], ")"];
Print["Analytic expression:"];
Print[analyticExpr];
Print["Half analytic expression:"];
Print[analyticExpr/2];
Print["Two-loop comparison from exported files:"];
Print[comparison];
Print["Comparison against half of the analytic expression:"];
Print[halfAnalyticComparison];
Print["doubleCutNormalization used: ", doubleCutNormalization];
Print["Inferred leading-pole normalization raw/analytic: ", inferredLeadingNormalization];

Export[
  FileNameJoin[{resultsDir, "compare_twoloop_from_files_result.wl"}],
  <|
    "metadata" -> <|
      "description" -> "Post-processing comparison from exported two-loop sign pieces.",
      "x" -> xval,
      "y" -> yval,
      "s" -> sval,
      "ml2" -> ml2val,
      "Ml2" -> Ml2val,
      "mk2" -> mk2val,
      "Mk2" -> Mk2val,
      "Delta0" -> delta0,
      "Delta1" -> delta1,
      "Delta0Sign" -> deltaSign["Delta0", delta0],
      "Delta1Sign" -> deltaSign["Delta1", delta1],
      "doubleCutNormalization" -> doubleCutNormalization,
      "inferredLeadingNormalization" -> inferredLeadingNormalization
    |>,
    "pieceExpressions" -> pieceExprs,
    "rawStoredPieceCombination" -> rawPieceCombination,
    "alternateSignedPieceCombination" -> altSignedPieceCombination,
    "doubleDiscontinuityPrefactor" -> discontinuityPrefactor,
    "rawDoubleDiscontinuityExpression" -> rawDiscExpr,
    "doubleDiscontinuityExpression" -> discExpr,
    "analyticExpression" -> analyticExpr,
    "comparison" -> comparison,
    "halfAnalyticComparison" -> halfAnalyticComparison
  |>
];

Export[
  FileNameJoin[{resultsDir, "compare_twoloop_from_files_result.txt"}],
  ToString[
    <|
      "pieceExpressions" -> pieceExprs,
      "rawDoubleDiscontinuityExpression" -> rawDiscExpr,
      "doubleDiscontinuityExpression" -> discExpr,
      "analyticExpression" -> analyticExpr,
      "comparison" -> comparison,
      "halfAnalyticComparison" -> halfAnalyticComparison
    |>,
    InputForm
  ],
  "Text"
];
