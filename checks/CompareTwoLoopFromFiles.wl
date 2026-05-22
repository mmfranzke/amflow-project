(* checks/CompareTwoLoopFromFiles.wl
   Post-processes previously exported PP, PM, MP, MM two-loop results.
   Use compare-twoloop for a fresh calculation and comparison. *)

ClearAll["Global`*"];

Print["Warning: this comparison reads previously exported two-loop result files."];
Print["Regenerate PP, PM, MP, MM after changing p2, masses, or epsOrder."];

projectDir =
  "/Users/FranzkeMM/Library/Mobile Documents/com~apple~CloudDocs/Documents/Dokumente privat/Studium/ETH Zürich/Master Thesis/amflow-project";

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

rawDiscExpr = Total[Values[pieceExprs]]/(2 Pi I)^2;

(* Same numerical point as TwoLoopKernelUncutFamilies.wl. *)
xval = 3/10;
yval = 1/4;
sval = 0;
ml2val = 49/100;
Ml2val = 4;
mk2val = 1/4;
Mk2val = 81/100;

Xval = 1 - xval;
lamval = yval/Xval;

delta0 =
  (1 - xval) ml2val + xval Ml2val - xval (1 - xval) sval;

delta1 =
  -((1 - lamval) mk2val + lamval Mk2val)
  + lamval (1 - lamval) Ml2val
  + lamval (1 - lamval xval) sval;

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

Print["Two-loop comparison from exported files:"];
Print[comparison];
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
      "doubleCutNormalization" -> doubleCutNormalization,
      "inferredLeadingNormalization" -> inferredLeadingNormalization
    |>,
    "pieceExpressions" -> pieceExprs,
    "rawDoubleDiscontinuityExpression" -> rawDiscExpr,
    "doubleDiscontinuityExpression" -> discExpr,
    "analyticExpression" -> analyticExpr,
    "comparison" -> comparison
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
      "comparison" -> comparison
    |>,
    InputForm
  ],
  "Text"
];
