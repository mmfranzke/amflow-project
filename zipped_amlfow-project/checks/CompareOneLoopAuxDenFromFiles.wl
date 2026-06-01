(* checks/CompareOneLoopAuxDenFromFiles.wl
   Post-processes exported one-loop auxiliary-denominator F+ and J- results.
   Use compare-oneloop-aux for a fresh calculation and comparison. *)

ClearAll["Global`*"];

Print["Warning: this comparison reads previously exported auxiliary +L and -L result files."];
Print["For a fresh one-command calculation and comparison, use ./run.sh compare-oneloop-aux."];

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

resultsDir = FileNameJoin[{projectDir, "results"}];

(* These must be regenerated after changing kinematics or epsOrder. *)
plusFile = FileNameJoin[{resultsDir, "oneloop_kernel_aux_plus_result.wl"}];
minusFile = FileNameJoin[{resultsDir, "oneloop_kernel_aux_minus_result.wl"}];

If[! FileExistsQ[plusFile],
  Print["Missing auxiliary +L result file: ", plusFile];
  Abort[];
];

If[! FileExistsQ[minusFile],
  Print["Missing auxiliary -L result file: ", minusFile];
  Abort[];
];

plusRules = Get[plusFile];
minusRules = Get[minusFile];

(* D4 is the auxiliary denominator and was requested with exponent 0. *)
Fplus = j[oneloopkernelauxplus, 1, 1, 1, 0] /. plusRules;
Jminus = j[oneloopkernelauxminus, 1, 1, 1, 0] /. minusRules;

(* If Jminus = integral 1/(-L+i0), then integral 1/(L-i0) = -Jminus. *)
discExpr = (-Jminus - Fplus)/(2 Pi I);

xval = 1/3;
m2val = 1/10;
M2val = 1/5;
sval = 0;

Aval = (1 - xval) m2val + xval M2val - xval (1 - xval) sval;

analyticExpr =
  N[Normal[Series[Gamma[1 + eps]/eps * Aval^(-eps), {eps, 0, 3}]], 30];

coeffPowers = {-1, 0, 1, 2, 3};

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

Print["Auxiliary-denominator comparison from files:"];
Print[comparison];

Export[
  FileNameJoin[{resultsDir, "compare_oneloop_aux_from_files_result.wl"}],
  <|
    "discontinuityExpression" -> discExpr,
    "analyticExpression" -> analyticExpr,
    "comparison" -> comparison
  |>
];

Export[
  FileNameJoin[{resultsDir, "compare_oneloop_aux_from_files_result.txt"}],
  ToString[
    <|
      "discontinuityExpression" -> discExpr,
      "analyticExpression" -> analyticExpr,
      "comparison" -> comparison
    |>,
    InputForm
  ],
  "Text"
];
