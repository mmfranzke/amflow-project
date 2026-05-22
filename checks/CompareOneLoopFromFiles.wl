(* checks/CompareOneLoopFromFiles.wl
   Post-processes previously exported F+ and J- results.
   Use compare-oneloop for a fresh calculation and comparison. *)

ClearAll["Global`*"];

Print["Warning: this comparison reads previously exported +L and -L result files."];
Print["For a fresh one-command calculation and comparison, use ./run.sh compare-oneloop."];

projectDir =
  "/Users/FranzkeMM/Library/Mobile Documents/com~apple~CloudDocs/Documents/Dokumente privat/Studium/ETH Zürich/Master Thesis/amflow-project";

resultsDir = FileNameJoin[{projectDir, "results"}];

(* These must be regenerated after changing kinematics or epsOrder. *)
plusFile = FileNameJoin[{resultsDir, "oneloop_kernel_uncut_plus_result.wl"}];
minusFile = FileNameJoin[{resultsDir, "oneloop_kernel_uncut_minus_result.wl"}];

If[! FileExistsQ[plusFile],
  Print["Missing +L result file: ", plusFile];
  Abort[];
];

If[! FileExistsQ[minusFile],
  Print["Missing -L result file: ", minusFile];
  Abort[];
];

plusRules = Get[plusFile];
minusRules = Get[minusFile];

(* Heads must match the family symbols used when exporting the results. *)
Fplus = j[oneloopkerneluncutplus, 1, 1, 1] /. plusRules;
Jminus = j[oneloopkerneluncutminus, 1, 1, 1] /. minusRules;

(* If Jminus = integral 1/(-L+i0), then integral 1/(L-i0) = -Jminus. *)
discCandidate = (-Jminus - Fplus)/(2 Pi I);

(* Also print sign/flipped variants to diagnose AMFlow conventions. *)
discAltSign = -discCandidate;
discSimpleDiff = (Jminus - Fplus)/(2 Pi I);
discSimpleSum = (Jminus + Fplus)/(2 Pi I);

xval = 1/3;
m2val = 1/10;
M2val = 1/5;
sval = 0;

(* Same routing as CompareOneLoop.wl and the family Numeric rules. *)
Aval = (1 - xval) m2val + xval M2val - xval (1 - xval) sval;

analyticExpr =
  Normal[Series[Gamma[1 + eps]/eps * Aval^(-eps), {eps, 0, 3}]];

analyticExprN = N[analyticExpr, 30];

Print["Fplus = integral with L = n.l + x:"];
Print[Fplus];
Print["Jminus = integral with denominator -L:"];
Print[Jminus];
Print["Analytic expression for current routing and s=0:"];
Print[analyticExprN];

coeffPowers = {-1, 0, 1, 2, 3};

(* Compare one candidate discontinuity against the analytic series. *)
compare[name_, expr_] :=
  <|
    "candidate" -> name,
    "coefficients" ->
      Table[
        With[
          {
            cAM = N[Coefficient[expr, eps, p], 30],
            cAN = N[Coefficient[analyticExprN, eps, p], 30]
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
      ]
  |>;

comparison = {
  compare["(-Jminus - Fplus)/(2 Pi I)", discCandidate],
  compare["(Jminus + Fplus)/(2 Pi I)", discSimpleSum],
  compare["(Jminus - Fplus)/(2 Pi I)", discSimpleDiff],
  compare["-(-Jminus - Fplus)/(2 Pi I)", discAltSign]
};

Print["Discontinuity comparison candidates:"];
Print[comparison];

Export[
  FileNameJoin[{resultsDir, "compare_oneloop_from_files_result.wl"}],
  <|
    "metadata" -> <|
      "description" -> "Post-processing comparison from exported F+ and J- files.",
      "x" -> xval,
      "m2" -> m2val,
      "M2" -> M2val,
      "s" -> sval,
      "A" -> Aval
    |>,
    "comparison" -> comparison
  |>
];

Export[
  FileNameJoin[{resultsDir, "compare_oneloop_from_files_result.txt"}],
  ToString[
    <|
      "metadata" -> <|
        "description" -> "Post-processing comparison from exported F+ and J- files.",
        "x" -> xval,
        "m2" -> m2val,
        "M2" -> M2val,
        "s" -> sval,
        "A" -> Aval
      |>,
      "comparison" -> comparison
    |>,
    InputForm
  ],
  "Text"
];
