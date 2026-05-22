(* checks/CompareOneLoopDirectCut.wl
   Legacy comparison for the direct cut-linear attempt.
   Not part of the recommended one-loop workflow. *)

ClearAll["Global`*"];

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

resultsDir = FileNameJoin[{projectDir, "results"}];

resultFile = FileNameJoin[{resultsDir, "oneloop_kernel_direct_cut_result.wl"}];

(* This file exists only if the diagnostic direct-cut target succeeded. *)
If[! FileExistsQ[resultFile],
  Print["Missing result file: ", resultFile];
  Abort[];
];

amflowRules = Get[resultFile];

Print["Raw AMFlow result:"];
Print[amflowRules];

amflowExpr = j[oneloopkerneldirectcut, 1, 1, 1] /. amflowRules;

Print["AMFlow expression:"];
Print[amflowExpr];

(* Same numerical point as in OneLoopKernelDirectCutFamily.wl. *)
xval = 1/3;
m2val = 1/10;
M2val = 1/5;

(* Direct-cut comparison kept for normalization diagnostics only. *)
Aval = xval m2val + (1 - xval) M2val;

(* Analytic result:
     Gamma[1+eps]/eps * Aval^(-eps)
   for 0<x<1.
*)
analyticExpr =
  Normal[Series[Gamma[1 + eps]/eps * Aval^(-eps), {eps, 0, 2}]];

analyticExprN = N[analyticExpr, 30];

Print["Analytic expression:"];
Print[analyticExprN];

(* AMFlow cut normalization may differ.
   Start with 1. If the ratio is constant, use that to infer the convention.
*)
cutNormalization = 1;

normalizedAMFlowExpr = amflowExpr/cutNormalization;

coeffPowers = {-1, 0, 1, 2};

coeffAMFlow[p_] := N[Coefficient[normalizedAMFlowExpr, eps, p], 30];
coeffAnalytic[p_] := N[Coefficient[analyticExprN, eps, p], 30];

comparison =
  Table[
    <|
      "power" -> p,
      "amflow" -> coeffAMFlow[p],
      "analytic" -> coeffAnalytic[p],
      "difference" -> N[coeffAMFlow[p] - coeffAnalytic[p], 30],
      "ratio" -> N[coeffAMFlow[p]/coeffAnalytic[p], 30]
    |>,
    {p, coeffPowers}
  ];

Print["Coefficient comparison:"];
Print[comparison];

Export[
  FileNameJoin[{resultsDir, "compare_oneloop_direct_cut_result.wl"}],
  comparison
];

Export[
  FileNameJoin[{resultsDir, "compare_oneloop_direct_cut_result.txt"}],
  ToString[comparison, InputForm],
  "Text"
];
