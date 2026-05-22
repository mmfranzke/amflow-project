(* targets/RunOneLoopKernelUncutMinus.wl
   Computes J- with the uncut -L+i0 linear denominator.
   This is one piece of the discontinuity reconstruction, not the delta kernel. *)

Get["/Users/FranzkeMM/Library/Mobile Documents/com~apple~CloudDocs/Documents/Dokumente privat/Studium/ETH Zürich/Master Thesis/amflow-project/config/LoadAMFlow.wl"];
Get[FileNameJoin[{$ProjectDirectory, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{$ProjectDirectory, "families", "OneLoopKernelUncutMinusFamily.wl"}]];

SetBasicAMFlowOptions[4];
DefineOneLoopKernelUncutMinusFamily[];

targets = {
  j[oneloopkerneluncutminus, 1, 1, 1]
};

(* GaugeLink order convention differs from plain SolveIntegrals. *)
precisionGoal = 20;
epsOrder = 5;

Print["Targets: ", targets];
Print["Precision goal: ", precisionGoal];
Print["Epsilon order: ", epsOrder];
Print["Note: J- is uncut. Combine it with F+ to reconstruct the delta constraint."];

(* Use GaugeLink because ordinary SolveIntegrals does not handle -L directly. *)
result = SolveIntegralsGaugeLink[targets, precisionGoal, epsOrder];

Print["Result:"];
Print[result];

(* These files are consumed only by compare-oneloop-from-files. *)
Export[
  FileNameJoin[{$ResultsDirectory, "oneloop_kernel_uncut_minus_result.wl"}],
  result
];

Export[
  FileNameJoin[{$ResultsDirectory, "oneloop_kernel_uncut_minus_result.txt"}],
  ToString[result, InputForm],
  "Text"
];
