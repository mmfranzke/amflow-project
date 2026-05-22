(* targets/RunOneLoopKernelUncutPlus.wl
   Computes F+ with the uncut L+i0 linear denominator.
   This is one piece of the discontinuity reconstruction, not the delta kernel. *)

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "OneLoopKernelUncutPlusFamily.wl"}]];

SetBasicAMFlowOptions[4];
DefineOneLoopKernelUncutPlusFamily[];

targets = {
  j[oneloopkerneluncutplus, 1, 1, 1]
};

(* GaugeLink order convention differs from plain SolveIntegrals. *)
precisionGoal = 20;
epsOrder = 5;

Print["Targets: ", targets];
Print["Precision goal: ", precisionGoal];
Print["Epsilon order: ", epsOrder];
Print["Note: F+ is uncut. Combine it with J- to reconstruct the delta constraint."];

(* Use GaugeLink because ordinary SolveIntegrals does not handle L directly. *)
result = SolveIntegralsGaugeLink[targets, precisionGoal, epsOrder];

Print["Result:"];
Print[result];

(* These files are consumed only by compare-oneloop-from-files. *)
Export[
  FileNameJoin[{$ResultsDirectory, "oneloop_kernel_uncut_plus_result.wl"}],
  result
];

Export[
  FileNameJoin[{$ResultsDirectory, "oneloop_kernel_uncut_plus_result.txt"}],
  ToString[result, InputForm],
  "Text"
];
