(* targets/RunOneLoopKernelDirectCut.wl
   Diagnostic direct cut-linear attempt.
   Kept to document the route that currently fails in AMFlow. *)

projectDir = DirectoryName[DirectoryName[$InputFileName]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "OneLoopKernelDirectCutFamily.wl"}]];

SetBasicAMFlowOptions[4];
DefineOneLoopKernelDirectCutFamily[];

targets = {
  j[oneloopkerneldirectcut, 1, 1, 1]
};

(* Kept low because this target is diagnostic and currently expected to fail. *)
precisionGoal = 20;
epsOrder = 4;

Print["Targets: ", targets];
Print["Precision goal: ", precisionGoal];
Print["Epsilon order: ", epsOrder];

(* GaugeLink is needed because the family contains a linear denominator. *)
result = SolveIntegralsGaugeLink[targets, precisionGoal, epsOrder];

Print["Result:"];
Print[result];

Export[
  FileNameJoin[{$ResultsDirectory, "oneloop_kernel_direct_cut_result.wl"}],
  result
];

Export[
  FileNameJoin[{$ResultsDirectory, "oneloop_kernel_direct_cut_result.txt"}],
  ToString[result, InputForm],
  "Text"
];
