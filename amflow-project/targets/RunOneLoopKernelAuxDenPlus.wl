(* targets/RunOneLoopKernelAuxDenPlus.wl
   Computes the one-loop F+ test with an auxiliary denominator at exponent 0. *)

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "OneLoopKernelAuxDenFamilies.wl"}]];

SetBasicAMFlowOptions[4];
DefineOneLoopKernelAuxDenPlusFamily[];

targets = {
  j[oneloopkernelauxplus, 1, 1, 1, 0]
};

precisionGoal = 20;
epsOrder = 5;

Print["Targets: ", targets];
Print["Precision goal: ", precisionGoal];
Print["Epsilon order: ", epsOrder];
Print["Note: D4 is auxiliary and has exponent 0."];

result = SolveIntegralsGaugeLink[targets, precisionGoal, epsOrder];

Print["Result:"];
Print[result];

Export[
  FileNameJoin[{$ResultsDirectory, "oneloop_kernel_aux_plus_result.wl"}],
  result
];

Export[
  FileNameJoin[{$ResultsDirectory, "oneloop_kernel_aux_plus_result.txt"}],
  ToString[result, InputForm],
  "Text"
];
