(* targets/RunOneLoopKernelAuxDenMinus.wl
   Computes the one-loop J- test with an auxiliary denominator at exponent 0. *)

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "OneLoopKernelAuxDenFamilies.wl"}]];

SetBasicAMFlowOptions[4];
DefineOneLoopKernelAuxDenMinusFamily[];

targets = {
  j[oneloopkernelauxminus, 1, 1, 1, 0]
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
  FileNameJoin[{$ResultsDirectory, "oneloop_kernel_aux_minus_result.wl"}],
  result
];

Export[
  FileNameJoin[{$ResultsDirectory, "oneloop_kernel_aux_minus_result.txt"}],
  ToString[result, InputForm],
  "Text"
];
