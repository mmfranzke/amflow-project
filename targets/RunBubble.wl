(* targets/RunBubble.wl
   Ordinary bubble sanity check for the AMFlow installation. *)

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "BubbleFamily.wl"}]];

SetBasicAMFlowOptions[4];
DefineBubbleFamily[];

targets = {
  j[bubble, 1, 1]
};

(* epsOrder = 3 is enough for the setup check and finite terms. *)
precisionGoal = 20;
epsOrder = 3;

Print["Targets: ", targets];
Print["Precision goal: ", precisionGoal];
Print["Epsilon order: ", epsOrder];

result = SolveIntegrals[targets, precisionGoal, epsOrder];

Print["Result:"];
Print[result];

(* Export both a reloadable Mathematica expression and a plain text snapshot. *)
Export[
  FileNameJoin[{$ResultsDirectory, "bubble_result.wl"}],
  result
];

Export[
  FileNameJoin[{$ResultsDirectory, "bubble_result.txt"}],
  ToString[result, InputForm],
  "Text"
];
