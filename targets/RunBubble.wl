(* targets/RunBubble.wl
   Ordinary bubble sanity check for the AMFlow installation. *)

Get["/Users/FranzkeMM/Library/Mobile Documents/com~apple~CloudDocs/Documents/Dokumente privat/Studium/ETH Zürich/Master Thesis/amflow-project/config/LoadAMFlow.wl"];
Get[FileNameJoin[{$ProjectDirectory, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{$ProjectDirectory, "families", "BubbleFamily.wl"}]];

SetBasicAMFlowOptions[4];
DefineBubbleFamily[];

targets = {
  j[bubble, 1, 1]
};

(* epsOrder = 3 is enough for the setup check and finite terms. *)
precisionGoal = 20;
epsOrder = 4;

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
