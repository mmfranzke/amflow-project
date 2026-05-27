(* targets/RunOuterEq29Minus.wl
   Computes the Eq. (29) outer-integral -L+i0 GaugeLink piece.
   This target intentionally uses an epsilon-dependent D3 index. *)

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "OuterEq29Families.wl"}]];

SetBasicAMFlowOptions[4];
DefineOuterEq29MinusFamily[];

OuterEq29EnvInteger[name_, default_] := Module[
  {value},

  value = Environment[name];
  If[StringQ[value] && StringMatchQ[value, DigitCharacter ..],
    ToExpression[value],
    default
  ]
];

precisionGoal = OuterEq29EnvInteger["AMFLOW_PRECISION_GOAL", 10];
epsOrder = OuterEq29EnvInteger["AMFLOW_EPS_ORDER", 4];

targets = {
  j[outereq29minus, 1, 1, eps, 1]
};

Print["Outer Eq. (29) -L+i0 target:"];
Print[targets];
Print["Precision goal: ", precisionGoal];
Print["Epsilon order: ", epsOrder];
Print["This target tests whether AMFlow supports epsilon-dependent index eps on D3."];

result = SolveIntegralsGaugeLink[targets, precisionGoal, epsOrder];

Print["Result:"];
Print[result];

Export[
  FileNameJoin[{$ResultsDirectory, "outer_eq29_minus_result.wl"}],
  result
];

Export[
  FileNameJoin[{$ResultsDirectory, "outer_eq29_minus_result.txt"}],
  ToString[result, InputForm],
  "Text"
];
