(* config/TwoLoopKinematicPoints.wl
   Shared point definitions for the two-loop double-collinear checks. *)

ClearAll[
  TwoLoopKinematicPoints,
  TwoLoopKinematicPointNames,
  TwoLoopKinematicPoint,
  TwoLoopSelectedPointName,
  TwoLoopSelectedPoint,
  TwoLoopPointRules,
  TwoLoopPrintPoint
];

TwoLoopKinematicPoints[] := <|
  "branch_safe_rational" -> <|
    "name" -> "branch_safe_rational",
    "pPlus" -> 1,
    "pMinus" -> 0,
    "pPerp2" -> 0,
    "p2" -> 0,
    "x" -> 3/10,
    "y" -> 1/4,
    "ml2" -> 49/100,
    "Ml2" -> 4,
    "mk2" -> 1/4,
    "Mk2" -> 81/100
  |>,
  "equal_mass_onshell_branch" -> <|
    "name" -> "equal_mass_onshell_branch",
    "pPlus" -> 1,
    "pMinus" -> 0,
    "pPerp2" -> 0,
    "p2" -> 0,
    "x" -> 1/4,
    "y" -> 1/4,
    "ml2" -> 1,
    "Ml2" -> 1,
    "mk2" -> 1,
    "Mk2" -> 1
  |>,
  "equal_mass_offshell_positive" -> <|
    "name" -> "equal_mass_offshell_positive",
    "pPlus" -> 1,
    "pMinus" -> 5,
    "pPerp2" -> 0,
    "p2" -> 5,
    "x" -> 1/4,
    "y" -> 1/4,
    "ml2" -> 1,
    "Ml2" -> 1,
    "mk2" -> 1,
    "Mk2" -> 1
  |>,
  "equal_mass_offshell_nonsymmetric" -> <|
    "name" -> "equal_mass_offshell_nonsymmetric",
    "pPlus" -> 1,
    "pMinus" -> 5,
    "pPerp2" -> 0,
    "p2" -> 5,
    "x" -> 1/5,
    "y" -> 3/10,
    "ml2" -> 1,
    "Ml2" -> 1,
    "mk2" -> 1,
    "Mk2" -> 1
  |>,
  "equal_mass_onshell_nonsymmetric" -> <|
    "name" -> "equal_mass_onshell_nonsymmetric",
    "pPlus" -> 1,
    "pMinus" -> 0,
    "pPerp2" -> 0,
    "p2" -> 0,
    "x" -> 1/5,
    "y" -> 3/10,
    "ml2" -> 1,
    "Ml2" -> 1,
    "mk2" -> 1,
    "Mk2" -> 1
  |>,
  "four_mass_offshell_positive_A" -> <|
    "name" -> "four_mass_offshell_positive_A",
    "pPlus" -> 1,
    "pMinus" -> 5,
    "pPerp2" -> 0,
    "p2" -> 5,
    "x" -> 1/5,
    "y" -> 3/10,
    "ml2" -> 49/100,
    "Ml2" -> 4,
    "mk2" -> 1/4,
    "Mk2" -> 81/100
  |>,
  "four_mass_onshell_nonsymmetric_A" -> <|
    "name" -> "four_mass_onshell_nonsymmetric_A",
    "pPlus" -> 1,
    "pMinus" -> 0,
    "pPerp2" -> 0,
    "p2" -> 0,
    "x" -> 1/5,
    "y" -> 3/10,
    "ml2" -> 49/100,
    "Ml2" -> 4,
    "mk2" -> 1/4,
    "Mk2" -> 81/100
  |>
|>;

TwoLoopKinematicPointNames[] := Keys[TwoLoopKinematicPoints[]];

TwoLoopKinematicPoint[name_String] := Module[{points = TwoLoopKinematicPoints[]},
  If[! KeyExistsQ[points, name],
    Print["Unknown two-loop point: ", name];
    Print["Known points: ", TwoLoopKinematicPointNames[]];
    Abort[];
  ];
  points[name]
];

TwoLoopSelectedPointName[default_: Automatic] := Module[{env, fallback},
  fallback =
    If[default === Automatic,
      If[StringQ[$TwoLoopDefaultPointName], $TwoLoopDefaultPointName, "branch_safe_rational"],
      default
    ];
  env = Environment["TWO_LOOP_POINT"];
  If[StringQ[env] && StringLength[env] > 0, env, fallback]
];

TwoLoopSelectedPoint[default_: Automatic] :=
  TwoLoopKinematicPoint[TwoLoopSelectedPointName[default]];

TwoLoopPointRules[point_Association] := {
  s -> point["p2"],
  ml2 -> point["ml2"],
  Ml2 -> point["Ml2"],
  mk2 -> point["mk2"],
  Mk2 -> point["Mk2"],
  x -> point["x"],
  y -> point["y"]
};

TwoLoopPrintPoint[point_Association] := Module[{},
  Print[
    "Two-loop point: ", point["name"],
    ", pPlus=", point["pPlus"],
    ", pMinus=", point["pMinus"],
    ", pPerp2=", point["pPerp2"],
    ", p2=", point["p2"],
    ", x=", point["x"],
    ", y=", point["y"]
  ];
  Print[
    "Mass squares: ml2=", point["ml2"],
    ", Ml2=", point["Ml2"],
    ", mk2=", point["mk2"],
    ", Mk2=", point["Mk2"]
  ];
];
