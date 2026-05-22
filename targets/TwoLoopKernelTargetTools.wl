(* targets/TwoLoopKernelTargetTools.wl
   Shared helpers for the two-loop GaugeLink targets and comparisons. *)

ClearAll[
  TwoLoopKernelPieces,
  DefineTwoLoopKernelPiece,
  TwoLoopKernelFamilySymbol,
  TwoLoopKernelTarget,
  TwoLoopKernelResultBase,
  SolveTwoLoopKernelPiece
];

TwoLoopKernelPieces[] := {"PP", "PM", "MP", "MM"};

DefineTwoLoopKernelPiece["PP"] := DefineTwoLoopKernelUncutPPFamily[];
DefineTwoLoopKernelPiece["PM"] := DefineTwoLoopKernelUncutPMFamily[];
DefineTwoLoopKernelPiece["MP"] := DefineTwoLoopKernelUncutMPFamily[];
DefineTwoLoopKernelPiece["MM"] := DefineTwoLoopKernelUncutMMFamily[];

TwoLoopKernelFamilySymbol["PP"] := twoloopkerneluncutpp;
TwoLoopKernelFamilySymbol["PM"] := twoloopkerneluncutpm;
TwoLoopKernelFamilySymbol["MP"] := twoloopkerneluncutmp;
TwoLoopKernelFamilySymbol["MM"] := twoloopkerneluncutmm;

TwoLoopKernelTarget[piece_] :=
  j[TwoLoopKernelFamilySymbol[piece], 1, 1, 1, 1, 1, 1, 0];

TwoLoopKernelResultBase[piece_] :=
  "twoloop_kernel_uncut_" <> ToLowerCase[piece] <> "_result";

SolveTwoLoopKernelPiece[piece_, precisionGoal_, epsOrder_] := Module[
  {target, result, base},

  If[! MemberQ[TwoLoopKernelPieces[], piece],
    Print["Unknown two-loop piece: ", piece];
    Abort[];
  ];

  DefineTwoLoopKernelPiece[piece];
  target = TwoLoopKernelTarget[piece];

  Print["Running two-loop GaugeLink piece ", piece, ": ", target];
  Print["Precision goal: ", precisionGoal];
  Print["Epsilon order: ", epsOrder];

  (* First two-loop pass: keep the order modest until reduction is stable. *)
  result = SolveIntegralsGaugeLink[{target}, precisionGoal, epsOrder];

  Print["Result for piece ", piece, ":"];
  Print[result];

  base = TwoLoopKernelResultBase[piece];

  (* Split targets and compare-twoloop-from-files read these exports. *)
  Export[
    FileNameJoin[{$ResultsDirectory, base <> ".wl"}],
    result
  ];

  Export[
    FileNameJoin[{$ResultsDirectory, base <> ".txt"}],
    ToString[result, InputForm],
    "Text"
  ];

  result
];
