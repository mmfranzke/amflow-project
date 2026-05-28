(* config/LoadAMFlow.wl
   Shared project bootstrap: loads local paths, AMFlow, and output folders. *)

ClearAll[projectDirFromInput];

projectDirFromInput =
  DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

$ProjectDirectory = projectDirFromInput;

$AMFlowDirectory =
  If[StringQ[Environment["AMFLOW_DIR"]] && Environment["AMFLOW_DIR"] =!= "",
    ExpandFileName[Environment["AMFLOW_DIR"]],
    FileNameJoin[{$HomeDirectory, "physics", "tools", "AMFlow"}]
  ];

$AMFlowPathsFile =
  If[StringQ[Environment["AMFLOW_PATHS_FILE"]] && Environment["AMFLOW_PATHS_FILE"] =!= "",
    ExpandFileName[Environment["AMFLOW_PATHS_FILE"]],
    FileNameJoin[{$HomeDirectory, "physics", "tools", "amflow-paths.wl"}]
  ];

If[! DirectoryQ[$ProjectDirectory],
  Print["Project directory not found: ", $ProjectDirectory];
  Abort[];
];

If[! DirectoryQ[$AMFlowDirectory],
  Print["AMFlow directory not found: ", $AMFlowDirectory];
  Abort[];
];

If[! FileExistsQ[$AMFlowPathsFile],
  Print["AMFlow paths file not found: ", $AMFlowPathsFile];
  Print["Set AMFLOW_PATHS_FILE if it lives somewhere else."];
  Abort[];
];

If[! MemberQ[$Path, $AMFlowDirectory],
  AppendTo[$Path, $AMFlowDirectory]
];

(* Adds FiniteFlow, LiteIBP, LiteRed, and dynamic-library paths. *)
Get[$AMFlowPathsFile];

Print["Dependency check before loading AMFlow:"];
Print["  FiniteFlow.m: ", FindFile["FiniteFlow.m"]];
Print["  LiteIBP.m: ", FindFile["LiteIBP.m"]];
Print["  LiteRed.m: ", FindFile["LiteRed.m"]];
Print["  fflowmlink.so: ", FileNames["fflowmlink.so", $LibraryPath]];
Print["  fflowmlink.dylib: ", FileNames["fflowmlink.dylib", $LibraryPath]];

Get["AMFlow.m"];

(* Re-apply after AMFlow loading in case AMFlow changed $Path/$LibraryPath. *)
Get[$AMFlowPathsFile];

Print["Dependency check after loading AMFlow:"];
Print["  FiniteFlow.m: ", FindFile["FiniteFlow.m"]];
Print["  LiteIBP.m: ", FindFile["LiteIBP.m"]];
Print["  LiteRed.m: ", FindFile["LiteRed.m"]];
Print["  fflowmlink.so: ", FileNames["fflowmlink.so", $LibraryPath]];
Print["  fflowmlink.dylib: ", FileNames["fflowmlink.dylib", $LibraryPath]];

Print["Loaded AMFlow from: ", $AMFlowDirectory];
Print["Project directory: ", $ProjectDirectory];

$ResultsDirectory = FileNameJoin[{$ProjectDirectory, "results"}];
$LogsDirectory = FileNameJoin[{$ProjectDirectory, "logs"}];

If[! DirectoryQ[$ResultsDirectory], CreateDirectory[$ResultsDirectory]];
If[! DirectoryQ[$LogsDirectory], CreateDirectory[$LogsDirectory]];