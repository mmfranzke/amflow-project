(* config/LoadAMFlow.wl
   Shared project bootstrap: loads local paths, AMFlow, and output folders. *)

ClearAll["Global`*"];

(* Keep this explicit. The iCloud path contains spaces and unicode. *)
$ProjectDirectory =
  "/Users/FranzkeMM/Library/Mobile Documents/com~apple~CloudDocs/Documents/Dokumente privat/Studium/ETH Zürich/Master Thesis/amflow-project";

(* AMFlow itself stays outside iCloud. *)
$AMFlowDirectory =
  FileNameJoin[{$HomeDirectory, "physics", "tools", "AMFlow"}];

If[! DirectoryQ[$ProjectDirectory],
  Print["Project directory not found: ", $ProjectDirectory];
  Abort[];
];

If[! DirectoryQ[$AMFlowDirectory],
  Print["AMFlow directory not found: ", $AMFlowDirectory];
  Abort[];
];

AppendTo[$Path, $AMFlowDirectory];

(* Adds FiniteFlow, LiteIBP, LiteRed, and dynamic-library paths. *)
Get["/Users/FranzkeMM/physics/tools/amflow-paths.wl"];

Get["AMFlow.m"];

Print["Loaded AMFlow from: ", $AMFlowDirectory];
Print["Project directory: ", $ProjectDirectory];

$ResultsDirectory = FileNameJoin[{$ProjectDirectory, "results"}];
$LogsDirectory = FileNameJoin[{$ProjectDirectory, "logs"}];

(* Targets export machine-readable results and text snapshots here. *)
If[! DirectoryQ[$ResultsDirectory], CreateDirectory[$ResultsDirectory]];
If[! DirectoryQ[$LogsDirectory], CreateDirectory[$LogsDirectory]];
