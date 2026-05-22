(* targets/RunTwoLoopKernelUncutPP.wl
   Computes the PP uncut two-loop GaugeLink piece alone. *)

projectDir = DirectoryName[DirectoryName[$InputFileName]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "TwoLoopKernelUncutFamilies.wl"}]];
Get[FileNameJoin[{projectDir, "targets", "TwoLoopKernelTargetTools.wl"}]];

SetBasicAMFlowOptions[4];

(* epsOrder = 4 should include the finite coefficient for this two-loop test. *)
SolveTwoLoopKernelPiece["PP", 10, 4];
