(* targets/RunTwoLoopKernelUncutPM.wl
   Computes the PM uncut two-loop GaugeLink piece alone. *)

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "TwoLoopKernelUncutFamilies.wl"}]];
Get[FileNameJoin[{projectDir, "targets", "TwoLoopKernelTargetTools.wl"}]];

SetBasicAMFlowOptions[4];

(* epsOrder = 4 should include the finite coefficient for this two-loop test. *)
SolveTwoLoopKernelPiece["PM", 10, 4];
