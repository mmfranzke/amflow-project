(* targets/RunTwoLoopKernelUncutPM.wl
   Computes the PM uncut two-loop GaugeLink piece alone. *)

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "TwoLoopKernelUncutFamilies.wl"}]];
Get[FileNameJoin[{projectDir, "targets", "TwoLoopKernelTargetTools.wl"}]];

SetBasicAMFlowOptions[4];

(* Override with AMFLOW_PRECISION_GOAL and AMFLOW_EPS_ORDER for quick tests. *)
SolveTwoLoopKernelPiece["PM", TwoLoopPrecisionGoal[10], TwoLoopEpsOrder[4]];
