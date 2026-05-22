(* targets/RunTwoLoopKernelUncutPP.wl
   Computes the PP uncut two-loop GaugeLink piece alone. *)

Get["/Users/FranzkeMM/Library/Mobile Documents/com~apple~CloudDocs/Documents/Dokumente privat/Studium/ETH Zürich/Master Thesis/amflow-project/config/LoadAMFlow.wl"];
Get[FileNameJoin[{$ProjectDirectory, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{$ProjectDirectory, "families", "TwoLoopKernelUncutFamilies.wl"}]];
Get[FileNameJoin[{$ProjectDirectory, "targets", "TwoLoopKernelTargetTools.wl"}]];

SetBasicAMFlowOptions[4];

(* epsOrder = 4 should include the finite coefficient for this two-loop test. *)
SolveTwoLoopKernelPiece["PP", 10, 4];
