(* config/AMFlowOptions.wl
   Small shared option helper used by all target scripts. *)

ClearAll[SetBasicAMFlowOptions];

SetBasicAMFlowOptions[nThreads_: 4] := Module[{},
  (* AMFlow reads this global option when launching reduction jobs. *)
  AMFlowInfo["NThread"] = nThreads;
  Print["AMFlow threads: ", nThreads];
];
