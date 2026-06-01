(* families/BubbleFamily.wl
   Ordinary massive bubble used as the first setup sanity check. *)

ClearAll[DefineBubbleFamily];

DefineBubbleFamily[] := Module[{},
  ClearAll[l, p, s, m2];

  AMFlowInfo["Family"] = bubble;
  AMFlowInfo["Loop"] = {l};
  AMFlowInfo["Leg"] = {p};

  AMFlowInfo["Conservation"] = {};

  AMFlowInfo["Replacement"] = {
    p^2 -> s
  };

  (* Ordinary quadratic denominators. No linear propagators here. *)
  AMFlowInfo["Propagator"] = {
    l^2 - m2,
    (l + p)^2 - m2
  };

  (* Euclidean point chosen only to verify the AMFlow setup. *)
  AMFlowInfo["Numeric"] = {
    s -> -1,
    m2 -> 1/10
  };

  AMFlowInfo["NThread"] = 4;

  Print["Defined family: bubble"];
];
