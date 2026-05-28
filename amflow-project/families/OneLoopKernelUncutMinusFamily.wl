(* families/OneLoopKernelUncutMinusFamily.wl
   GaugeLink family for J- with denominator -L+i0, where L=n.l+x. *)

ClearAll[DefineOneLoopKernelUncutMinusFamily];

DefineOneLoopKernelUncutMinusFamily[] := Module[{},
  ClearAll[l, p, n, x, s, m2, M2];

  (* This family name determines AMFlow cache and result heads. *)
  AMFlowInfo["Family"] = oneloopkerneluncutminus;

  AMFlowInfo["Loop"] = {l};
  AMFlowInfo["Leg"] = {p, n};

  AMFlowInfo["Conservation"] = {};

  AMFlowInfo["Replacement"] = {
    p^2 -> s,
    n^2 -> 0,
    p n -> 1
  };

  (* The denominator is -L. This realizes 1/(L-i0) as -1/(-L+i0). *)
  AMFlowInfo["Propagator"] = {
    l^2 - m2,
    (l + p)^2 - M2,
    -(n l + x)
  };

  (* Uncut linear denominator: this is J-, not the delta constraint. *)
  AMFlowInfo["Cut"] = {
    0, 0, 0
  };

  AMFlowInfo["Prescription"] = {
    1
  };

  (* Same light-like point as the analytic one-loop comparison. *)
  AMFlowInfo["Numeric"] = {
    s -> 0,
    m2 -> 1/10,
    M2 -> 1/5,
    x -> 1/3
  };

  AMFlowInfo["NThread"] = 4;

  Print["Defined uncut minus-linear-propagator family: oneloopkerneluncutminus"];
];
