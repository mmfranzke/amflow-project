(* families/OneLoopKernelDirectCutFamily.wl
   Direct cut-linear diagnostic. The working route uses F+ and J- instead. *)

ClearAll[DefineOneLoopKernelDirectCutFamily];

DefineOneLoopKernelDirectCutFamily[] := Module[{},
  ClearAll[l, p, n, x, s, m2, M2];

  (* Do not use eta or j here. AMFlow reserves both names. *)
  AMFlowInfo["Family"] = oneloopkerneldirectcut;

  AMFlowInfo["Loop"] = {l};

  AMFlowInfo["Leg"] = {p, n};

  AMFlowInfo["Conservation"] = {};

  AMFlowInfo["Replacement"] = {
    p^2 -> s,
    n^2 -> 0,
    p n -> 1
  };

  (* Denominators:
     D1 = l^2 - m^2
     D2 = (l+p)^2 - M^2
     D3 = n.l + x

     D3 is the cut linear denominator representing delta(x+n.l).
  *)
  AMFlowInfo["Propagator"] = {
    l^2 - m2,
    (l + p)^2 - M2,
    n l + x
  };

  (* Mark only the linear denominator as cut. *)
  AMFlowInfo["Cut"] = {
    0, 0, 1
  };

  (* For this cut-linear test, use an insensitive prescription as in AMFlow's
     phase-space examples. *)
  AMFlowInfo["Prescription"] = {
    0
  };

  (* Light-like collinear test point, inside support 0<x<1. *)
  AMFlowInfo["Numeric"] = {
    s -> 0,
    m2 -> 1/10,
    M2 -> 1/5,
    x -> 1/3
  };

  AMFlowInfo["NThread"] = 4;

  Print["Defined one-loop cut-linear collinear kernel family: oneloopkerneldirectcut"];
];
