(* families/TwoLoopKernelUncutFamilies.wl
   Four uncut GaugeLink families for the two-loop double-delta kernel.
   P means denominator L+i0. M means denominator -L+i0. *)

ClearAll[
  DefineTwoLoopKernelUncutFamily,
  DefineTwoLoopKernelUncutPPFamily,
  DefineTwoLoopKernelUncutPMFamily,
  DefineTwoLoopKernelUncutMPFamily,
  DefineTwoLoopKernelUncutMMFamily
];

DefineTwoLoopKernelUncutFamily[familySymbol_, sx_, sy_] := Module[{},
  ClearAll[l, k, p, n, x, y, s, ml2, Ml2, mk2, Mk2];

  (* Do not use eta or j here. AMFlow reserves both names. *)
  AMFlowInfo["Family"] = familySymbol;

  AMFlowInfo["Loop"] = {l, k};
  AMFlowInfo["Leg"] = {p, n};

  AMFlowInfo["Conservation"] = {};

  AMFlowInfo["Replacement"] = {
    p^2 -> s,
    n^2 -> 0,
    p n -> 1
  };

  (* Routing from derivations/snippets/dckernel/main.tex:
     D1 = l^2 - ml^2
     D2 = (l-p)^2 - Ml^2
     D3 = k^2 - mk^2
     D4 = (l+k-p)^2 - Mk^2
     Lx = x - n.l
     Ly = y - n.k

     D7 = (k-p)^2 is an auxiliary denominator. It completes the
     two-loop scalar-product basis and is used with exponent 0.
  *)
  AMFlowInfo["Propagator"] = {
    l^2 - ml2,
    (l - p)^2 - Ml2,
    k^2 - mk2,
    (l + k - p)^2 - Mk2,
    sx (x - n l),
    sy (y - n k),
    (k - p)^2
  };

  (* These are uncut GaugeLink denominators. The delta constraints are
     reconstructed later by a double discontinuity. *)
  AMFlowInfo["Cut"] = {
    0, 0, 0, 0, 0, 0, 0
  };

  (* Both loop momenta use +i0. Negative L denominators implement L-i0
     through 1/(L-i0) = -1/(-L+i0). *)
  AMFlowInfo["Prescription"] = {
    1, 1
  };

  (* Light-like, phase-safe point for debugging the discontinuity normalization.
     Ml2 is deliberately large enough that Delta0 and Delta1 are both positive. *)
  AMFlowInfo["Numeric"] = {
    s -> 0,
    ml2 -> 49/100,
    Ml2 -> 4,
    mk2 -> 1/4,
    Mk2 -> 81/100,
    x -> 3/10,
    y -> 1/4
  };

  AMFlowInfo["NThread"] = 4;

  Print[
    "Defined two-loop uncut GaugeLink family: ",
    familySymbol,
    " with signs sx=", sx,
    ", sy=", sy
  ];
];

DefineTwoLoopKernelUncutPPFamily[] :=
  DefineTwoLoopKernelUncutFamily[twoloopkerneluncutpp, 1, 1];

DefineTwoLoopKernelUncutPMFamily[] :=
  DefineTwoLoopKernelUncutFamily[twoloopkerneluncutpm, 1, -1];

DefineTwoLoopKernelUncutMPFamily[] :=
  DefineTwoLoopKernelUncutFamily[twoloopkerneluncutmp, -1, 1];

DefineTwoLoopKernelUncutMMFamily[] :=
  DefineTwoLoopKernelUncutFamily[twoloopkerneluncutmm, -1, -1];
