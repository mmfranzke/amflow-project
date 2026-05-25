(* families/OneLoopKernelAuxDenFamilies.wl
   One-loop GaugeLink test with one extra denominator at exponent 0.
   This checks whether a zero-power auxiliary denominator changes the result. *)

ClearAll[
  DefineOneLoopKernelAuxDenFamily,
  DefineOneLoopKernelAuxDenPlusFamily,
  DefineOneLoopKernelAuxDenMinusFamily
];

DefineOneLoopKernelAuxDenFamily[familySymbol_, sign_] := Module[{},
  ClearAll[l, p, n, x, s, m2, M2];

  AMFlowInfo["Family"] = familySymbol;

  AMFlowInfo["Loop"] = {l};
  AMFlowInfo["Leg"] = {p, n};

  AMFlowInfo["Conservation"] = {};

  AMFlowInfo["Replacement"] = {
    p^2 -> s,
    n^2 -> 0,
    p n -> 1
  };

  Print[
    "Unsupported check: AMFlow's GaugeLink completeness check counts the ",
    "zero-power auxiliary denominator after linear-denominator squaring. ",
    "For this one-loop family there are only three independent scalar products, ",
    "so the four-denominator auxiliary test aborts before solving."
  ];
  Print[
    "Use the ordinary uncut plus/minus targets instead: ",
    "oneloop-kernel-uncut-plus, oneloop-kernel-uncut-minus, or compare-oneloop."
  ];
  Abort[];

  (* Kept below as documentation of the attempted unsupported check.
     D4 is auxiliary and would be called with exponent 0:
     j[family, 1, 1, 1, 0].
     AMFlow still includes it in CheckCompleteness. *)
  AMFlowInfo["Propagator"] = {
    l^2 - m2,
    (l + p)^2 - M2,
    sign (n l + x),
    (l - p)^2
  };

  AMFlowInfo["Cut"] = {
    0, 0, 0, 0
  };

  AMFlowInfo["Prescription"] = {
    1
  };

  (* Same point as the ordinary one-loop GaugeLink comparison. *)
  AMFlowInfo["Numeric"] = {
    s -> 0,
    m2 -> 1/10,
    M2 -> 1/5,
    x -> 1/3
  };

  AMFlowInfo["NThread"] = 4;

  Print[
    "Defined one-loop auxiliary-denominator family: ",
    familySymbol,
    " with sign=", sign
  ];
];

DefineOneLoopKernelAuxDenPlusFamily[] :=
  DefineOneLoopKernelAuxDenFamily[oneloopkernelauxplus, 1];

DefineOneLoopKernelAuxDenMinusFamily[] :=
  DefineOneLoopKernelAuxDenFamily[oneloopkernelauxminus, -1];
