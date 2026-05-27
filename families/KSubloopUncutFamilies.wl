(* families/KSubloopUncutFamilies.wl
   GaugeLink families for the original k-subloop kernel K(q;y).
   P means denominator L+i0. M means denominator -L+i0. *)

ClearAll[
  DefineKSubloopUncutFamily,
  DefineKSubloopUncutPlusFamily,
  DefineKSubloopUncutMinusFamily
];

DefineKSubloopUncutFamily[familySymbol_, sign_] := Module[{},
  ClearAll[k, q, n, X, y, q2, mk2, Mk2];

  AMFlowInfo["Family"] = familySymbol;
  AMFlowInfo["Loop"] = {k};
  AMFlowInfo["Leg"] = {q, n};
  AMFlowInfo["Conservation"] = {};

  AMFlowInfo["Replacement"] = {
    q^2 -> q2,
    n^2 -> 0,
    q n -> X
  };

  (* D1 = k^2 - mk2, D2 = (k-q)^2 - Mk2, L = y - n.k. *)
  AMFlowInfo["Propagator"] = {
    k^2 - mk2,
    (k - q)^2 - Mk2,
    sign (y - n k)
  };

  AMFlowInfo["Cut"] = {0, 0, 0};
  AMFlowInfo["Prescription"] = {1};

  AMFlowInfo["Numeric"] = {
    X -> 7/10,
    y -> 1/4,
    q2 -> 1,
    mk2 -> 1/4,
    Mk2 -> 81/100
  };

  AMFlowInfo["NThread"] = 4;

  Print["Defined k-subloop uncut GaugeLink family: ", familySymbol];
  Print["Denominators: D1=k^2-mk2, D2=(k-q)^2-Mk2, L=", sign (y - n k)];
  Print["Numeric point: X=7/10, y=1/4, q2=1, mk2=1/4, Mk2=81/100"];
];

DefineKSubloopUncutPlusFamily[] :=
  DefineKSubloopUncutFamily[ksubloopuncutplus, 1];

DefineKSubloopUncutMinusFamily[] :=
  DefineKSubloopUncutFamily[ksubloopuncutminus, -1];
