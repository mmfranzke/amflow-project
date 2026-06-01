(* families/OuterEq29Families.wl
   GaugeLink families for the Eq. (29) reduced outer integral.
   P means denominator L+i0. M means denominator -L+i0. *)

ClearAll[
  DefineOuterEq29Family,
  DefineOuterEq29PlusFamily,
  DefineOuterEq29MinusFamily
];

DefineOuterEq29Family[familySymbol_, sign_] := Module[{},
  ClearAll[l, p, n, x, y, s, ml2, Ml2, mk2, Mk2];

  AMFlowInfo["Family"] = familySymbol;
  AMFlowInfo["Loop"] = {l};
  AMFlowInfo["Leg"] = {p, n};
  AMFlowInfo["Conservation"] = {};

  AMFlowInfo["Replacement"] = {
    p^2 -> s,
    n^2 -> 0,
    p n -> 1
  };

  (* Eq. (29) denominators:
     D1 = l^2 - ml2
     D2 = (l-p)^2 - Ml2
     D3 = C[l] = muKlambda - a (p-l)^2 - b p^2
     L  = x - n.l

     Here X, lambda, a, b, and muKlambda are written directly in terms of
     x, y, and masses so that AMFlow sees one family with the standard
     replacement rules. *)
  AMFlowInfo["Propagator"] = {
    l^2 - ml2,
    (l - p)^2 - Ml2,
    ((1 - y/(1 - x)) mk2 + (y/(1 - x)) Mk2)
      - (y/(1 - x)) (1 - y/(1 - x)) (p - l)^2
      - (y/(1 - x)) (1 - (y/(1 - x)) x) p^2,
    sign (x - n l)
  };

  AMFlowInfo["Cut"] = {0, 0, 0, 0};
  AMFlowInfo["Prescription"] = {1};

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

  Print["Defined Outer Eq. (29) GaugeLink family: ", familySymbol];
  Print["D1 = l^2 - ml2"];
  Print["D2 = (l-p)^2 - Ml2"];
  Print["D3 = muKlambda - a (p-l)^2 - b p^2"];
  Print["L = ", sign (x - n l)];
  Print["Numeric point: x=3/10, y=1/4, s=0, ml2=49/100, Ml2=4, mk2=1/4, Mk2=81/100"];
];

DefineOuterEq29PlusFamily[] :=
  DefineOuterEq29Family[outereq29plus, 1];

DefineOuterEq29MinusFamily[] :=
  DefineOuterEq29Family[outereq29minus, -1];
