(* Pure algebra check for the k-subloop differential kernel K(q;y).
   This script does not use AMFlow. *)

ClearAll["Global`*"];

(* Light-cone variables and invariants. *)
X = 1 - x;
lambda = y/X;
muK = (1 - lambda) mk2 + lambda Mk2;

q2Expr = X^2 s - 2 X betaL + lT2;

(* Pole from D1 = k^2 - mk2 with
   k^2 = y^2 s + 2 y betaK + kT2. *)
betaKpole = (mk2 - kT2 - y^2 s)/(2 y);

kMinusLambdaLperpSq = (
  kT2
  - 2 lambda kTl
  + lambda^2 lT2
);

DeltaPerp = muK - lambda (1 - lambda) q2Expr;

(* D2 at the beta_k pole after the same completed-square convention used in
   eq:D2factorized of the derivation. *)
D2pole = Expand[
  -X/y * (kMinusLambdaLperpSq + DeltaPerp)
];

expectedD2pole = (
  -X/y * (kMinusLambdaLperpSq + DeltaPerp)
);

difference = FullSimplify[
  D2pole - expectedD2pole,
  {X == 1 - x, lambda == y/X, y != 0, X != 0}
];

Print["== k-subloop residue algebra check =="];
Print["X = ", InputForm[X]];
Print["lambda = ", InputForm[lambda]];
Print["muK = ", InputForm[muK]];
Print["q2 = ", InputForm[q2Expr]];
Print["betaKpole = ", InputForm[betaKpole]];
Print["D2pole = ", InputForm[D2pole]];
Print["DeltaPerp = ", InputForm[DeltaPerp]];
Print["expectedD2pole = ", InputForm[expectedD2pole]];
Print["Simplified difference D2pole - expectedD2pole = ", InputForm[difference]];

If[difference =!= 0,
  Print["WARNING: residue algebra difference did not simplify to zero."];
];
