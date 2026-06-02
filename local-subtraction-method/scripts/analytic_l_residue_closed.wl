ClearAll["Global`*"];

pointName = Environment["TWO_LOOP_POINT"];
If[! StringQ[pointName] || pointName === "", pointName = "equal_mass_offshell_positive"];

points = <|
  "equal_mass_offshell_positive" -> <|
    "x" -> 1/4, "y" -> 1/4, "pPlus" -> 1, "pMinus" -> 5, "pPerp2" -> 0,
    "ml2" -> 1, "Ml2" -> 1, "mk2" -> 1, "Mk2" -> 1|>,
  "equal_mass_offshell_nonsymmetric" -> <|
    "x" -> 1/5, "y" -> 3/10, "pPlus" -> 1, "pMinus" -> 5, "pPerp2" -> 0,
    "ml2" -> 1, "Ml2" -> 1, "mk2" -> 1, "Mk2" -> 1|>,
  "equal_mass_onshell_nonsymmetric" -> <|
    "x" -> 1/5, "y" -> 3/10, "pPlus" -> 1, "pMinus" -> 0, "pPerp2" -> 0,
    "ml2" -> 1, "Ml2" -> 1, "mk2" -> 1, "Mk2" -> 1|>,
  "equal_mass_onshell_branch" -> <|
    "x" -> 1/4, "y" -> 1/4, "pPlus" -> 1, "pMinus" -> 0, "pPerp2" -> 0,
    "ml2" -> 1, "Ml2" -> 1, "mk2" -> 1, "Mk2" -> 1|>
|>;

If[! KeyExistsQ[points, pointName],
  Print["Unknown point: ", pointName];
  Print["Known points: ", Keys[points]];
  Exit[1];
];

pt = points[pointName];
x = pt["x"]; y = pt["y"]; pPlus = pt["pPlus"]; pMinus = pt["pMinus"]; pPerp2 = pt["pPerp2"];
p2 = pPlus pMinus - pPerp2;
ml2 = pt["ml2"]; Ml2 = pt["Ml2"]; mk2 = pt["mk2"]; Mk2 = pt["Mk2"];

If[pPlus =!= 1 || pPerp2 =!= 0,
  Print["l-residue closed form currently implemented for pPlus=1, pPerp2=0."];
  Exit[1];
];

X = 1 - x;
lambda = y/X;
kappa = lambda (1 - lambda);
muK = (1 - lambda) mk2 + lambda Mk2;
Delta0 = (1 - x) ml2 + x Ml2 - x (1 - x) p2;
A1Short = muK - kappa X p2 + kappa X/x ml2;
extraOffshellTerm = -lambda (1 - lambda x) p2;
A1 = A1Short + extraOffshellTerm;
B1 = kappa/x;
z = 1 - B1 Delta0/A1;
zShort = 1 - B1 Delta0/A1Short;

expr = Gamma[2 eps]/(X eps) (Delta0 A1)^(-eps) Hypergeometric2F1[eps, 1 - eps, 1 + eps, z];
series = Quiet[Normal@Series[expr, {eps, 0, 2}], Series::ztest1];

Print["Analytic l-residue closed candidate. Method12 defaults are unchanged."];
Print["point = ", pointName];
Print["x = ", InputForm[x]];
Print["y = ", InputForm[y]];
Print["X = ", InputForm[X]];
Print["lambda = ", InputForm[lambda]];
Print["kappa = ", InputForm[kappa]];
Print["Delta0 = ", InputForm[Delta0]];
Print["C_full_includes_extra_offshell_term = True"];
Print["extra_offshell_term = ", InputForm[extraOffshellTerm]];
Print["A1_full = ", InputForm[A1]];
Print["A1_short_without_extra_term = ", InputForm[A1Short]];
Print["A1_full - A1_short_without_extra_term = ", InputForm[A1 - A1Short]];
Print["B1 = ", InputForm[B1]];
Print["z = ", InputForm[z]];
Print["z_short_diagnostic = ", InputForm[zShort]];
Print["expected c[-2] = ", InputForm[1/(2 X)]];
Print["series = ", InputForm[series]];

Do[
  coeff = Coefficient[series, eps, power];
  Print["L_RESIDUE_CLOSED_COEFF power=", power, " value=", InputForm[coeff]];
  Print["L_RESIDUE_CLOSED_COEFF_NUMERIC power=", power, " value=", InputForm[N[coeff, 30]]],
  {power, -2, 2}
];
