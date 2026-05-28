(* targets/CompareKSubloop.wl
   Prepared one-loop AMFlow comparison for the original k-subloop kernel K(q;y).
   Do not run this target unless an AMFlow check is explicitly desired. *)

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "KSubloopUncutFamilies.wl"}]];

SetBasicAMFlowOptions[4];

KSubloopEnvInteger[name_, default_] := Module[
  {value},

  value = Environment[name];
  If[StringQ[value] && StringMatchQ[value, DigitCharacter ..],
    ToExpression[value],
    default
  ]
];

precisionGoal = KSubloopEnvInteger["AMFLOW_PRECISION_GOAL", 10];
epsOrder = KSubloopEnvInteger["AMFLOW_EPS_ORDER", 4];
coeffPowers = {-1, 0, 1};

Xval = 7/10;
yval = 1/4;
q2val = 1;
mk2val = 1/4;
Mk2val = 81/100;

lambdaVal = yval/Xval;
muKlambdaVal = (1 - lambdaVal) mk2val + lambdaVal Mk2val;
CKval = muKlambdaVal - lambdaVal (1 - lambdaVal) q2val;

etaReg = 10^-30;
KclosedExpr =
  N[
    Normal[
      Series[
        Gamma[eps]/Xval * (CKval - I etaReg)^(-eps),
        {eps, 0, 1}
      ]
    ],
    30
  ];


Print["k-subloop GaugeLink comparison"];
Print["This target runs AMFlow if executed."];
Print["Precision goal: ", precisionGoal];
Print["Epsilon order passed to SolveIntegralsGaugeLink: ", epsOrder];
Print["Compared coefficient powers: ", coeffPowers];
Print["Denominators:"];
Print["  D1 = k^2 - mk2"];
Print["  D2 = (k - q)^2 - Mk2"];
Print["  Lplus = y - n.k"];
Print["  Lminus = -(y - n.k)"];
Print["Kinematic point: X=", Xval, ", y=", yval, ", lambda=", lambdaVal, ", q2=", q2val];
Print["Mass squares: mk2=", mk2val, ", Mk2=", Mk2val];
Print["muKlambda = ", muKlambdaVal];
Print["CK = ", CKval, " (Sign=", Sign[CKval], ")"];
Print["Kclosed expression:"];
Print[KclosedExpr];
Print["Kclosed coefficients:"];
Print[
  Table[
    <|
      "power" -> p,
      "coefficient" -> N[Coefficient[KclosedExpr, eps, p], 30]
    |>,
    {p, coeffPowers}
  ]
];

DefineKSubloopUncutPlusFamily[];
plusTarget = j[ksubloopuncutplus, 1, 1, 1];
Print["Running +L GaugeLink integral: ", plusTarget];
plusResult = SolveIntegralsGaugeLink[{plusTarget}, precisionGoal, epsOrder];
Fplus = plusTarget /. plusResult;
Print["+L result:"];
Print[plusResult];

DefineKSubloopUncutMinusFamily[];
minusTarget = j[ksubloopuncutminus, 1, 1, 1];
Print["Running -L GaugeLink integral: ", minusTarget];
minusResult = SolveIntegralsGaugeLink[{minusTarget}, precisionGoal, epsOrder];
Jminus = minusTarget /. minusResult;
Print["-L result:"];
Print[minusResult];

(* Same GaugeLink convention as targets/CompareOneLoop.wl:
   Fplus = integral with 1/(L+i0), Jminus = integral with 1/(-L+i0),
   and 1/(L-i0) = -1/(-L+i0). *)
discExpr = (-Jminus - Fplus)/(2 Pi I);

Print["Discontinuity expression (-Jminus - Fplus)/(2 Pi I):"];
Print[discExpr];

comparison =
  Table[
    With[
      {
        cAM = N[Coefficient[discExpr, eps, p], 30],
        cAN = N[Coefficient[KclosedExpr, eps, p], 30]
      },
      <|
        "power" -> p,
        "amflow" -> cAM,
        "analytic" -> cAN,
        "difference" -> N[cAM - cAN, 30],
        "ratio" -> N[cAM/cAN, 30]
      |>
    ],
    {p, coeffPowers}
  ];

Print["Coefficient comparison:"];
Print[comparison];

fullResult =
  <|
    "metadata" -> <|
      "description" -> "Prepared k-subloop GaugeLink comparison.",
      "precisionGoal" -> precisionGoal,
      "epsOrder" -> epsOrder,
      "coeffPowers" -> coeffPowers,
      "X" -> Xval,
      "y" -> yval,
      "lambda" -> lambdaVal,
      "q2" -> q2val,
      "mk2" -> mk2val,
      "Mk2" -> Mk2val,
      "muKlambda" -> muKlambdaVal,
      "CK" -> CKval,
      "GaugeLinkConvention" -> "discExpr = (-Jminus - Fplus)/(2 Pi I), same as CompareOneLoop.wl"
    |>,
    "plusResult" -> plusResult,
    "minusResult" -> minusResult,
    "discontinuityExpression" -> discExpr,
    "analyticExpression" -> KclosedExpr,
    "comparison" -> comparison
  |>;

Export[
  FileNameJoin[{$ResultsDirectory, "compare_ksubloop_result.wl"}],
  fullResult
];

Export[
  FileNameJoin[{$ResultsDirectory, "compare_ksubloop_result.txt"}],
  ToString[fullResult, InputForm],
  "Text"
];
