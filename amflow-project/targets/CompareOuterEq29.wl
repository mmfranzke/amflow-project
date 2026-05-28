(* targets/CompareOuterEq29.wl
   Direct AMFlow/GaugeLink check of Eq. (29) against the analytic Eq. (51)/(53).
   This target runs AMFlow if executed. *)

projectDir = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];

Get[FileNameJoin[{projectDir, "config", "LoadAMFlow.wl"}]];
Get[FileNameJoin[{projectDir, "config", "AMFlowOptions.wl"}]];
Get[FileNameJoin[{projectDir, "families", "OuterEq29Families.wl"}]];

SetBasicAMFlowOptions[4];

OuterEq29EnvInteger[name_, default_] := Module[
  {value},

  value = Environment[name];
  If[StringQ[value] && StringMatchQ[value, DigitCharacter ..],
    ToExpression[value],
    default
  ]
];

precisionGoal = OuterEq29EnvInteger["AMFLOW_PRECISION_GOAL", 10];
epsOrder = OuterEq29EnvInteger["AMFLOW_EPS_ORDER", 4];
coeffPowers = {-2, -1, 0};

xval = 3/10;
yval = 1/4;
sval = 0;
ml2val = 49/100;
Ml2val = 4;
mk2val = 1/4;
Mk2val = 81/100;

Xval = 1 - xval;
lambdaVal = yval/Xval;
aVal = lambdaVal (1 - lambdaVal);
bVal = lambdaVal (1 - lambdaVal xval);
muLxVal = (1 - xval) ml2val + xval Ml2val;
muKlambdaVal = (1 - lambdaVal) mk2val + lambdaVal Mk2val;
Delta0Val = muLxVal - xval (1 - xval) sval;
DeltaWVal = (
  ((1 + aVal) xval (1 - xval) + bVal) sval
  - muKlambdaVal
  + aVal Ml2val
  - (1 + aVal) muLxVal
);
Delta1Val = (1 + aVal) Delta0Val + DeltaWVal;
zVal = -DeltaWVal/((1 + aVal) Delta0Val);

etaReg = 10^-30;

I51 =
  N[
    Normal[
      Series[
        Gamma[eps]^2/Xval *
          (1 + aVal)^(-eps) *
          (Delta0Val - I etaReg)^(-2 eps) *
          Hypergeometric2F1[2 eps, eps, 2 eps, zVal],
        {eps, 0, 0}
      ]
    ],
    30
  ];

I53 =
  N[
    Normal[
      Series[
        Gamma[eps]^2/Xval *
          (Delta0Val - I etaReg)^(-eps) *
          (Delta1Val - I etaReg)^(-eps),
        {eps, 0, 0}
      ]
    ],
    30
  ];

Print["Outer Eq. (29) direct GaugeLink comparison"];
Print["This target runs AMFlow if executed."];
Print["Precision goal: ", precisionGoal];
Print["Epsilon order passed to SolveIntegralsGaugeLink: ", epsOrder];
Print["Compared coefficient powers: ", coeffPowers];
Print["Denominators:"];
Print["  D1 = l^2 - ml2"];
Print["  D2 = (l - p)^2 - Ml2"];
Print["  D3 = C[l] = muKlambda - a (p-l)^2 - b p^2"];
Print["  Lplus = x - n.l"];
Print["  Lminus = -(x - n.l)"];
Print["Kinematic point: x=", xval, ", y=", yval, ", s=", sval];
Print["Mass squares: ml2=", ml2val, ", Ml2=", Ml2val, ", mk2=", mk2val, ", Mk2=", Mk2val];
Print["X=", Xval, ", lambda=", lambdaVal, ", a=", aVal, ", b=", bVal];
Print["muKlambda=", muKlambdaVal];
Print["Delta0=", Delta0Val];
Print["DeltaW=", DeltaWVal];
Print["Delta1=", Delta1Val];
Print["z=", zVal];
Print["Analytic Eq. (51):"];
Print[I51];
Print["Analytic Eq. (53):"];
Print[I53];
Print["Eq. (51) - Eq. (53):"];
Print[N[I51 - I53, 30]];
Print["Expected leading eps^-2 coefficient: ", 1/Xval];
Print["The AMFlow target uses j[..., 1, 1, eps, 1]. If AMFlow rejects this, epsilon-dependent indices are unsupported in this setup."];

DefineOuterEq29PlusFamily[];
plusTarget = j[outereq29plus, 1, 1, eps, 1];
Print["Running +L GaugeLink integral: ", plusTarget];
plusResult = SolveIntegralsGaugeLink[{plusTarget}, precisionGoal, epsOrder];
Fplus = plusTarget /. plusResult;
Print["+L result:"];
Print[plusResult];

DefineOuterEq29MinusFamily[];
minusTarget = j[outereq29minus, 1, 1, eps, 1];
Print["Running -L GaugeLink integral: ", minusTarget];
minusResult = SolveIntegralsGaugeLink[{minusTarget}, precisionGoal, epsOrder];
Jminus = minusTarget /. minusResult;
Print["-L result:"];
Print[minusResult];

(* Same GaugeLink convention as the successful one-loop check:
   Fplus = integral with 1/(L+i0), Jminus = integral with 1/(-L+i0),
   and 1/(L-i0) = -1/(-L+i0). *)
rawDeltaIntegral = (-Jminus - Fplus)/(2 Pi I);
reconstructedEq29 = Gamma[eps]/Xval * rawDeltaIntegral;

Print["Raw delta-reconstructed l integral (-Jminus - Fplus)/(2 Pi I):"];
Print[rawDeltaIntegral];
Print["Reconstructed Eq. (29), Gamma[eps]/X times raw integral:"];
Print[reconstructedEq29];

comparison =
  Table[
    With[
      {
        cAM = N[Coefficient[reconstructedEq29, eps, p], 30],
        cAN = N[Coefficient[I53, eps, p], 30]
      },
      <|
        "power" -> p,
        "amflowEq29" -> cAM,
        "analyticEq53" -> cAN,
        "difference" -> N[cAM - cAN, 30],
        "ratio" -> N[cAM/cAN, 30]
      |>
    ],
    {p, coeffPowers}
  ];

Print["Coefficient comparison against Eq. (53):"];
Print[comparison];

fullResult =
  <|
    "metadata" -> <|
      "description" -> "Direct Eq. (29) GaugeLink comparison against Eq. (51)/(53).",
      "precisionGoal" -> precisionGoal,
      "epsOrder" -> epsOrder,
      "coeffPowers" -> coeffPowers,
      "x" -> xval,
      "y" -> yval,
      "s" -> sval,
      "ml2" -> ml2val,
      "Ml2" -> Ml2val,
      "mk2" -> mk2val,
      "Mk2" -> Mk2val,
      "X" -> Xval,
      "lambda" -> lambdaVal,
      "a" -> aVal,
      "b" -> bVal,
      "muKlambda" -> muKlambdaVal,
      "Delta0" -> Delta0Val,
      "DeltaW" -> DeltaWVal,
      "Delta1" -> Delta1Val,
      "z" -> zVal,
      "GaugeLinkConvention" -> "rawDeltaIntegral = (-Jminus - Fplus)/(2 Pi I); reconstructedEq29 = Gamma[eps]/X rawDeltaIntegral",
      "targetIndexPattern" -> "j[..., 1, 1, eps, 1]"
    |>,
    "plusResult" -> plusResult,
    "minusResult" -> minusResult,
    "rawDeltaIntegral" -> rawDeltaIntegral,
    "reconstructedEq29" -> reconstructedEq29,
    "analyticEq51" -> I51,
    "analyticEq53" -> I53,
    "comparison" -> comparison
  |>;

Export[
  FileNameJoin[{$ResultsDirectory, "compare_outer_eq29_result.wl"}],
  fullResult
];

Export[
  FileNameJoin[{$ResultsDirectory, "compare_outer_eq29_result.txt"}],
  ToString[fullResult, InputForm],
  "Text"
];
