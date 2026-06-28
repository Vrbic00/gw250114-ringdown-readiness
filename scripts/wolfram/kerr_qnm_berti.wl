(* ::Package:: *)

(* Kerr QNM benchmark from Berti, Cardoso & Will, PRD 73, 064030 (2006).
   Usage:
     wolframscript -file scripts/wolfram/kerr_qnm_berti.wl MfSolar spin
   Example:
     wolframscript -file scripts/wolfram/kerr_qnm_berti.wl 70 0.7
*)

ClearAll["Global`*"];

solarMassTimeSeconds = 4.925490947*10^-6;

qnmCoefficients = <|
  "220" -> <|"l" -> 2, "m" -> 2, "n" -> 0,
    "f1" -> 1.5251, "f2" -> -1.1568, "f3" -> 0.1292,
    "q1" -> 0.7000, "q2" -> 1.4187, "q3" -> -0.4990|>,
  "221" -> <|"l" -> 2, "m" -> 2, "n" -> 1,
    "f1" -> 1.3673, "f2" -> -1.0260, "f3" -> 0.1628,
    "q1" -> 0.1000, "q2" -> 0.5436, "q3" -> -0.4731|>,
  "222" -> <|"l" -> 2, "m" -> 2, "n" -> 2,
    "f1" -> 1.3223, "f2" -> -1.0257, "f3" -> 0.1860,
    "q1" -> -0.1000, "q2" -> 0.4206, "q3" -> -0.4256|>,
  "330" -> <|"l" -> 3, "m" -> 3, "n" -> 0,
    "f1" -> 1.8956, "f2" -> -1.3043, "f3" -> 0.1818,
    "q1" -> 0.9000, "q2" -> 2.3430, "q3" -> -0.4810|>,
  "440" -> <|"l" -> 4, "m" -> 4, "n" -> 0,
    "f1" -> 2.3000, "f2" -> -1.5056, "f3" -> 0.2244,
    "q1" -> 1.1929, "q2" -> 3.1191, "q3" -> -0.4825|>
|>;

dimensionlessAngularFrequency[mode_String, spin_?NumericQ] := Module[
  {c = qnmCoefficients[mode]},
  c["f1"] + c["f2"] (1 - spin)^c["f3"]
];

qualityFactor[mode_String, spin_?NumericQ] := Module[
  {c = qnmCoefficients[mode]},
  c["q1"] + c["q2"] (1 - spin)^c["q3"]
];

kerrQNM[mode_String, massSolar_?NumericQ, spin_?NumericQ] := Module[
  {mSeconds, omegaDimless, omegaSI, fHz, q, tauSeconds},
  If[! KeyExistsQ[qnmCoefficients, mode],
    Return[Failure["UnknownMode", <|"Mode" -> mode|>]]
  ];
  If[massSolar <= 0,
    Return[Failure["InvalidMass", <|"MassSolar" -> massSolar|>]]
  ];
  If[spin < 0 || spin >= 1,
    Return[Failure["InvalidSpin", <|"Spin" -> spin|>]]
  ];

  mSeconds = massSolar solarMassTimeSeconds;
  omegaDimless = dimensionlessAngularFrequency[mode, spin];
  omegaSI = omegaDimless/mSeconds;
  fHz = omegaSI/(2 Pi);
  q = qualityFactor[mode, spin];
  tauSeconds = 2 q/omegaSI;

  <|
    "mode" -> mode,
    "l" -> qnmCoefficients[mode]["l"],
    "m" -> qnmCoefficients[mode]["m"],
    "n" -> qnmCoefficients[mode]["n"],
    "Momega" -> omegaDimless,
    "f_Hz" -> fHz,
    "Q" -> q,
    "tau_s" -> tauSeconds,
    "tau_ms" -> 1000 tauSeconds
  |>
];

deformedQNM[
  mode_String,
  massSolar_?NumericQ,
  spin_?NumericQ,
  deltaF_: 0,
  deltaTau_: 0
] := Module[{base = kerrQNM[mode, massSolar, spin]},
  If[FailureQ[base], Return[base]];
  Join[
    base,
    <|
      "delta_f" -> deltaF,
      "delta_tau" -> deltaTau,
      "f_deformed_Hz" -> base["f_Hz"] (1 + deltaF),
      "tau_deformed_s" -> base["tau_s"] (1 + deltaTau),
      "tau_deformed_ms" -> base["tau_ms"] (1 + deltaTau)
    |>
  ]
];

formatNumber[x_?NumericQ, spec_List] := ToString[NumberForm[N[x], spec], OutputForm];

formatRow[row_Association] := StringRiffle[
  {
    row["mode"],
    formatNumber[row["Momega"], {8, 5}],
    formatNumber[row["f_Hz"], {9, 3}],
    formatNumber[row["Q"], {8, 4}],
    formatNumber[row["tau_ms"], {8, 3}]
  },
  "\t"
];

args = Rest[$ScriptCommandLine];

If[Length[args] < 2,
  Print["Usage: wolframscript -file scripts/wolfram/kerr_qnm_berti.wl MfSolar spin"];
  Print["Example: wolframscript -file scripts/wolfram/kerr_qnm_berti.wl 70 0.7"];
  Exit[1];
];

massSolar = ToExpression[args[[1]]];
spin = ToExpression[args[[2]]];
modes = {"220", "221", "222", "330", "440"};
rows = kerrQNM[#, massSolar, spin] & /@ modes;

If[AnyTrue[rows, FailureQ],
  Print[rows];
  Exit[1];
];

Print["Kerr QNM benchmark using Berti-Cardoso-Will fitting formulas"];
Print["Mf [Msun] = ", massSolar, ", spin j = ", spin];
Print["mode\tMomega\tf_Hz\tQ\ttau_ms"];
Scan[Print[formatRow[#]] &, rows];
