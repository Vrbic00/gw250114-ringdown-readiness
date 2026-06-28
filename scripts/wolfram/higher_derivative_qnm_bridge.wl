(* ::Package:: *)

(* Theory-backed rotating-QNM bridge using the polynomial frequency shifts
   published in Cano et al., arXiv:2307.07431.

   Usage:
     wolframscript -file scripts/wolfram/higher_derivative_qnm_bridge.wl

   Optional:
     wolframscript -file scripts/wolfram/higher_derivative_qnm_bridge.wl config/custom.wl results/custom_dir
*)

ClearAll["Global`*"];

solarMassTimeSeconds = 4.925490947*10^-6;

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[{Directory[], "config", "higher_derivative_qnm_bridge.wl"}]
];
outputDir = If[
  Length[args] >= 2,
  args[[2]],
  FileNameJoin[{Directory[], "results", "higher_derivative_qnm_bridge"}]
];

If[! FileExistsQ[configPath],
  Print["Configuration not found: ", configPath];
  Exit[1];
];

Get[configPath];
If[! AssociationQ[eftQNMConfig],
  Print["Configuration must define association eftQNMConfig."];
  Exit[1];
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

coefficientPath = eftQNMConfig["coefficient_path"];
validationPath = eftQNMConfig["validation_path"];
kerrNumericPath = Lookup[
  eftQNMConfig,
  "kerr_numeric_path",
  Missing["NotConfigured"]
];
event = eftQNMConfig["event"];
exampleAlpha = N[eftQNMConfig["example_alpha"]];
scenarios = eftQNMConfig["scenarios"];

If[! FileExistsQ[coefficientPath] || ! FileExistsQ[validationPath],
  Print["Coefficient or validation data file not found."];
  Exit[1];
];

importAssociations[path_] := Module[{raw = Import[path, "CSV"]},
  AssociationThread[First[raw], #] & /@ Rest[raw]
];

coefficientRows = importAssociations[coefficientPath];
validationReference = importAssociations[validationPath];
kerrNumericRows = If[
  StringQ[kerrNumericPath] && FileExistsQ[kerrNumericPath],
  importAssociations[kerrNumericPath],
  {}
];

qnmFits = <|
  "220" -> <|
    "l" -> 2, "m" -> 2, "n" -> 0,
    "f1" -> 1.5251, "f2" -> -1.1568, "f3" -> 0.1292,
    "q1" -> 0.7000, "q2" -> 1.4187, "q3" -> -0.4990
  |>,
  "330" -> <|
    "l" -> 3, "m" -> 3, "n" -> 0,
    "f1" -> 1.8956, "f2" -> -1.3043, "f3" -> 0.1818,
    "q1" -> 0.9000, "q2" -> 2.3430, "q3" -> -0.4810
  |>
|>;

num[value_] := N[If[NumericQ[value], value, ToExpression[ToString[value]]]];
txt[value_] := ToString[value];

kerrMomega[mode_String, spin_?NumericQ] := Module[{c, omegaR, quality},
  If[! KeyExistsQ[qnmFits, mode],
    Return[Failure["UnknownMode", <|"mode" -> mode|>]]
  ];
  If[spin < 0 || spin >= 1,
    Return[Failure["InvalidSpin", <|"spin" -> spin|>]]
  ];
  c = qnmFits[mode];
  omegaR = c["f1"] + c["f2"] (1 - spin)^c["f3"];
  quality = c["q1"] + c["q2"] (1 - spin)^c["q3"];
  omegaR - I omegaR/(2 quality)
];

massSolar = num[event["mass_detector_msun"]];
spin = num[event["spin"]];

numericBaselineCompatibleQ =
  Length[kerrNumericRows] > 0 &&
  AllTrue[
    kerrNumericRows,
    Abs[num[#["mass_detector_msun"]] - massSolar] < 10^-9 &&
      Abs[num[#["spin"]] - spin] < 10^-9 &
  ];

numericalKerrByMode = If[
  numericBaselineCompatibleQ,
  Association[
    (txt[#["mode"]] ->
      (num[#["qnm_Momega_real"]] + I num[#["qnm_Momega_imag"]])) & /@
      kerrNumericRows
  ],
  <||>
];

eventKerrMomega[mode_String] := Lookup[
  numericalKerrByMode,
  mode,
  kerrMomega[mode, spin]
];

eventKerrSource[mode_String] := If[
  KeyExistsQ[numericalKerrByMode, mode],
  "Python qnm numerical Kerr sequence",
  "Berti-Cardoso-Will Kerr fit"
];

physicalObservables[momega_?NumericQ, massSolar_?NumericQ] := Module[
  {mTime = massSolar solarMassTimeSeconds, omegaR, omegaI},
  omegaR = Re[momega];
  omegaI = Im[momega];
  If[omegaR <= 0 || omegaI >= 0,
    Return[Failure["UnstableOrInvalidMode", <|"Momega" -> momega|>]]
  ];
  <|
    "Momega_re" -> omegaR,
    "Momega_im" -> omegaI,
    "f_Hz" -> omegaR/(2 Pi mTime),
    "tau_ms" -> -1000 mTime/omegaI,
    "Q" -> -omegaR/(2 omegaI)
  |>
];

coefficientVector[row_Association] := Table[
  num[row["c" <> ToString[k] <> "_re"]] +
    I num[row["c" <> ToString[k] <> "_im"]],
  {k, 0, 12}
];

shiftPolynomialFromRow[row_Association, spin_?NumericQ] :=
  Sum[coefficientVector[row][[k + 1]] spin^k, {k, 0, 12}];

findCoefficientRow[
  mode_String,
  operator_String,
  branch_String
] := Module[{matches},
  matches = Select[
    coefficientRows,
    txt[#["mode"]] == mode &&
      txt[#["operator"]] == operator &&
      txt[#["branch"]] == branch &
  ];
  If[
    Length[matches] == 1,
    First[matches],
    Failure[
      "CoefficientLookup",
      <|"mode" -> mode, "operator" -> operator, "branch" -> branch|>
    ]
  ]
];

shiftPolynomial[
  mode_String,
  operator_String,
  branch_String,
  spin_?NumericQ
] := Module[{row = findCoefficientRow[mode, operator, branch]},
  If[FailureQ[row], row, shiftPolynomialFromRow[row, spin]]
];

coupling[couplings_Association, operator_String] :=
  num[Lookup[couplings, operator, 0]];

combinedModeShifts[
  mode_String,
  spin_?NumericQ,
  couplings_Association
] := Module[
  {
    preservingOperators = {"lambda_ev", "epsilon1", "epsilon2"},
    breakingOperators = {"lambda_odd", "epsilon3"},
    plusShift, minusShift, breakingShift, averageShift, split
  },
  plusShift = Total[
    coupling[couplings, #] shiftPolynomial[mode, #, "plus", spin] & /@
      preservingOperators
  ];
  minusShift = Total[
    coupling[couplings, #] shiftPolynomial[mode, #, "minus", spin] & /@
      preservingOperators
  ];
  breakingShift = Total[
    coupling[couplings, #] shiftPolynomial[mode, #, "break", spin] & /@
      breakingOperators
  ];

  Which[
    Chop[breakingShift] == 0,
      {
        <|"polarization" -> "plus", "delta_Momega" -> plusShift|>,
        <|"polarization" -> "minus", "delta_Momega" -> minusShift|>
      },
    Chop[plusShift] == 0 && Chop[minusShift] == 0,
      {
        <|"polarization" -> "plus", "delta_Momega" -> breakingShift|>,
        <|"polarization" -> "minus", "delta_Momega" -> -breakingShift|>
      },
    True,
      averageShift = (plusShift + minusShift)/2;
      split = Sqrt[((plusShift - minusShift)/2)^2 + breakingShift^2];
      {
        <|"polarization" -> "eigen_plus", "delta_Momega" -> averageShift + split|>,
        <|"polarization" -> "eigen_minus", "delta_Momega" -> averageShift - split|>
      }
  ]
];

validationRows = Map[
  Function[reference,
    Module[{row, calculated, paper, difference},
      row = findCoefficientRow[
        txt[reference["mode"]],
        txt[reference["operator"]],
        txt[reference["branch"]]
      ];
      calculated = shiftPolynomialFromRow[row, num[reference["spin"]]];
      paper = num[reference["paper_shift_re"]] +
        I num[reference["paper_shift_im"]];
      difference = calculated - paper;
      <|
        "mode" -> txt[reference["mode"]],
        "operator" -> txt[reference["operator"]],
        "branch" -> txt[reference["branch"]],
        "spin" -> num[reference["spin"]],
        "calculated_re" -> Re[calculated],
        "calculated_im" -> Im[calculated],
        "paper_re" -> Re[paper],
        "paper_im" -> Im[paper],
        "absolute_complex_difference" -> Abs[difference],
        "paper_lower_order_difference_percent" ->
          num[reference["lower_order_difference_percent"]],
        "source_table" -> txt[reference["source_table"]]
      |>
    ]
  ],
  validationReference
];

maxValidationDifference = Max[
  Lookup[validationRows, "absolute_complex_difference"]
];
validationPassed = maxValidationDifference < 0.002;

modes = {"220", "330"};

baselineRows = Map[
  Function[mode,
    Module[{base = eventKerrMomega[mode], fit = kerrMomega[mode, spin]},
      Join[
        <|
          "mode" -> mode,
          "baseline_source" -> eventKerrSource[mode],
          "berti_minus_baseline_re_percent" ->
            100 (Re[fit]/Re[base] - 1),
          "berti_minus_baseline_im_percent" ->
            100 (Im[fit]/Im[base] - 1)
        |>,
        physicalObservables[base, massSolar]
      ]
    ]
  ],
  modes
];

independentShiftRows = Flatten[
  Map[
    Function[row,
      Module[
        {
          mode = txt[row["mode"]],
          operator = txt[row["operator"]],
          parity = txt[row["parity"]],
          storedBranch = txt[row["branch"]],
          base, shift, branches
        },
        base = eventKerrMomega[mode];
        shift = shiftPolynomialFromRow[row, spin];
        branches = If[
          parity == "breaking",
          {{"plus", shift}, {"minus", -shift}},
          {{storedBranch, shift}}
        ];
        Map[
          Function[branchData,
            Module[
              {
                polarization = branchData[[1]],
                unitShift = branchData[[2]],
                corrected, baseObs, correctedObs, frequencyDerivative,
                dampingDerivative
              },
              corrected = base + exampleAlpha unitShift;
              baseObs = physicalObservables[base, massSolar];
              correctedObs = physicalObservables[corrected, massSolar];
              frequencyDerivative = Re[unitShift]/Re[base];
              dampingDerivative = -Im[unitShift]/Im[base];
              <|
                "mode" -> mode,
                "operator" -> operator,
                "parity" -> parity,
                "polarization" -> polarization,
                "example_alpha" -> exampleAlpha,
                "kerr_Momega_re" -> Re[base],
                "kerr_Momega_im" -> Im[base],
                "shift_per_alpha_re" -> Re[unitShift],
                "shift_per_alpha_im" -> Im[unitShift],
                "dln_frequency_dalpha" -> frequencyDerivative,
                "dln_tau_dalpha" -> dampingDerivative,
                "example_frequency_shift_percent" ->
                  100 (correctedObs["f_Hz"]/baseObs["f_Hz"] - 1),
                "example_tau_shift_percent" ->
                  100 (correctedObs["tau_ms"]/baseObs["tau_ms"] - 1),
                "kerr_f_Hz" -> baseObs["f_Hz"],
                "corrected_f_Hz" -> correctedObs["f_Hz"],
                "kerr_tau_ms" -> baseObs["tau_ms"],
                "corrected_tau_ms" -> correctedObs["tau_ms"]
              |>
            ]
          ],
          branches
        ]
      ]
    ],
    coefficientRows
  ],
  2
];

scenarioRows = Flatten[
  Table[
    Module[{scenario = scenarioItem, mode = modeItem, base, shifts, baseObs},
      base = eventKerrMomega[mode];
      baseObs = physicalObservables[base, massSolar];
      shifts = combinedModeShifts[mode, spin, scenario["couplings"]];
      Map[
        Function[shiftRow,
          Module[{corrected, correctedObs},
            corrected = base + shiftRow["delta_Momega"];
            correctedObs = physicalObservables[corrected, massSolar];
            <|
              "scenario" -> scenario["name"],
              "mode" -> mode,
              "polarization" -> shiftRow["polarization"],
              "couplings" -> ToString[scenario["couplings"], InputForm],
              "delta_Momega_re" -> Re[shiftRow["delta_Momega"]],
              "delta_Momega_im" -> Im[shiftRow["delta_Momega"]],
              "frequency_shift_percent" ->
                100 (correctedObs["f_Hz"]/baseObs["f_Hz"] - 1),
              "tau_shift_percent" ->
                100 (correctedObs["tau_ms"]/baseObs["tau_ms"] - 1),
              "f_Hz" -> correctedObs["f_Hz"],
              "tau_ms" -> correctedObs["tau_ms"],
              "Q" -> correctedObs["Q"]
            |>
          ]
        ],
        shifts
      ]
    ],
    {scenarioItem, scenarios},
    {modeItem, modes}
  ],
  3
];

exportAssociationCSV[path_, rows_List] := Module[{fields},
  fields = Keys[First[rows]];
  Export[path, Prepend[Lookup[#, fields] & /@ rows, fields], "CSV"]
];

validationCsvPath = FileNameJoin[{outputDir, "cano_polynomial_validation.csv"}];
baselineCsvPath = FileNameJoin[{outputDir, "gw250114_kerr_220_330_baseline.csv"}];
sensitivityCsvPath = FileNameJoin[{outputDir, "gw250114_eft_qnm_sensitivities.csv"}];
scenarioCsvPath = FileNameJoin[{outputDir, "gw250114_eft_qnm_scenarios.csv"}];
plotPath = FileNameJoin[{outputDir, "gw250114_eft_qnm_fractional_shifts.png"}];
reportPath = FileNameJoin[{outputDir, "higher_derivative_qnm_bridge_report.md"}];

exportAssociationCSV[validationCsvPath, validationRows];
exportAssociationCSV[baselineCsvPath, baselineRows];
exportAssociationCSV[sensitivityCsvPath, independentShiftRows];
exportAssociationCSV[scenarioCsvPath, scenarioRows];

shortLabel[row_Association] := StringJoin[
  row["operator"],
  " ",
  row["polarization"]
];

modePlot[mode_String, field_String, label_String] := Module[
  {rows, values, labels},
  rows = Select[independentShiftRows, #["mode"] == mode &];
  values = Lookup[rows, field];
  labels = shortLabel /@ rows;
  BarChart[
    values,
    BarOrigin -> Left,
    ChartLabels -> Placed[labels, Before],
    Frame -> True,
    FrameLabel -> {None, None},
    PlotLabel -> mode <> ": " <> label <> ", spin = " <> ToString[spin],
    GridLines -> {None, {0}},
    ImagePadding -> {{180, 25}, {55, 35}},
    ImageSize -> 700,
    PlotRange -> All
  ]
];

plot = GraphicsGrid[
  {
    {
      modePlot[
        "220",
        "example_frequency_shift_percent",
        "Delta f / f [%], alpha = " <> ToString[exampleAlpha]
      ],
      modePlot[
        "220",
        "example_tau_shift_percent",
        "Delta tau / tau [%], alpha = " <> ToString[exampleAlpha]
      ]
    },
    {
      modePlot[
        "330",
        "example_frequency_shift_percent",
        "Delta f / f [%], alpha = " <> ToString[exampleAlpha]
      ],
      modePlot[
        "330",
        "example_tau_shift_percent",
        "Delta tau / tau [%], alpha = " <> ToString[exampleAlpha]
      ]
    }
  },
  Spacings -> {0.2, 0.25},
  ImageSize -> 1500
];
Export[plotPath, plot, ImageResolution -> 144];

fmt[value_?NumericQ, digits_Integer: 4] :=
  ToString[NumberForm[N[value], {14, digits}], OutputForm];

baselineTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["mode"],
        row["baseline_source"],
        fmt[row["Momega_re"], 6],
        fmt[row["Momega_im"], 6],
        fmt[row["f_Hz"], 3],
        fmt[row["tau_ms"], 3],
        fmt[row["Q"], 4]
      },
      " | "
    ] <> " |"
  ],
  baselineRows
];

sensitivityTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["mode"],
        row["operator"],
        row["polarization"],
        fmt[row["shift_per_alpha_re"], 4],
        fmt[row["shift_per_alpha_im"], 4],
        fmt[row["example_frequency_shift_percent"], 3],
        fmt[row["example_tau_shift_percent"], 3]
      },
      " | "
    ] <> " |"
  ],
  independentShiftRows
];

report = StringRiffle[
  Join[
    {
      "# Higher-Derivative Rotating-QNM Bridge",
      "",
      "This report implements the gravitational QNM shifts of Cano et al. (2023), arXiv:2307.07431.",
      "",
      "## Scope",
      "",
      "- Full EFT action and corrected gravitational Teukolsky equations are supplied by the source theory.",
      "- Published polynomial shifts are available for the fundamental 220 and 330 modes.",
      "- The calculation is first order in the dimensionless EFT couplings `alpha_q`.",
      "- This is a theory-backed spectral bridge, not yet a fit to GW strain or LVK posterior samples.",
      "",
      "## Polynomial Reproduction",
      "",
      "- Validation spin: `chi = 0.7`.",
      "- Maximum complex difference from the rounded values in Tables V-VI: `" <>
        fmt[maxValidationDifference, 6] <> "`.",
      "- Validation status: `" <> If[validationPassed, "PASS", "FAIL"] <> "`.",
      "- The source estimates roughly 0.9%-4.7% internal differences at `chi = 0.7` between successive spin-expansion orders, depending on mode and operator.",
      "",
      "## GW250114 Kerr Baseline",
      "",
      "- Detector-frame remnant mass: `" <> ToString[massSolar] <> " Msun`.",
      "- Dimensionless remnant spin: `" <> ToString[spin] <> "`.",
      "- Event source: " <> event["source"],
      "",
      "| mode | Kerr baseline | Re(Momega) | Im(Momega) | f [Hz] | tau [ms] | Q |",
      "| --- | --- | ---: | ---: | ---: | ---: | ---: |"
    },
    baselineTableRows,
    {
      "",
      "## Independent Operator Sensitivities",
      "",
      "The last two columns use the illustrative value `alpha_q = " <>
        ToString[exampleAlpha] <> "`. They are not measured constraints.",
      "",
      "| mode | operator | polarization | Re(delta omega)/alpha | Im(delta omega)/alpha | Delta f/f [%] | Delta tau/tau [%] |",
      "| --- | --- | --- | ---: | ---: | ---: | ---: |"
    },
    sensitivityTableRows,
    {
      "",
      "## Interpretation Guardrails",
      "",
      "- `lambda_ev`, `epsilon1`, and `epsilon2` preserve parity and split polar/axial branches.",
      "- `lambda_odd` and `epsilon3` break parity; their two eigenfrequencies receive opposite shifts.",
      "- When preserving and parity-breaking terms coexist, shifts do not combine by naive addition. The scenario evaluator uses the paper's polarization combination rule.",
      "- The EFT requires `|alpha_q| << 1`; large shifts must not be extrapolated outside this regime.",
      "- The published fit is designed for spins through about `chi = 0.7`; use beyond that range is extrapolation.",
      "- Only 220 and 330 fundamental modes are supplied. The 221 overtone central to the current GW250114 spectroscopy result is not yet available in this framework.",
      "- Mode amplitudes, excitation, merger dynamics, and a time-domain waveform are outside this spectral calculation.",
      "- The event table uses the numerical Python `qnm` Kerr sequence when its cached mass and spin match the configuration. Berti fits remain the WL fallback and differ here by roughly 0.5%-0.8%.",
      "",
      "## Generated Files",
      "",
      "- `cano_polynomial_validation.csv`",
      "- `gw250114_kerr_220_330_baseline.csv`",
      "- `gw250114_eft_qnm_sensitivities.csv`",
      "- `gw250114_eft_qnm_scenarios.csv`",
      "- `gw250114_eft_qnm_fractional_shifts.png`",
      "",
      "## Scientific Status",
      "",
      "- Reproduction: polynomial QNM shifts from Cano et al.",
      "- Numerical check: independent evaluation at `chi = 0.7` against Tables V-VI.",
      "- Interpretation: derivatives and illustrative shifts around the GW250114 Kerr remnant.",
      "- Not yet done: event-level likelihood or exclusion on EFT couplings."
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated higher-derivative rotating-QNM bridge"];
Print["Validation: ", If[validationPassed, "PASS", "FAIL"],
  ", max complex difference = ", N[maxValidationDifference]];
Print["Report: ", reportPath];
Print["Plot: ", plotPath];
Print["mode\tf_Kerr_Hz\ttau_Kerr_ms"];
Scan[
  Print[
    StringRiffle[
      {
        #["mode"],
        fmt[#["f_Hz"], 3],
        fmt[#["tau_ms"], 3]
      },
      "\t"
    ]
  ] &,
  baselineRows
];

If[! validationPassed, Exit[2]];
