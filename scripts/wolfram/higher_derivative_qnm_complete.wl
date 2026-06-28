(* ::Package:: *)

(* Complete selected-mode higher-derivative QNM spectrum.

   Source:
     Cano et al., arXiv:2409.04517
     https://github.com/pacmn91/BeyondKerrQNM

   Usage:
     wolframscript -file scripts/wolfram/higher_derivative_qnm_complete.wl
*)

ClearAll["Global`*"];

solarMassTimeSeconds = 4.925490947*10^-6;

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[{Directory[], "config", "higher_derivative_qnm_complete.wl"}]
];
outputDir = If[
  Length[args] >= 2,
  args[[2]],
  FileNameJoin[{Directory[], "results", "higher_derivative_qnm_complete"}]
];

If[! FileExistsQ[configPath],
  Print["Configuration not found: ", configPath];
  Exit[1];
];

Get[configPath];
If[! AssociationQ[completeEFTQNMConfig],
  Print["Configuration must define association completeEFTQNMConfig."];
  Exit[1];
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

event = completeEFTQNMConfig["event"];
fitPath = completeEFTQNMConfig["fit_path"];
kerrPath = completeEFTQNMConfig["kerr_numeric_path"];
modes = completeEFTQNMConfig["modes"];
exampleAlpha = N[completeEFTQNMConfig["example_alpha"]];

If[! FileExistsQ[fitPath] || ! FileExistsQ[kerrPath],
  Print["Fit or numerical Kerr input is missing."];
  Exit[1];
];

importAssociations[path_] := Module[{raw = Import[path, "CSV"]},
  AssociationThread[First[raw], #] & /@ Rest[raw]
];

fitRows = importAssociations[fitPath];
kerrRows = importAssociations[kerrPath];

num[value_] := N[If[NumericQ[value], value, ToExpression[ToString[value]]]];
txt[value_] := ToString[value];

massSolar = num[event["mass_detector_msun"]];
spin = num[event["spin"]];

compatibleKerrRows = Select[
  kerrRows,
  Abs[num[#["mass_detector_msun"]] - massSolar] < 10^-9 &&
    Abs[num[#["spin"]] - spin] < 10^-9 &
];

kerrByMode = Association[
  (txt[#["mode"]] ->
    (num[#["qnm_Momega_real"]] + I num[#["qnm_Momega_imag"]])) & /@
    compatibleKerrRows
];

missingModes = Select[modes, ! KeyExistsQ[kerrByMode, #] &];
If[Length[missingModes] > 0,
  Print["Numerical Kerr baseline missing modes: ", missingModes];
  Exit[1];
];

operatorOrder = {
  "lambda_ev", "lambda_odd", "epsilon1", "epsilon2", "epsilon3"
};
branchOrder = {"plus", "minus"};

shiftPolynomial[
  mode_String,
  operator_String,
  branch_String,
  chi_?NumericQ
] := Module[{rows},
  rows = Select[
    fitRows,
    txt[#["mode"]] == mode &&
      txt[#["operator"]] == operator &&
      txt[#["branch"]] == branch &
  ];
  If[Length[rows] == 0,
    Return[
      Failure[
        "MissingFit",
        <|"mode" -> mode, "operator" -> operator, "branch" -> branch|>
      ]
    ]
  ];
  Total[
    (num[#["coefficient_re"]] + I num[#["coefficient_im"]]) *
      chi^Round[num[#["k"]]] & /@ rows
  ]
];

physicalObservables[momega_?NumericQ] := Module[
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

sourceCommit = First[fitRows]["source_commit"];
sourceRepository = First[fitRows]["source_repository"];
sourcePaper = First[fitRows]["source_paper"];
sourceLicense = First[fitRows]["license"];

baselineRows = Map[
  Function[mode,
    Join[
      <|"mode" -> mode, "source" -> "Python qnm numerical Kerr sequence"|>,
      physicalObservables[kerrByMode[mode]]
    ]
  ],
  modes
];

sensitivityRows = Flatten[
  Table[
    Module[
      {
        base = kerrByMode[mode],
        shift = shiftPolynomial[mode, operator, branch, spin],
        baseObs, correctedObs, corrected, fitSubset
      },
      baseObs = physicalObservables[base];
      corrected = base + exampleAlpha shift;
      correctedObs = physicalObservables[corrected];
      fitSubset = Select[
        fitRows,
        txt[#["mode"]] == mode &&
          txt[#["operator"]] == operator &&
          txt[#["branch"]] == branch &
      ];
      <|
        "mode" -> mode,
        "operator" -> operator,
        "parity" -> txt[First[fitSubset]["parity"]],
        "polarization" -> branch,
        "maximum_fit_order" -> Max[num /@ Lookup[fitSubset, "maximum_fit_order"]],
        "example_alpha" -> exampleAlpha,
        "kerr_Momega_re" -> Re[base],
        "kerr_Momega_im" -> Im[base],
        "shift_per_alpha_re" -> Re[shift],
        "shift_per_alpha_im" -> Im[shift],
        "relative_complex_sensitivity_per_alpha" -> Abs[shift/base],
        "dln_frequency_dalpha" -> Re[shift]/Re[base],
        "dln_tau_dalpha" -> -Im[shift]/Im[base],
        "example_frequency_shift_percent" ->
          100 (correctedObs["f_Hz"]/baseObs["f_Hz"] - 1),
        "example_tau_shift_percent" ->
          100 (correctedObs["tau_ms"]/baseObs["tau_ms"] - 1),
        "corrected_f_Hz" -> correctedObs["f_Hz"],
        "corrected_tau_ms" -> correctedObs["tau_ms"]
      |>
    ],
    {mode, modes},
    {operator, operatorOrder},
    {branch, branchOrder}
  ],
  2
];

modeSummaryRows = Map[
  Function[mode,
    Module[{rows, largest},
      rows = Select[sensitivityRows, #["mode"] == mode &];
      largest = First[
        ReverseSortBy[rows, #["relative_complex_sensitivity_per_alpha"] &]
      ];
      <|
        "mode" -> mode,
        "largest_operator" -> largest["operator"],
        "largest_polarization" -> largest["polarization"],
        "maximum_relative_complex_sensitivity_per_alpha" ->
          largest["relative_complex_sensitivity_per_alpha"],
        "example_frequency_shift_percent_at_max" ->
          largest["example_frequency_shift_percent"],
        "example_tau_shift_percent_at_max" ->
          largest["example_tau_shift_percent"]
      |>
    ]
  ],
  modes
];

referenceSensitivity[mode_, operator_, branch_] := First[
  Select[
    sensitivityRows,
    #["mode"] == mode &&
      #["operator"] == operator &&
      #["polarization"] == branch &
  ]
]["relative_complex_sensitivity_per_alpha"];

overtoneRows = Flatten[
  Table[
    Module[
      {
        s220 = referenceSensitivity["220", operator, branch],
        s221 = referenceSensitivity["221", operator, branch],
        s222 = referenceSensitivity["222", operator, branch]
      },
      <|
        "operator" -> operator,
        "polarization" -> branch,
        "relative_sensitivity_220" -> s220,
        "relative_sensitivity_221" -> s221,
        "relative_sensitivity_222" -> s222,
        "amplification_221_over_220" -> s221/s220,
        "amplification_222_over_220" -> s222/s220
      |>
    ],
    {operator, operatorOrder},
    {branch, branchOrder}
  ],
  1
];

exportAssociationCSV[path_, rows_List] := Module[{fields = Keys[First[rows]]},
  Export[path, Prepend[Lookup[#, fields] & /@ rows, fields], "CSV"]
];

baselineCsvPath = FileNameJoin[{outputDir, "gw250114_complete_kerr_baseline.csv"}];
sensitivityCsvPath = FileNameJoin[{outputDir, "gw250114_complete_eft_sensitivities.csv"}];
overtoneCsvPath = FileNameJoin[{outputDir, "gw250114_overtone_amplification.csv"}];
summaryCsvPath = FileNameJoin[{outputDir, "gw250114_mode_sensitivity_summary.csv"}];
plotPath = FileNameJoin[{outputDir, "gw250114_complete_eft_sensitivity_heatmap.png"}];
reportPath = FileNameJoin[{outputDir, "higher_derivative_qnm_complete_report.md"}];

exportAssociationCSV[baselineCsvPath, baselineRows];
exportAssociationCSV[sensitivityCsvPath, sensitivityRows];
exportAssociationCSV[overtoneCsvPath, overtoneRows];
exportAssociationCSV[summaryCsvPath, modeSummaryRows];

columnKeys = Flatten[
  Table[{operator, branch}, {operator, operatorOrder}, {branch, branchOrder}],
  1
];
columnLabels = (#[[1]] <> " " <> #[[2]]) & /@ columnKeys;

heatmapData = Table[
  Log10[
    First[
      Select[
        sensitivityRows,
        #["mode"] == mode &&
          #["operator"] == key[[1]] &&
          #["polarization"] == key[[2]] &
      ]
    ]["relative_complex_sensitivity_per_alpha"]
  ],
  {mode, Reverse[modes]},
  {key, columnKeys}
];

heatmap = ArrayPlot[
  heatmapData,
  Frame -> True,
  FrameTicks -> {
    Thread[{Range[Length[modes]], Reverse[modes]}],
    Thread[{Range[Length[columnLabels]], Rotate[#, Pi/3] & /@ columnLabels}]
  },
  ColorFunction -> "SolarColors",
  PlotLegends -> BarLegend[
    Automatic,
    LegendLabel -> "log10 |delta omega / omega| per unit alpha"
  ],
  PlotLabel -> "GW250114 higher-derivative QNM sensitivity at spin 0.68",
  ImagePadding -> {{70, 190}, {180, 60}},
  ImageSize -> 1200
];
Export[plotPath, heatmap, ImageResolution -> 144];

fmt[value_?NumericQ, digits_Integer: 4] :=
  ToString[NumberForm[N[value], {14, digits}], OutputForm];

baselineTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["mode"],
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

modeSummaryTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["mode"],
        row["largest_operator"],
        row["largest_polarization"],
        fmt[row["maximum_relative_complex_sensitivity_per_alpha"], 3],
        fmt[row["example_frequency_shift_percent_at_max"], 3],
        fmt[row["example_tau_shift_percent_at_max"], 3]
      },
      " | "
    ] <> " |"
  ],
  modeSummaryRows
];

overtoneTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["operator"],
        row["polarization"],
        fmt[row["amplification_221_over_220"], 3],
        fmt[row["amplification_222_over_220"], 3]
      },
      " | "
    ] <> " |"
  ],
  overtoneRows
];

report = StringRiffle[
  Join[
    {
      "# Complete Higher-Derivative QNM Spectrum for GW250114",
      "",
      "This report evaluates the public polynomial fits associated with Cano et al. (2024), arXiv:2409.04517.",
      "",
      "## Data Provenance",
      "",
      "- Repository: " <> sourceRepository,
      "- Imported commit: `" <> sourceCommit <> "`",
      "- License: `" <> sourceLicense <> "`",
      "- Selected modes: `" <> StringRiffle[modes, ", "] <> "`",
      "- The repository commit is newer than the paper release; the exact commit is recorded to keep the calculation reproducible.",
      "",
      "## Physical Model",
      "",
      "- The source supplies a higher-derivative EFT action, perturbative rotating backgrounds, modified gravitational Teukolsky equations, and numerical QNM fits.",
      "- The fits contain both polarizations for all five EFT operators and include the 221 and 222 overtones.",
      "- Frequencies obey `M omega = M omega_Kerr + alpha_q delta omega_q + O(alpha_q^2)`.",
      "- The Kerr baseline below uses the numerical Python `qnm` sequences, not the approximate Berti fit.",
      "",
      "## GW250114 Kerr Baseline",
      "",
      "| mode | Re(Momega) | Im(Momega) | f [Hz] | tau [ms] | Q |",
      "| --- | ---: | ---: | ---: | ---: | ---: |"
    },
    baselineTableRows,
    {
      "",
      "## Mode Sensitivity Summary",
      "",
      "The final two columns use the illustrative coupling `alpha_q = " <>
        ToString[exampleAlpha] <> "` for the most sensitive operator/polarization of each mode.",
      "",
      "| mode | largest operator | polarization | max |delta omega/omega| per alpha | Delta f/f [%] | Delta tau/tau [%] |",
      "| --- | --- | --- | ---: | ---: | ---: |"
    },
    modeSummaryTableRows,
    {
      "",
      "## Overtone Amplification",
      "",
      "Ratios below compare the magnitude `|delta omega/omega|` with the 220 fundamental at the same spin.",
      "",
      "| operator | polarization | 221 / 220 | 222 / 220 |",
      "| --- | --- | ---: | ---: |"
    },
    overtoneTableRows,
    {
      "",
      "## Scientific Reading",
      "",
      "- This removes the main spectral mismatch of the earlier bridge: the theory now supplies both 220 and 221, the two modes central to the current GW250114 spectroscopy result.",
      "- Overtones are often substantially more sensitive than the fundamental, but this also means their first-order EFT regime can be narrower.",
      "- At `chi = 0.68`, the relevant corotating modes lie in the paper's intended accuracy range near `chi approximately 0.7`.",
      "- A QNM spectrum is not yet a ringdown likelihood. Mode amplitudes, excitation, start time, detector noise, remnant priors, and polarization content still have to be modeled.",
      "- Couplings must be tested one at a time first. Mixed parity-preserving and parity-breaking terms require the full polarization eigenvalue rule rather than naive addition.",
      "",
      "## Next Defensible Step",
      "",
      "Build a GW250114 spectral likelihood for the 220+221 content with Kerr remnant mass and spin as nuisance parameters, then map each one-at-a-time EFT coupling into the complex mode frequencies. Start with a synthetic injection/recovery and only then connect to strain or published posterior products.",
      "",
      "## Generated Files",
      "",
      "- `gw250114_complete_kerr_baseline.csv`",
      "- `gw250114_complete_eft_sensitivities.csv`",
      "- `gw250114_overtone_amplification.csv`",
      "- `gw250114_mode_sensitivity_summary.csv`",
      "- `gw250114_complete_eft_sensitivity_heatmap.png`"
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated complete higher-derivative QNM spectrum"];
Print["Source commit: ", sourceCommit];
Print["Report: ", reportPath];
Print["Plot: ", plotPath];
Print["mode\tmax_relative_sensitivity_per_alpha\toperator\tpolarization"];
Scan[
  Print[
    StringRiffle[
      {
        #["mode"],
        fmt[#["maximum_relative_complex_sensitivity_per_alpha"], 4],
        #["largest_operator"],
        #["largest_polarization"]
      },
      "\t"
    ]
  ] &,
  modeSummaryRows
];
