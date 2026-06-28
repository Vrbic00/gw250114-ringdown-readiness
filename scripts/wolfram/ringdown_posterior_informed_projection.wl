(* ::Package:: *)

(* Posterior-informed synthetic EFT spectral projection.

   Usage:
     wolframscript -file scripts/wolfram/ringdown_posterior_informed_projection.wl
*)

ClearAll["Global`*"];

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[
    {Directory[], "config", "ringdown_posterior_informed_projection.wl"}
  ]
];
outputDir = If[
  Length[args] >= 2,
  args[[2]],
  FileNameJoin[
    {Directory[], "results", "ringdown_posterior_informed_projection"}
  ]
];

If[! FileExistsQ[configPath],
  Print["Configuration not found: ", configPath];
  Exit[1];
];

Get[configPath];
If[! AssociationQ[posteriorInformedProjectionConfig],
  Print["Configuration must define posteriorInformedProjectionConfig."];
  Exit[1];
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

config = posteriorInformedProjectionConfig;
event = config["event"];
fitPath = config["fit_path"];
kerrPath = config["kerr_numeric_path"];
posteriorSamplesPath = config["posterior_samples_path"];
modeSets = config["mode_sets"];
operatorOrder = config["operators"];
branchOrder = config["polarizations"];
scenarios = config["uncertainty_scenarios"];
priorScenarios = config["remnant_prior_scenarios"];
alphaReference = N[config["alpha_reference"]];
spinStep = N[config["spin_derivative_step"]];
mass0 = N[event["mass_detector_msun"]];
spin0 = N[event["spin"]];

If[
  ! FileExistsQ[fitPath] || ! FileExistsQ[kerrPath] ||
    ! FileExistsQ[posteriorSamplesPath],
  Print["Input file missing."];
  Exit[1];
];

importAssociations[path_] := Module[{raw = Import[path, "CSV"]},
  AssociationThread[First[raw], #] & /@ Rest[raw]
];

fitRows = importAssociations[fitPath];
kerrRows = importAssociations[kerrPath];
posteriorRows = importAssociations[posteriorSamplesPath];

num[value_] := N[If[NumericQ[value], value, ToExpression[ToString[value]]]];
txt[value_] := ToString[value];

qnmCoefficients = <|
  "220" -> <|"f1" -> 1.5251, "f2" -> -1.1568, "f3" -> 0.1292,
    "q1" -> 0.7000, "q2" -> 1.4187, "q3" -> -0.4990|>,
  "221" -> <|"f1" -> 1.3673, "f2" -> -1.0260, "f3" -> 0.1628,
    "q1" -> 0.1000, "q2" -> 0.5436, "q3" -> -0.4731|>,
  "222" -> <|"f1" -> 1.3223, "f2" -> -1.0257, "f3" -> 0.1860,
    "q1" -> -0.1000, "q2" -> 0.4206, "q3" -> -0.4256|>,
  "330" -> <|"f1" -> 1.8956, "f2" -> -1.3043, "f3" -> 0.1818,
    "q1" -> 0.9000, "q2" -> 2.3430, "q3" -> -0.4810|>,
  "440" -> <|"f1" -> 2.3000, "f2" -> -1.5056, "f3" -> 0.2244,
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

bertiMomega[mode_String, spin_?NumericQ] := Module[{omegaR, q},
  omegaR = dimensionlessAngularFrequency[mode, spin];
  q = qualityFactor[mode, spin];
  omegaR - I omegaR/(2 q)
];

compatibleKerrRows = Select[
  kerrRows,
  Abs[num[#["mass_detector_msun"]] - mass0] < 10^-9 &&
    Abs[num[#["spin"]] - spin0] < 10^-9 &
];

kerrExactByMode = Association[
  (txt[#["mode"]] ->
    (num[#["qnm_Momega_real"]] + I num[#["qnm_Momega_imag"]])) & /@
    compatibleKerrRows
];

allConfiguredModes = Union[Flatten[Lookup[modeSets, "modes"]]];
missingKerrModes = Select[allConfiguredModes, ! KeyExistsQ[kerrExactByMode, #] &];
If[Length[missingKerrModes] > 0,
  Print["Numerical Kerr baseline missing modes: ", missingKerrModes];
  Exit[1];
];

kerrAnchorByMode = Association[
  (# -> (kerrExactByMode[#] - bertiMomega[#, spin0])) & /@
    allConfiguredModes
];

kerrMomega[mode_String, spin_?NumericQ] :=
  bertiMomega[mode, spin] + kerrAnchorByMode[mode];

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

observableVector[modes_List, spin_?NumericQ] := Flatten[
  Table[
    Module[{momega = kerrMomega[mode, spin], omegaR, omegaI},
      omegaR = Re[momega];
      omegaI = Im[momega];
      {Log[omegaR], -Log[-omegaI]}
    ],
    {mode, modes}
  ]
];

spinColumn[modes_List] :=
  (observableVector[modes, spin0 + spinStep] -
    observableVector[modes, spin0 - spinStep])/(2 spinStep);

lnMassColumn[modes_List] :=
  Flatten[ConstantArray[{-1., 1.}, Length[modes]]];

alphaColumn[modes_List, operator_String, branch_String] := Flatten[
  Table[
    Module[
      {
        base = kerrMomega[mode, spin0],
        shift = shiftPolynomial[mode, operator, branch, spin0]
      },
      {Re[shift]/Re[base], -Im[shift]/Im[base]}
    ],
    {mode, modes}
  ]
];

sigmaVector[modes_List, scenario_Association] :=
  Flatten[
    ConstantArray[
      {num[scenario["sigma_lnf"]], num[scenario["sigma_lntau"]]},
      Length[modes]
    ]
  ];

weightedProfileInfo[
  alphaVec_List,
  nuisanceColumns_List,
  sigma_List,
  priorPrecision_?MatrixQ
] := Module[
  {w, nMat, raw, gram, rhs, coeffs, profiled},
  w = DiagonalMatrix[1/sigma^2];
  nMat = Transpose[nuisanceColumns];
  raw = alphaVec.w.alphaVec;
  gram = Transpose[nMat].w.nMat + priorPrecision;
  rhs = Transpose[nMat].w.alphaVec;
  coeffs = PseudoInverse[gram].rhs;
  profiled = Max[0, raw - rhs.coeffs];
  <|"raw" -> raw, "profiled" -> profiled, "coefficients" -> coeffs|>
];

safeLimit[info_?NumericQ, deltaChi2_?NumericQ] :=
  If[info <= 10^-8, Infinity, Sqrt[deltaChi2/info]];

posteriorNuisanceValues = ({Log[num[#["final_mass"]]/mass0],
      num[#["final_spin"]] - spin0} &) /@ posteriorRows;
remnantCovariance = Covariance[posteriorNuisanceValues];
remnantSigmas = Sqrt[Diagonal[remnantCovariance]];
zeroPriorPrecision = ConstantArray[0., {2, 2}];

priorPrecision[prior_Association] := Module[
  {scale = prior["scale"]},
  If[scale === Infinity || ToString[scale] == "Infinity",
    zeroPriorPrecision,
    Inverse[(num[scale]^2) remnantCovariance]
  ]
];

comboRows = Flatten[
  Table[
    <|"operator" -> operator, "polarization" -> branch|>,
    {operator, operatorOrder},
    {branch, branchOrder}
  ],
  1
];

projectionRows = Flatten[
  Table[
    Module[
      {
        modes = modeSet["modes"], sigma = sigmaVector[modeSet["modes"],
          scenario], nuisanceColumns, alphaVec, info, alpha1, alpha2,
        alpha3, coeffs, retained
      },
      nuisanceColumns = {lnMassColumn[modes], spinColumn[modes]};
      alphaVec = alphaColumn[modes, combo["operator"], combo["polarization"]];
      info = weightedProfileInfo[
        alphaVec,
        nuisanceColumns,
        sigma,
        priorPrecision[prior]
      ];
      alpha1 = safeLimit[info["profiled"], 1.];
      alpha2 = safeLimit[info["profiled"], 4.];
      alpha3 = safeLimit[info["profiled"], 9.];
      coeffs = info["coefficients"];
      retained = If[info["raw"] <= 0, 0, info["profiled"]/info["raw"]];
      <|
        "mode_set" -> modeSet["name"],
        "modes" -> StringRiffle[modes, "+"],
        "measurement_scenario" -> scenario["name"],
        "remnant_prior" -> prior["name"],
        "prior_scale" -> If[prior["scale"] === Infinity, "inf", prior["scale"]],
        "sigma_lnf" -> num[scenario["sigma_lnf"]],
        "sigma_lntau" -> num[scenario["sigma_lntau"]],
        "operator" -> combo["operator"],
        "polarization" -> combo["polarization"],
        "raw_information_per_alpha2" -> info["raw"],
        "profiled_information_per_alpha2" -> info["profiled"],
        "retained_information_fraction" -> retained,
        "alpha_1sigma_profiled" -> alpha1,
        "alpha_2sigma_profiled" -> alpha2,
        "alpha_3sigma_profiled" -> alpha3,
        "alpha_reference" -> alphaReference,
        "kerr_only_sigma_at_alpha_reference" ->
          alphaReference Sqrt[info["profiled"]],
        "best_fit_delta_lnM_for_alpha_reference" ->
          alphaReference coeffs[[1]],
        "best_fit_delta_chi_for_alpha_reference" ->
          alphaReference coeffs[[2]]
      |>
    ],
    {modeSet, modeSets},
    {scenario, scenarios},
    {prior, priorScenarios},
    {combo, comboRows}
  ],
  3
];

summaryRows = Flatten[
  Table[
    Module[
      {
        rows = Select[
          projectionRows,
          #["mode_set"] == modeSet["name"] &&
            #["measurement_scenario"] == scenario["name"] &&
            #["remnant_prior"] == prior["name"] &
        ],
        finiteAlpha1, best, median, meanRetained, identifiable
      },
      finiteAlpha1 = Select[Lookup[rows, "alpha_1sigma_profiled"], NumericQ];
      best = If[Length[finiteAlpha1] == 0, Infinity, Min[finiteAlpha1]];
      median = If[Length[finiteAlpha1] == 0, Infinity, Median[finiteAlpha1]];
      meanRetained = Mean[Lookup[rows, "retained_information_fraction"]];
      identifiable = Count[
        Lookup[rows, "kerr_only_sigma_at_alpha_reference"],
        x_?NumericQ /; x >= 3
      ];
      <|
        "mode_set" -> modeSet["name"],
        "modes" -> StringRiffle[modeSet["modes"], "+"],
        "measurement_scenario" -> scenario["name"],
        "remnant_prior" -> prior["name"],
        "best_alpha_1sigma_profiled" -> best,
        "median_alpha_1sigma_profiled" -> median,
        "mean_retained_information_fraction" -> meanRetained,
        "operators_above_3sigma_at_alpha_reference" -> identifiable,
        "operator_count" -> Length[rows]
      |>
    ],
    {modeSet, modeSets},
    {scenario, scenarios},
    {prior, priorScenarios}
  ],
  2
];

csvValue[x_] := Which[
  x === Infinity, "inf",
  x === -Infinity, "-inf",
  Head[x] === DirectedInfinity, "inf",
  NumericQ[x], N[x],
  True, x
];

exportAssociationCSV[path_, rows_List] := Module[{fields = Keys[First[rows]]},
  Export[path, Prepend[csvValue /@ Lookup[#, fields] & /@ rows, fields], "CSV"]
];

fmt[Infinity, digits_Integer : 4] := "inf";
fmt[value_?NumericQ, digits_Integer : 4] := Module[
  {x = N[value], ax, exponent, mantissa},
  If[x == 0, Return["0"]];
  ax = Abs[x];
  If[10^-3 <= ax < 10^6,
    Return[
      ToString[
        NumberForm[x, {14, digits}, ExponentFunction -> (Null &)],
        OutputForm
      ]
    ]
  ];
  exponent = Floor[Log10[ax]];
  mantissa = x/10.^exponent;
  ToString[
    NumberForm[mantissa, {10, digits}, ExponentFunction -> (Null &)],
    OutputForm
  ] <> "e" <> ToString[exponent]
];
fmt[value_, digits_Integer : 4] := ToString[value];

projectionCsvPath = FileNameJoin[
  {outputDir, "posterior_informed_spectral_projection.csv"}
];
summaryCsvPath = FileNameJoin[
  {outputDir, "posterior_informed_spectral_summary.csv"}
];
priorCsvPath = FileNameJoin[
  {outputDir, "nrSur7dq4_remnant_prior_covariance.csv"}
];
plotPath = FileNameJoin[
  {outputDir, "posterior_informed_alpha_1sigma.png"}
];
reportPath = FileNameJoin[
  {outputDir, "posterior_informed_projection_report.md"}
];

priorRows = Flatten[
  Table[
    <|
      "row" -> {"delta_lnM", "delta_chi"}[[i]],
      "column" -> {"delta_lnM", "delta_chi"}[[j]],
      "covariance" -> remnantCovariance[[i, j]],
      "correlation" -> remnantCovariance[[i, j]]/
        (remnantSigmas[[i]] remnantSigmas[[j]])
    |>,
    {i, 2},
    {j, 2}
  ],
  1
];

exportAssociationCSV[projectionCsvPath, projectionRows];
exportAssociationCSV[summaryCsvPath, summaryRows];
exportAssociationCSV[priorCsvPath, priorRows];

plotRows = Select[
  summaryRows,
  #["mode_set"] == "220_221" &&
    #["measurement_scenario"] == "moderate_2pct_f_10pct_tau" &
];
plotLabels = Lookup[plotRows, "remnant_prior"];
plotValues = Log10 /@ Lookup[plotRows, "best_alpha_1sigma_profiled"];

plot = BarChart[
  plotValues,
  ChartLabels -> Placed[plotLabels, Below],
  Frame -> True,
  Axes -> False,
  FrameLabel -> {None, "log10 best profiled 1-sigma coupling scale"},
  PlotLabel ->
    "220+221 EFT projection with NRSur7dq4 remnant prior comparison",
  GridLines -> {None, Automatic},
  ImagePadding -> {{70, 40}, {90, 60}},
  ImageSize -> 1000
];
Export[plotPath, plot, ImageResolution -> 144];

summaryTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["mode_set"],
        row["measurement_scenario"],
        row["remnant_prior"],
        fmt[row["best_alpha_1sigma_profiled"], 5],
        fmt[row["median_alpha_1sigma_profiled"], 5],
        fmt[row["mean_retained_information_fraction"], 4],
        ToString[row["operators_above_3sigma_at_alpha_reference"]] <>
          "/" <> ToString[row["operator_count"]]
      },
      " | "
    ] <> " |"
  ],
  Select[
    summaryRows,
    #["measurement_scenario"] == "moderate_2pct_f_10pct_tau" &
  ]
];

bestRows = TakeSmallestBy[
  Select[
    projectionRows,
    #["mode_set"] == "220_221" &&
      #["measurement_scenario"] == "moderate_2pct_f_10pct_tau" &&
      #["remnant_prior"] == "nrSur7dq4_imr_prior" &&
      NumericQ[#["alpha_1sigma_profiled"]] &
  ],
  #["alpha_1sigma_profiled"] &,
  UpTo[5]
];

bestTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["operator"],
        row["polarization"],
        fmt[row["alpha_1sigma_profiled"], 6],
        fmt[row["retained_information_fraction"], 4],
        fmt[row["kerr_only_sigma_at_alpha_reference"], 3]
      },
      " | "
    ] <> " |"
  ],
  bestRows
];

report = StringRiffle[
  Join[
    {
      "# Posterior-Informed Synthetic EFT Projection",
      "",
      "This report is generated by `scripts/wolfram/ringdown_posterior_informed_projection.wl`.",
      "",
      "## Scope",
      "",
      "- This is not an EFT observation or a ringdown-only likelihood.",
      "- The NRSur7dq4 full-IMR posterior is used only as a Gaussian prior on local remnant nuisance parameters `{delta ln M_f, delta chi}`.",
      "- The measurement widths for `{log f, log tau}` remain synthetic scenarios.",
      "- The comparison shows how strongly an external GR-informed remnant prior can break the mass/spin degeneracy.",
      "",
      "## Remnant Prior Widths",
      "",
      "- `sigma(delta ln M_f) = " <> fmt[remnantSigmas[[1]], 6] <> "`",
      "- `sigma(delta chi) = " <> fmt[remnantSigmas[[2]], 6] <> "`",
      "- `corr(delta ln M_f, delta chi) = " <>
        fmt[priorRows[[2]]["correlation"], 4] <> "`",
      "",
      "## Moderate Measurement Scenario Summary",
      "",
      "| mode set | measurement scenario | remnant prior | best alpha 1sigma | median alpha 1sigma | mean retained information | above 3sigma at alpha_ref |",
      "| --- | --- | --- | ---: | ---: | ---: | ---: |"
    },
    summaryTableRows,
    {
      "",
      "## Most Detectable 220+221 Fingerprints With NRSur7dq4 Prior",
      "",
      "| operator | polarization | alpha 1sigma | retained information | Kerr-only sigma at alpha_ref |",
      "| --- | --- | ---: | ---: | ---: |"
    },
    bestTableRows,
    {
      "",
      "## Interpretation",
      "",
      "- Free remnant profiling is the conservative synthetic ringdown limit.",
      "- Adding the NRSur7dq4 IMR prior makes even a 220-only spectrum partly informative, but this imports GR waveform information and should not be advertised as an independent no-hair test.",
      "- The 220+221 case remains the physically cleaner path because the second mode adds an internal spectral consistency lever arm.",
      "",
      "## Next Defensible Step",
      "",
      "Look for a public ringdown-specific posterior or mode-frequency posterior for GW250114. If none is available below the data threshold, the next publishable route is to present the current result as a reproducible theory-to-synthetic-likelihood pipeline and clearly mark the observational calibration as IMR-prior-informed.",
      "",
      "## Generated Files",
      "",
      "- `posterior_informed_spectral_projection.csv`",
      "- `posterior_informed_spectral_summary.csv`",
      "- `nrSur7dq4_remnant_prior_covariance.csv`",
      "- `posterior_informed_alpha_1sigma.png`"
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated posterior-informed EFT projection"];
Print["Report: ", reportPath];
Print["sigma_delta_lnM: ", fmt[remnantSigmas[[1]], 6]];
Print["sigma_delta_chi: ", fmt[remnantSigmas[[2]], 6]];
Print["Moderate 220+221 best alpha by prior:"];
Scan[
  Print[
    StringRiffle[
      {
        #["remnant_prior"],
        fmt[#["best_alpha_1sigma_profiled"], 6],
        fmt[#["median_alpha_1sigma_profiled"], 6]
      },
      "\t"
    ]
  ] &,
  Select[
    summaryRows,
    #["mode_set"] == "220_221" &&
      #["measurement_scenario"] == "moderate_2pct_f_10pct_tau" &
  ]
];
