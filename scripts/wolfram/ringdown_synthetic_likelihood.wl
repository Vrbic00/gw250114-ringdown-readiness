(* ::Package:: *)

(* Synthetic ringdown identifiability for the higher-derivative QNM bridge.

   Usage:
     wolframscript -file scripts/wolfram/ringdown_synthetic_likelihood.wl
     wolframscript -file scripts/wolfram/ringdown_synthetic_likelihood.wl config/custom.wl results/custom_dir
*)

ClearAll["Global`*"];

solarMassTimeSeconds = 4.925490947*10^-6;

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[{Directory[], "config", "ringdown_synthetic_likelihood.wl"}]
];
outputDir = If[
  Length[args] >= 2,
  args[[2]],
  FileNameJoin[{Directory[], "results", "ringdown_synthetic_likelihood"}]
];

If[! FileExistsQ[configPath],
  Print["Configuration not found: ", configPath];
  Exit[1];
];

Get[configPath];
If[! AssociationQ[syntheticRingdownLikelihoodConfig],
  Print["Configuration must define association syntheticRingdownLikelihoodConfig."];
  Exit[1];
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

config = syntheticRingdownLikelihoodConfig;
event = config["event"];
fitPath = config["fit_path"];
kerrPath = config["kerr_numeric_path"];
modeSets = config["mode_sets"];
operatorOrder = config["operators"];
branchOrder = config["polarizations"];
scenarios = config["uncertainty_scenarios"];
alphaReference = N[config["alpha_reference"]];
spinStep = N[config["spin_derivative_step"]];
alphaStep = N[config["alpha_derivative_step"]];

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

mass0 = num[event["mass_detector_msun"]];
spin0 = num[event["spin"]];

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
  (# -> (kerrExactByMode[#] - bertiMomega[#, spin0])) & /@ allConfiguredModes
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

modeObservables[mode_, massSolar_, spin_] :=
  modeObservables[mode, massSolar, spin, "kerr", "none", 0.];

modeObservables[mode_, massSolar_, spin_, operator_, branch_, alpha_] :=
 Module[
  {
    modeString = txt[mode], operatorString = txt[operator],
    branchString = txt[branch], massValue = num[massSolar],
    spinValue = num[spin], alphaValue = num[alpha],
    momega, shift, omegaR, omegaI, fHz, tauMs
  },
  momega = kerrMomega[modeString, spinValue];
  If[operatorString =!= "kerr" && alphaValue != 0,
    shift = shiftPolynomial[modeString, operatorString, branchString, spinValue];
    If[FailureQ[shift], Return[shift]];
    momega = momega + alphaValue shift;
  ];
  omegaR = Re[momega];
  omegaI = Im[momega];
  If[omegaR <= 0 || omegaI >= 0 || massValue <= 0,
    Return[
      Failure[
        "InvalidMode",
        <|"mode" -> modeString, "Momega" -> momega|>
      ]
    ]
  ];
  fHz = omegaR/(2 Pi massValue solarMassTimeSeconds);
  tauMs = -1000 massValue solarMassTimeSeconds/omegaI;
  {Log[fHz], Log[tauMs]}
];

observableVector[modes_List, massSolar_, spin_] :=
  observableVector[modes, massSolar, spin, "kerr", "none", 0.];

observableVector[modes_List, massSolar_, spin_, operator_, branch_, alpha_] :=
  Flatten[
    modeObservables[#, massSolar, spin, operator, branch, alpha] & /@
      modes
  ];

observableLabels[modes_List] :=
  Flatten[({"ln_f_" <> #, "ln_tau_" <> #} & /@ modes)];

sigmaVector[modes_List, scenario_Association] :=
  Flatten[
    ConstantArray[
      {num[scenario["sigma_lnf"]], num[scenario["sigma_lntau"]]},
      Length[modes]
    ]
  ];

lnMassColumn[modes_List] :=
  Flatten[ConstantArray[{-1., 1.}, Length[modes]]];

spinColumn[modes_List] :=
  (
    observableVector[modes, mass0, spin0 + spinStep] -
      observableVector[modes, mass0, spin0 - spinStep]
  )/(2 spinStep);

alphaColumn[modes_List, operator_String, branch_String] :=
  (
    observableVector[modes, mass0, spin0, operator, branch, alphaStep] -
      observableVector[modes, mass0, spin0, operator, branch, -alphaStep]
  )/(2 alphaStep);

weightedFit[target_List, columns_List, sigma_List] := Module[
  {weightedTarget, weightedColumns, coeffs, residual},
  weightedTarget = target/sigma;
  If[Length[columns] == 0,
    Return[
      <|
        "coefficients" -> {},
        "residual" -> weightedTarget,
        "chi2" -> weightedTarget.weightedTarget
      |>
    ]
  ];
  weightedColumns = Transpose[(#/sigma) & /@ columns];
  coeffs = LeastSquares[weightedColumns, weightedTarget];
  residual = weightedTarget - weightedColumns.coeffs;
  <|
    "coefficients" -> coeffs,
    "residual" -> residual,
    "chi2" -> residual.residual
  |>
];

safeLimit[info_?NumericQ, deltaChi2_?NumericQ] :=
  If[info <= 10^-12, Infinity, Sqrt[deltaChi2/info]];

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

comboRows = Flatten[
  Table[
    <|"operator" -> operator, "polarization" -> branch|>,
    {operator, operatorOrder},
    {branch, branchOrder}
  ],
  1
];

detectabilityRows = Flatten[
  Table[
    Module[
      {
        modeSetName = modeSet["name"], modes = modeSet["modes"],
        sigma = sigmaVector[modeSet["modes"], scenario],
        columns, rawInfo, profiled, profiledInfo, retained,
        nuisanceFit, coeffs, alpha1, alpha2, alpha3, alphaCol
      },
      columns = {lnMassColumn[modes], spinColumn[modes]};
      alphaCol = alphaColumn[modes, combo["operator"], combo["polarization"]];
      rawInfo = (alphaCol/sigma).(alphaCol/sigma);
      profiled = weightedFit[alphaCol, columns, sigma];
      profiledInfo = profiled["chi2"];
      retained = If[rawInfo <= 0, 0, profiledInfo/rawInfo];
      nuisanceFit = weightedFit[alphaReference alphaCol, columns, sigma];
      coeffs = nuisanceFit["coefficients"];
      alpha1 = safeLimit[profiledInfo, 1.];
      alpha2 = safeLimit[profiledInfo, 4.];
      alpha3 = safeLimit[profiledInfo, 9.];
      <|
        "mode_set" -> modeSetName,
        "modes" -> StringRiffle[modes, "+"],
        "scenario" -> scenario["name"],
        "sigma_lnf" -> num[scenario["sigma_lnf"]],
        "sigma_lntau" -> num[scenario["sigma_lntau"]],
        "operator" -> combo["operator"],
        "polarization" -> combo["polarization"],
        "raw_information_per_alpha2" -> rawInfo,
        "profiled_information_per_alpha2" -> profiledInfo,
        "retained_information_fraction" -> retained,
        "alpha_1sigma_profiled" -> alpha1,
        "alpha_2sigma_profiled" -> alpha2,
        "alpha_3sigma_profiled" -> alpha3,
        "alpha_reference" -> alphaReference,
        "kerr_only_delta_chi2_at_alpha_reference" ->
          alphaReference^2 profiledInfo,
        "kerr_only_sigma_at_alpha_reference" ->
          alphaReference Sqrt[profiledInfo],
        "best_fit_delta_lnM_for_alpha_reference" ->
          If[Length[coeffs] >= 1, coeffs[[1]], 0],
        "best_fit_delta_chi_for_alpha_reference" ->
          If[Length[coeffs] >= 2, coeffs[[2]], 0]
      |>
    ],
    {modeSet, modeSets},
    {scenario, scenarios},
    {combo, comboRows}
  ],
  2
];

summaryRows = Flatten[
  Table[
    Module[
      {
        rows = Select[
          detectabilityRows,
          #["mode_set"] == modeSet["name"] &&
            #["scenario"] == scenario["name"] &
        ],
        finiteAlpha1, best, median, identifiable, meanRetained
      },
      finiteAlpha1 = Select[Lookup[rows, "alpha_1sigma_profiled"], NumericQ];
      best = If[Length[finiteAlpha1] == 0, Infinity, Min[finiteAlpha1]];
      median = If[Length[finiteAlpha1] == 0, Infinity, Median[finiteAlpha1]];
      identifiable = Count[
        Lookup[rows, "kerr_only_sigma_at_alpha_reference"],
        x_?NumericQ /; x >= 3
      ];
      meanRetained = Mean[Lookup[rows, "retained_information_fraction"]];
      <|
        "mode_set" -> modeSet["name"],
        "modes" -> StringRiffle[modeSet["modes"], "+"],
        "scenario" -> scenario["name"],
        "best_alpha_1sigma_profiled" -> best,
        "median_alpha_1sigma_profiled" -> median,
        "operators_above_3sigma_at_alpha_reference" -> identifiable,
        "operator_count" -> Length[rows],
        "mean_retained_information_fraction" -> meanRetained
      |>
    ],
    {modeSet, modeSets},
    {scenario, scenarios}
  ],
  1
];

confusionModeSet = First[
  Select[modeSets, #["name"] == config["confusion_mode_set"] &]
];
confusionScenario = First[
  Select[scenarios, #["name"] == config["confusion_scenario"] &]
];
confusionModes = confusionModeSet["modes"];
confusionSigma = sigmaVector[confusionModes, confusionScenario];
confusionNuisanceColumns = {
  lnMassColumn[confusionModes],
  spinColumn[confusionModes]
};

confusionRows = Flatten[
  Table[
    Module[
      {
        injCol = alphaColumn[
          confusionModes,
          injection["operator"],
          injection["polarization"]
        ],
        fitCol = alphaColumn[
          confusionModes,
          fit["operator"],
          fit["polarization"]
        ],
        target, columns, fitResult, coeffs
      },
      target = alphaReference injCol;
      columns = Join[confusionNuisanceColumns, {fitCol}];
      fitResult = weightedFit[target, columns, confusionSigma];
      coeffs = fitResult["coefficients"];
      <|
        "mode_set" -> confusionModeSet["name"],
        "scenario" -> confusionScenario["name"],
        "alpha_reference" -> alphaReference,
        "injected_operator" -> injection["operator"],
        "injected_polarization" -> injection["polarization"],
        "fit_operator" -> fit["operator"],
        "fit_polarization" -> fit["polarization"],
        "delta_chi2_after_profile" -> fitResult["chi2"],
        "sigma_after_profile" -> Sqrt[fitResult["chi2"]],
        "best_fit_delta_lnM" -> coeffs[[1]],
        "best_fit_delta_chi" -> coeffs[[2]],
        "best_fit_alpha" -> coeffs[[3]],
        "best_fit_alpha_per_injected_alpha" -> coeffs[[3]]/alphaReference
      |>
    ],
    {injection, comboRows},
    {fit, comboRows}
  ],
  1
];

detectabilityCsvPath = FileNameJoin[
  {outputDir, "synthetic_eft_detectability.csv"}
];
summaryCsvPath = FileNameJoin[
  {outputDir, "synthetic_mode_set_summary.csv"}
];
confusionCsvPath = FileNameJoin[
  {outputDir, "synthetic_fingerprint_confusion.csv"}
];
detectabilityPlotPath = FileNameJoin[
  {outputDir, "synthetic_alpha_1sigma_by_operator.png"}
];
confusionPlotPath = FileNameJoin[
  {outputDir, "synthetic_fingerprint_confusion_heatmap.png"}
];
reportPath = FileNameJoin[
  {outputDir, "ringdown_synthetic_likelihood_report.md"}
];

exportAssociationCSV[detectabilityCsvPath, detectabilityRows];
exportAssociationCSV[summaryCsvPath, summaryRows];
exportAssociationCSV[confusionCsvPath, confusionRows];

plotRows = Select[
  detectabilityRows,
  #["mode_set"] == "220_221" &&
    #["scenario"] == config["confusion_scenario"] &
];
plotLabels = (#["operator"] <> " " <> #["polarization"]) & /@ plotRows;
plotValues = Log10 /@ Lookup[plotRows, "alpha_1sigma_profiled"];

detectabilityPlot = BarChart[
  plotValues,
  ChartLabels -> Placed[Rotate[#, Pi/4] & /@ plotLabels, Below],
  Frame -> True,
  Axes -> False,
  FrameLabel -> {
    None,
    "log10 profiled 1-sigma coupling scale"
  },
  PlotLabel ->
    "Synthetic 220+221 EFT identifiability after M, chi profiling",
  GridLines -> {None, Automatic},
  ImagePadding -> {{70, 40}, {150, 50}},
  ImageSize -> 1100
];
Export[detectabilityPlotPath, detectabilityPlot, ImageResolution -> 144];

comboLabels = (#["operator"] <> " " <> #["polarization"]) & /@ comboRows;
confusionMatrix = Table[
  Module[
    {row = First[
       Select[
         confusionRows,
         #["injected_operator"] == injection["operator"] &&
           #["injected_polarization"] == injection["polarization"] &&
           #["fit_operator"] == fit["operator"] &&
           #["fit_polarization"] == fit["polarization"] &
       ]
     ]},
    Log10[Max[row["delta_chi2_after_profile"], 10^-14]]
  ],
  {injection, comboRows},
  {fit, comboRows}
];

confusionPlot = ArrayPlot[
  confusionMatrix,
  Frame -> True,
  FrameTicks -> {
    Thread[{Range[Length[comboLabels]], comboLabels}],
    Thread[{Range[Length[comboLabels]], Rotate[#, Pi/3] & /@ comboLabels}]
  },
  ColorFunction -> "SolarColors",
  PlotLegends -> BarLegend[
    Automatic,
    LegendLabel -> "log10 Delta chi2"
  ],
  PlotLabel ->
    "Synthetic EFT fingerprint confusion, alpha = " <>
      ToString[alphaReference],
  ImagePadding -> {{170, 140}, {190, 60}},
  ImageSize -> 1200
];
Export[confusionPlotPath, confusionPlot, ImageResolution -> 144];

summaryTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["mode_set"],
        row["scenario"],
        fmt[row["best_alpha_1sigma_profiled"], 5],
        fmt[row["median_alpha_1sigma_profiled"], 5],
        fmt[row["mean_retained_information_fraction"], 4],
        ToString[row["operators_above_3sigma_at_alpha_reference"]] <>
          "/" <> ToString[row["operator_count"]]
      },
      " | "
    ] <> " |"
  ],
  summaryRows
];

bestRows = TakeSmallestBy[
  Select[
    detectabilityRows,
    #["mode_set"] == "220_221" &&
      #["scenario"] == config["confusion_scenario"] &&
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
        fmt[row["kerr_only_sigma_at_alpha_reference"], 3],
        fmt[100 row["best_fit_delta_lnM_for_alpha_reference"], 4],
        fmt[row["best_fit_delta_chi_for_alpha_reference"], 5]
      },
      " | "
    ] <> " |"
  ],
  bestRows
];

report = StringRiffle[
  Join[
    {
      "# Synthetic 220+221 Ringdown Identifiability",
      "",
      "This report is generated by `scripts/wolfram/ringdown_synthetic_likelihood.wl`.",
      "",
      "## Scope",
      "",
      "- This is a synthetic spectral likelihood, not a fit to GW250114 strain.",
      "- The data vector is `{log f, log tau}` for each selected mode.",
      "- Remnant mass and spin are profiled as nuisance parameters.",
      "- One higher-derivative EFT coupling is enabled at a time.",
      "- The Kerr spin dependence uses Berti-Cardoso-Will fits anchored to the numerical Python `qnm` baseline at the GW250114 central spin.",
      "",
      "## Mode-Set Summary",
      "",
      "`operators_above_3sigma_at_alpha_reference` counts how many of the ten operator/polarization fingerprints would give at least a 3-sigma Kerr-only residual for the illustrative `alpha_reference = " <>
        ToString[alphaReference] <> "`.",
      "",
      "| mode set | scenario | best alpha 1sigma | median alpha 1sigma | mean retained information | above 3sigma at alpha_ref |",
      "| --- | --- | ---: | ---: | ---: | ---: |"
    },
    summaryTableRows,
    {
      "",
      "## Most Detectable 220+221 Fingerprints",
      "",
      "These are sorted by the profiled 1-sigma coupling scale for the moderate uncertainty scenario.",
      "",
      "| operator | polarization | alpha 1sigma | retained information | Kerr-only sigma at alpha_ref | best dlnM [%] | best dchi |",
      "| --- | --- | ---: | ---: | ---: | ---: | ---: |"
    },
    bestTableRows,
    {
      "",
      "## Interpretation",
      "",
      "- A single 220 mode has no robust spectral lever arm after profiling over mass and spin: its two observables can be absorbed by two nuisance directions.",
      "- Adding 221 makes the problem over-constrained and leaves finite EFT fingerprints in the residual subspace.",
      "- The retained information fraction is the part of an EFT spectral shift that cannot be mimicked by a small change of remnant mass and spin.",
      "- A matched synthetic EFT injection is recovered exactly in this linearized setup; the confusion matrix asks whether a different operator can mimic the same fingerprint after profiling.",
      "- The absolute `alpha` scales are illustrative until calibrated against real posteriors or strain-level noise.",
      "",
      "## Next Defensible Step",
      "",
      "Use the same projected-frequency machinery to build a toy time-domain likelihood with amplitudes, phases, start time, and a simple colored-noise model. Once that is stable, replace the synthetic widths with GW250114 posterior products or strain-level likelihoods.",
      "",
      "## Generated Files",
      "",
      "- `synthetic_eft_detectability.csv`",
      "- `synthetic_mode_set_summary.csv`",
      "- `synthetic_fingerprint_confusion.csv`",
      "- `synthetic_alpha_1sigma_by_operator.png`",
      "- `synthetic_fingerprint_confusion_heatmap.png`"
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated synthetic ringdown likelihood diagnostics"];
Print["Report: ", reportPath];
Print["Detectability CSV: ", detectabilityCsvPath];
Print["Confusion CSV: ", confusionCsvPath];
Print["mode_set\tscenario\tbest_alpha_1sigma\tmedian_alpha_1sigma\tmean_retained_info"];
Scan[
  Print[
    StringRiffle[
      {
        #["mode_set"],
        #["scenario"],
        fmt[#["best_alpha_1sigma_profiled"], 5],
        fmt[#["median_alpha_1sigma_profiled"], 5],
        fmt[#["mean_retained_information_fraction"], 4]
      },
      "\t"
    ]
  ] &,
  summaryRows
];
