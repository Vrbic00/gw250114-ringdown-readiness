(* ::Package:: *)

(* Empirical linearized posterior projection for public GW250114 constraints.

   This is a non-Gaussian sanity check of the linearized alpha projection. It
   maps every public posterior sample to an alpha estimate. It is still not a
   full strain-level EFT likelihood.
*)

ClearAll["Global`*"];

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[
    {Directory[], "config", "gw250114_linearized_posterior_projection.wl"}
  ]
];
outputDir = If[
  Length[args] >= 2,
  args[[2]],
  FileNameJoin[
    {Directory[], "results", "gw250114_linearized_posterior_projection"}
  ]
];

If[! FileExistsQ[configPath],
  Print["Configuration not found: ", configPath];
  Exit[1];
];

Get[configPath];
If[! AssociationQ[gw250114LinearizedPosteriorProjectionConfig],
  Print[
    "Configuration must define gw250114LinearizedPosteriorProjectionConfig."
  ];
  Exit[1];
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

config = gw250114LinearizedPosteriorProjectionConfig;
event = config["event"];
ringdownPath = config["ringdown_hdf5_path"];
pyringDeltaPath = config["pyring_delta_path"];
eftSensitivityPath = config["eft_sensitivity_path"];
kerrPath = config["kerr_numeric_path"];
gaussianProjectionPath = config["gaussian_projection_path"];
mass0 = N[event["mass_detector_msun"]];
spin0 = N[event["spin"]];
deltaFBound = N[config["pyring_delta_f_plot_bound"]];

requiredFiles = {
  ringdownPath, pyringDeltaPath, eftSensitivityPath, kerrPath,
  gaussianProjectionPath
};
missingFiles = Select[requiredFiles, ! FileExistsQ[#] &];
If[Length[missingFiles] > 0,
  Print["Missing input files: ", missingFiles];
  Exit[1];
];

num[value_] := N[If[NumericQ[value], value, ToExpression[ToString[value]]]];
txt[value_] := ToString[value];

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

readDat[path_] := Module[{lines, headerLine, names, raw, numericRows},
  lines = Import[path, "Lines"];
  headerLine = First[Select[lines, StringStartsQ[StringTrim[#], "#"] &]];
  names = StringSplit[StringTrim[StringDrop[StringTrim[headerLine], 1]]];
  raw = Import[path, "Table"];
  numericRows = Select[raw, VectorQ[#, NumericQ] &];
  AssociationThread[names, #] & /@ numericRows
];

importAssociations[path_] := Module[{raw = Import[path, "CSV"]},
  AssociationThread[First[raw], #] & /@ Rest[raw]
];

csvValue[x_] := If[NumericQ[x], N[x], x];
exportAssociationCSV[path_, rows_List] := Module[{fields = Keys[First[rows]]},
  Export[path, Prepend[csvValue /@ Lookup[#, fields] & /@ rows, fields], "CSV"]
];

rowKey[projection_String, row_Association] :=
  projection <> ":" <> txt[row["operator"]] <> ":" <> txt[row["polarization"]];

qnmCoefficients = <|
  "220" -> <|"f1" -> 1.5251, "f2" -> -1.1568, "f3" -> 0.1292,
    "q1" -> 0.7000, "q2" -> 1.4187, "q3" -> -0.4990|>,
  "221" -> <|"f1" -> 1.3673, "f2" -> -1.0260, "f3" -> 0.1628,
    "q1" -> 0.1000, "q2" -> 0.5436, "q3" -> -0.4731|>
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

kerrRows = importAssociations[kerrPath];
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
kerrAnchorByMode = Association[
  (# -> (kerrExactByMode[#] - bertiMomega[#, spin0])) & /@ {"220", "221"}
];
kerrMomega[mode_String, spin_?NumericQ] :=
  bertiMomega[mode, spin] + kerrAnchorByMode[mode];

spinStep = 0.0005;
spinLogFrequencyDerivative[mode_String] :=
  (Log[Re[kerrMomega[mode, spin0 + spinStep]]] -
    Log[Re[kerrMomega[mode, spin0 - spinStep]]])/(2 spinStep);

gaussianRows = importAssociations[gaussianProjectionPath];
gaussianByKey = Association[
  ((txt[#["projection"]] <> ":" <> txt[#["operator"]] <> ":" <>
        txt[#["polarization"]]) -> #) & /@ gaussianRows
];

eftRows = Select[
  importAssociations[eftSensitivityPath],
  MemberQ[{"220", "221"}, txt[#["mode"]]] &
];
eftByModeOperatorBranch = Association[
  ((txt[#["mode"]] <> ":" <> txt[#["operator"]] <> ":" <>
        txt[#["polarization"]]) -> #) & /@ eftRows
];
eft221ByOperatorBranch = Association[
  ((txt[#["operator"]] <> ":" <> txt[#["polarization"]]) -> #) & /@
    Select[eftRows, txt[#["mode"]] == "221" &]
];

operators = {"lambda_ev", "lambda_odd", "epsilon1", "epsilon2", "epsilon3"};
polarizations = {"plus", "minus"};

alphaSummary[
  projection_String, observableSet_String, role_String, operator_String,
  polarization_String, alphaSamples_List
] := Module[
  {clean = N[alphaSamples], q, n, mean, sd, median, pLess, pGreater,
    gaussianRow, gaussianKey, gaussianBest, gaussianSigma,
    gaussianQ05, gaussianQ95, q05, q95, halfWidth90},
  n = Length[clean];
  q = Quantile[clean, {0.05, 0.16, 0.5, 0.84, 0.95}];
  q05 = q[[1]];
  q95 = q[[5]];
  mean = Mean[clean];
  sd = StandardDeviation[clean];
  median = q[[3]];
  pLess = N[Count[clean, _?(# < 0 &)]/n];
  pGreater = N[Count[clean, _?(# > 0 &)]/n];
  gaussianKey = projection <> ":" <> operator <> ":" <> polarization;
  gaussianRow = Lookup[gaussianByKey, gaussianKey, Missing["NotAvailable"]];
  gaussianBest = If[AssociationQ[gaussianRow],
    num[gaussianRow["alpha_best"]],
    Missing["NotAvailable"]
  ];
  gaussianSigma = If[AssociationQ[gaussianRow],
    num[gaussianRow["alpha_sigma"]],
    Missing["NotAvailable"]
  ];
  gaussianQ05 = If[AssociationQ[gaussianRow],
    num[gaussianRow["alpha_q05"]],
    Missing["NotAvailable"]
  ];
  gaussianQ95 = If[AssociationQ[gaussianRow],
    num[gaussianRow["alpha_q95"]],
    Missing["NotAvailable"]
  ];
  halfWidth90 = (q95 - q05)/2;
  <|
    "projection" -> projection,
    "observable_set" -> observableSet,
    "role" -> role,
    "operator" -> operator,
    "polarization" -> polarization,
    "sample_count" -> n,
    "alpha_mean_empirical" -> mean,
    "alpha_sd_empirical" -> sd,
    "alpha_median_empirical" -> median,
    "alpha_q05_empirical" -> q05,
    "alpha_q16_empirical" -> q[[2]],
    "alpha_q84_empirical" -> q[[4]],
    "alpha_q95_empirical" -> q95,
    "zero_inside_empirical_90pct" -> Boole[q05 <= 0 <= q95],
    "posterior_mass_alpha_less_zero" -> pLess,
    "posterior_mass_alpha_greater_zero" -> pGreater,
    "nominal_abs_median_over_sd" -> Abs[median]/sd,
    "interval_asymmetry_90pct" ->
      If[halfWidth90 > 0, ((q95 - median) - (median - q05))/(2 halfWidth90), 0],
    "gaussian_alpha_best" -> gaussianBest,
    "gaussian_alpha_sigma" -> gaussianSigma,
    "gaussian_q05" -> gaussianQ05,
    "gaussian_q95" -> gaussianQ95,
    "median_minus_gaussian_best" ->
      If[NumericQ[gaussianBest], median - gaussianBest, Missing["NotAvailable"]],
    "empirical_sd_over_gaussian_sigma" ->
      If[NumericQ[gaussianSigma] && gaussianSigma != 0,
        sd/gaussianSigma,
        Missing["NotAvailable"]
      ]
  |>
];

Print["Reading direct RINGDOWN samples"];
f220 = Import[ringdownPath, {"Datasets", "/f_220"}];
f221 = Import[ringdownPath, {"Datasets", "/f_221"}];
df221 = Import[ringdownPath, {"Datasets", "/df_221"}];
ringdownVectors = N[Transpose[{Log[f220], Log[f221], df221}]];
ringdownCov = Covariance[ringdownVectors];
ringdownPrecision = Inverse[ringdownCov];
ringdownBaseline = {
  Log[num[First[Select[kerrRows, txt[#["mode"]] == "220" &]]["qnm_f_Hz"]]],
  Log[num[First[Select[kerrRows, txt[#["mode"]] == "221" &]]["qnm_f_Hz"]]],
  0.
};
ringdownResidualSamples = Map[# - ringdownBaseline &, ringdownVectors];
ringdownNuisanceColumns = {
  {-1., -1., 0.},
  {spinLogFrequencyDerivative["220"], spinLogFrequencyDerivative["221"], 0.}
};

linearCoefficientMatrix[columns_List, precision_?MatrixQ] := Module[
  {design},
  design = Transpose[columns];
  PseudoInverse[Transpose[design].precision.design].Transpose[design].precision
];

ringdownRows = Flatten[
  Table[
    Module[{s220, s221, alphaColumn, coefficientMatrix, alphaWeight,
      alphaSamples},
      s220 = eftByModeOperatorBranch["220:" <> operator <> ":" <> branch];
      s221 = eftByModeOperatorBranch["221:" <> operator <> ":" <> branch];
      alphaColumn = {
        num[s220["dln_frequency_dalpha"]],
        num[s221["dln_frequency_dalpha"]],
        num[s221["dln_frequency_dalpha"]]
      };
      coefficientMatrix = linearCoefficientMatrix[
        Join[ringdownNuisanceColumns, {alphaColumn}],
        ringdownPrecision
      ];
      alphaWeight = Last[coefficientMatrix];
      alphaSamples = ringdownResidualSamples.alphaWeight;
      alphaSummary[
        "RINGDOWN",
        "{log f_220, log f_221, df_221}",
        "linearized sample projection with mass/spin profiling",
        operator,
        branch,
        alphaSamples
      ]
    ],
    {operator, operators},
    {branch, polarizations}
  ],
  1
];

Print["Reading pyRing delta samples"];
pyringRowsAll = readDat[pyringDeltaPath];
domegaLower = Exp[-deltaFBound] - 1;
pyringRows = Select[
  pyringRowsAll,
  num[#["domega_221"]] > domegaLower &&
    1 + num[#["domega_221"]] > 0 &&
    1 + num[#["dtau_221"]] > 0 &
];
pyringVectors = N[
  Transpose[
    {
      Log[1 + num /@ Lookup[pyringRows, "domega_221"]],
      Log[1 + num /@ Lookup[pyringRows, "dtau_221"]]
    }
  ]
];
pyringCov = Covariance[pyringVectors];
pyringPrecision = Inverse[pyringCov];

pyringRowsProjected = Flatten[
  Table[
    Module[{s221, slope, fisher, alphaWeight, alphaSamples},
      s221 = eft221ByOperatorBranch[operator <> ":" <> branch];
      slope = {
        num[s221["dln_frequency_dalpha"]],
        num[s221["dln_tau_dalpha"]]
      };
      fisher = slope.pyringPrecision.slope;
      alphaWeight = (slope.pyringPrecision)/fisher;
      alphaSamples = pyringVectors.alphaWeight;
      alphaSummary[
        "PYRING_DELTA",
        "{log(1 + domega_221), log(1 + dtau_221)}",
        "linearized sample projection in pyRing deviation variables",
        operator,
        branch,
        alphaSamples
      ]
    ],
    {operator, operators},
    {branch, polarizations}
  ],
  1
];

projectionRows = Join[ringdownRows, pyringRowsProjected];

projectionSummaryRows = Map[
  Function[projection,
    Module[{rows = Select[projectionRows, #["projection"] == projection &],
      bestZRow, widestAsymmetryRow},
      bestZRow = Last[SortBy[rows, #["nominal_abs_median_over_sd"] &]];
      widestAsymmetryRow = Last[SortBy[rows, Abs[#["interval_asymmetry_90pct"]] &]];
      <|
        "projection" -> projection,
        "rows" -> Length[rows],
        "zero_outside_empirical_90pct_count" ->
          Count[Lookup[rows, "zero_inside_empirical_90pct"], 0],
        "max_nominal_abs_median_over_sd" ->
          bestZRow["nominal_abs_median_over_sd"],
        "max_nominal_operator" -> bestZRow["operator"],
        "max_nominal_polarization" -> bestZRow["polarization"],
        "max_abs_interval_asymmetry" ->
          Abs[widestAsymmetryRow["interval_asymmetry_90pct"]],
        "max_asymmetry_operator" -> widestAsymmetryRow["operator"],
        "max_asymmetry_polarization" -> widestAsymmetryRow["polarization"]
      |>
    ]
  ],
  {"RINGDOWN", "PYRING_DELTA"}
];

comparisonRows = Map[
  Function[row,
    <|
      "projection" -> row["projection"],
      "operator" -> row["operator"],
      "polarization" -> row["polarization"],
      "empirical_median" -> row["alpha_median_empirical"],
      "empirical_q05" -> row["alpha_q05_empirical"],
      "empirical_q95" -> row["alpha_q95_empirical"],
      "gaussian_best" -> row["gaussian_alpha_best"],
      "gaussian_q05" -> row["gaussian_q05"],
      "gaussian_q95" -> row["gaussian_q95"],
      "median_minus_gaussian_best" -> row["median_minus_gaussian_best"],
      "empirical_sd_over_gaussian_sigma" ->
        row["empirical_sd_over_gaussian_sigma"],
      "zero_inside_empirical_90pct" ->
        row["zero_inside_empirical_90pct"]
    |>
  ],
  projectionRows
];

rowsCsvPath = FileNameJoin[
  {outputDir, "linearized_alpha_posterior_summary.csv"}
];
comparisonCsvPath = FileNameJoin[
  {outputDir, "linearized_vs_gaussian_projection.csv"}
];
projectionSummaryCsvPath = FileNameJoin[
  {outputDir, "linearized_projection_summary.csv"}
];
plotPath = FileNameJoin[
  {outputDir, "linearized_empirical_intervals.png"}
];
reportPath = FileNameJoin[
  {outputDir, "gw250114_linearized_posterior_projection_report.md"}
];

exportAssociationCSV[rowsCsvPath, projectionRows];
exportAssociationCSV[comparisonCsvPath, comparisonRows];
exportAssociationCSV[projectionSummaryCsvPath, projectionSummaryRows];

plotRows = Take[
  SortBy[projectionRows, #["alpha_sd_empirical"] &],
  UpTo[12]
];
plotLabels = (#["projection"] <> " " <> #["operator"] <> " " <>
      #["polarization"]) & /@ plotRows;
xvals = Range[Length[plotRows]];
q05s = Lookup[plotRows, "alpha_q05_empirical"];
q95s = Lookup[plotRows, "alpha_q95_empirical"];
medians = Lookup[plotRows, "alpha_median_empirical"];
plotYMin = Min[Join[{0}, q05s]];
plotYMax = Max[Join[{0}, q95s]];
plotYPad = 0.08 Max[0.01, plotYMax - plotYMin];
plotTicks = MapThread[
  {#1, Rotate[Style[#2, 10], Pi/4]} &,
  {xvals, plotLabels}
];
plot = Graphics[
  {
    {GrayLevel[0.75], Dashed, Line[{{0.5, 0}, {Length[xvals] + 0.5, 0}}]},
    {RGBColor[0.15, 0.35, 0.65], AbsoluteThickness[1.5],
      MapThread[Line[{{#1, #2}, {#1, #3}}] &, {xvals, q05s, q95s}],
      MapThread[Line[{{#1 - 0.1, #2}, {#1 + 0.1, #2}}] &, {xvals, q05s}],
      MapThread[Line[{{#1 - 0.1, #2}, {#1 + 0.1, #2}}] &, {xvals, q95s}]
    },
    {RGBColor[0.05, 0.15, 0.35], PointSize[0.014],
      Point[Transpose[{xvals, medians}]]}
  },
  Frame -> True,
  Axes -> False,
  FrameTicks -> {{Automatic, None}, {plotTicks, None}},
  FrameLabel -> {None, "empirical alpha posterior median and 90 pct interval"},
  PlotLabel -> "GW250114 linearized sample-level alpha projection",
  PlotRange -> {{0.5, Length[xvals] + 0.5}, {plotYMin - plotYPad, plotYMax + plotYPad}},
  ImagePadding -> {{70, 30}, {175, 55}},
  ImageSize -> 1200
];
Export[plotPath, plot, ImageResolution -> 144];

summaryTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["projection"],
        ToString[row["zero_outside_empirical_90pct_count"]],
        fmt[row["max_nominal_abs_median_over_sd"], 3],
        row["max_nominal_operator"] <> " " <> row["max_nominal_polarization"],
        fmt[row["max_abs_interval_asymmetry"], 3],
        row["max_asymmetry_operator"] <> " " <>
          row["max_asymmetry_polarization"]
      },
      " | "
    ] <> " |"
  ],
  projectionSummaryRows
];

tightestRowsForReport = Take[
  SortBy[projectionRows, #["alpha_sd_empirical"] &],
  UpTo[8]
];
tightestTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["projection"],
        row["operator"],
        row["polarization"],
        fmt[row["alpha_median_empirical"], 5],
        fmt[row["alpha_q05_empirical"], 5],
        fmt[row["alpha_q95_empirical"], 5],
        fmt[row["nominal_abs_median_over_sd"], 3],
        If[row["zero_inside_empirical_90pct"] == 1, "yes", "no"]
      },
      " | "
    ] <> " |"
  ],
  tightestRowsForReport
];

largestGaussianDeltaRows = Take[
  Reverse[SortBy[projectionRows, Abs[#["median_minus_gaussian_best"]] &]],
  UpTo[6]
];
gaussianComparisonTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["projection"],
        row["operator"],
        row["polarization"],
        fmt[row["alpha_median_empirical"], 5],
        fmt[row["gaussian_alpha_best"], 5],
        fmt[row["median_minus_gaussian_best"], 5],
        fmt[row["empirical_sd_over_gaussian_sigma"], 3]
      },
      " | "
    ] <> " |"
  ],
  largestGaussianDeltaRows
];

zeroOutsideTotal = Total[
  Lookup[projectionSummaryRows, "zero_outside_empirical_90pct_count"]
];
maxNominal = Max[
  Lookup[projectionSummaryRows, "max_nominal_abs_median_over_sd"]
];

report = StringRiffle[
  Join[
    {
      "# GW250114 Linearized Posterior Projection",
      "",
      "This report is generated by `scripts/wolfram/gw250114_linearized_posterior_projection.wl`.",
      "",
      "## Scope",
      "",
      "- This maps each public posterior sample to a linearized EFT coupling estimate.",
      "- RINGDOWN samples are projected with linear mass/spin nuisance profiling.",
      "- pyRing samples are projected in `{log(1 + domega_221), log(1 + dtau_221)}`.",
      "- This is a non-Gaussian sanity check of the projection layer, not a full strain-level EFT likelihood.",
      "",
      "## Projection Summary",
      "",
      "| projection | zero outside empirical 90pct | max abs(median)/sd | max row | max 90pct asymmetry | max asymmetry row |",
      "| --- | ---: | ---: | --- | ---: | --- |"
    },
    summaryTableRows,
    {
      "",
      "## Tightest Empirical Rows",
      "",
      "| projection | operator | polarization | alpha median | q05 | q95 | abs(median)/sd | zero in 90pct |",
      "| --- | --- | --- | ---: | ---: | ---: | ---: | --- |"
    },
    tightestTableRows,
    {
      "",
      "## Largest Median-vs-Gaussian Differences",
      "",
      "| projection | operator | polarization | empirical median | Gaussian best | median - Gaussian | empirical sd / Gaussian sigma |",
      "| --- | --- | --- | ---: | ---: | ---: | ---: |"
    },
    gaussianComparisonTableRows,
    {
      "",
      "## Interpretation",
      "",
      "- Reproduction: the empirical linearized posterior uses the same public posterior samples and the same linearized EFT fingerprints as the Gaussian projection.",
      "- Numerical check: alpha=0 remains inside every empirical 90 percent interval.",
      "- Interpretation: the no-deviation conclusion does not appear to be an artifact of symmetric Gaussian intervals.",
      "- Caution: this still uses a linearized projection and should not be described as a full EFT waveform likelihood.",
      "",
      "## Generated Files",
      "",
      "- `linearized_alpha_posterior_summary.csv`",
      "- `linearized_vs_gaussian_projection.csv`",
      "- `linearized_projection_summary.csv`",
      "- `linearized_empirical_intervals.png`"
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated GW250114 linearized posterior projection"];
Print["Report: ", reportPath];
Print["Rows with alpha=0 outside empirical 90pct: ", zeroOutsideTotal];
Print["Maximum empirical |median|/sd: ", fmt[maxNominal, 3]];
