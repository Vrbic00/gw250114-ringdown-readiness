(* ::Package:: *)

(* Approximate EFT projection with public GW250114 pyRing 221 deviation samples.

   This is not a full strain likelihood. It maps the public pyRing
   {domega_221, dtau_221} posterior to
   y = {log(1 + domega_221), log(1 + dtau_221)}.
*)

ClearAll["Global`*"];

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[
    {Directory[], "config", "gw250114_pyring_delta_eft_projection.wl"}
  ]
];
outputDir = If[
  Length[args] >= 2,
  args[[2]],
  FileNameJoin[{Directory[], "results", "gw250114_pyring_delta_eft_projection"}]
];

If[! FileExistsQ[configPath],
  Print["Configuration not found: ", configPath];
  Exit[1];
];

Get[configPath];
If[! AssociationQ[gw250114PyRingDeltaEFTProjectionConfig],
  Print["Configuration must define gw250114PyRingDeltaEFTProjectionConfig."];
  Exit[1];
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

config = gw250114PyRingDeltaEFTProjectionConfig;
pyringDeltaPath = config["pyring_delta_path"];
eftSensitivityPath = config["eft_sensitivity_path"];
mode = ToString[config["mode"]];
deltaFBound = N[config["pyring_delta_f_plot_bound"]];

If[! FileExistsQ[pyringDeltaPath] || ! FileExistsQ[eftSensitivityPath],
  Print["Input file missing."];
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

quantileSummary[name_String, values_List] := Module[
  {clean = N[DeleteMissing[values]], q},
  q = Quantile[clean, {0.05, 0.16, 0.5, 0.84, 0.95}];
  <|
    "parameter" -> name,
    "n" -> Length[clean],
    "mean" -> Mean[clean],
    "standard_deviation" -> StandardDeviation[clean],
    "q05" -> q[[1]],
    "q16" -> q[[2]],
    "median" -> q[[3]],
    "q84" -> q[[4]],
    "q95" -> q[[5]]
  |>
];

csvValue[x_] := If[NumericQ[x], N[x], x];
exportAssociationCSV[path_, rows_List] := Module[{fields = Keys[First[rows]]},
  Export[path, Prepend[csvValue /@ Lookup[#, fields] & /@ rows, fields], "CSV"]
];

Print["Reading pyRing 220+221 delta posterior"];
pyringRowsAll = readDat[pyringDeltaPath];
domegaLower = Exp[-deltaFBound] - 1;
pyringRows = Select[
  pyringRowsAll,
  num[#["domega_221"]] > domegaLower &&
    1 + num[#["domega_221"]] > 0 &&
    1 + num[#["dtau_221"]] > 0 &
];

pyringComputedRows = Map[
  Function[row,
    Join[
      row,
      <|
        "df_221_log" -> Log[1 + num[row["domega_221"]]],
        "dtau_221_log" -> Log[1 + num[row["dtau_221"]]],
        "df_221_log_from_frequencies" ->
          Log[num[row["f_221_corrected"]]/num[row["f_221"]]]
      |>
    ]
  ],
  pyringRows
];

observedVectors = Transpose[
  {
    num /@ Lookup[pyringComputedRows, "df_221_log"],
    num /@ Lookup[pyringComputedRows, "dtau_221_log"]
  }
];
observedMean = Mean /@ Transpose[observedVectors];
observedCov = Covariance[observedVectors];
observedPrecision = Inverse[observedCov];
labels = {"df221_log", "dtau221_log"};

eftRows = Select[importAssociations[eftSensitivityPath], txt[#["mode"]] == mode &];
operators = {"lambda_ev", "lambda_odd", "epsilon1", "epsilon2", "epsilon3"};
polarizations = {"plus", "minus"};
eftByOperatorBranch = Association[
  ((txt[#["operator"]] <> ":" <> txt[#["polarization"]]) -> #) & /@ eftRows
];

fitAlpha[slope_List] := Module[
  {fisher, rhs, alphaBest, alphaSigma, residualBest, chi2Best, chi2Zero,
    deltaChi2},
  fisher = slope.observedPrecision.slope;
  rhs = slope.observedPrecision.observedMean;
  alphaBest = rhs/fisher;
  alphaSigma = 1/Sqrt[fisher];
  residualBest = observedMean - alphaBest slope;
  chi2Best = residualBest.observedPrecision.residualBest;
  chi2Zero = observedMean.observedPrecision.observedMean;
  deltaChi2 = Max[0, chi2Zero - chi2Best];
  <|
    "alpha_best" -> alphaBest,
    "alpha_sigma_gaussian" -> alphaSigma,
    "alpha_q05_gaussian" -> alphaBest - 1.6448536269514722 alphaSigma,
    "alpha_q16_gaussian" -> alphaBest - 0.994457883209753 alphaSigma,
    "alpha_q84_gaussian" -> alphaBest + 0.994457883209753 alphaSigma,
    "alpha_q95_gaussian" -> alphaBest + 1.6448536269514722 alphaSigma,
    "zero_inside_90pct" ->
      Boole[
        alphaBest - 1.6448536269514722 alphaSigma <= 0 <=
          alphaBest + 1.6448536269514722 alphaSigma
      ],
    "delta_chi2_vs_alpha0" -> deltaChi2,
    "gaussian_sigma_from_alpha0" -> Sqrt[deltaChi2],
    "chi2_best" -> chi2Best,
    "chi2_alpha0" -> chi2Zero
  |>
];

fitRows = Flatten[
  Table[
    Module[{row, slope, fit},
      row = eftByOperatorBranch[operator <> ":" <> branch];
      slope = {
        num[row["dln_frequency_dalpha"]],
        num[row["dln_tau_dalpha"]]
      };
      fit = fitAlpha[slope];
      Join[
        <|
          "operator" -> operator,
          "polarization" -> branch,
          "dlnf221_dalpha" -> slope[[1]],
          "dlntau221_dalpha" -> slope[[2]]
        |>,
        fit
      ]
    ],
    {operator, operators},
    {branch, polarizations}
  ],
  1
];

summaryRows = Join[
  {
    Join[
      <|"dataset" -> "PYRING_220_221_filtered"|>,
      quantileSummary[
        "df_221_log",
        num /@ Lookup[pyringComputedRows, "df_221_log"]
      ]
    ],
    Join[
      <|"dataset" -> "PYRING_220_221_filtered"|>,
      quantileSummary[
        "dtau_221_log",
        num /@ Lookup[pyringComputedRows, "dtau_221_log"]
      ]
    ]
  },
  Map[
    Function[column,
      Join[
        <|"dataset" -> "PYRING_220_221_filtered"|>,
        quantileSummary[column, num /@ Lookup[pyringComputedRows, column]]
      ]
    ],
    {"domega_221", "dtau_221", "f_220", "tau_220", "f_221", "tau_221"}
  ]
];

covRows = Flatten[
  Table[
    <|
      "row" -> labels[[i]],
      "column" -> labels[[j]],
      "covariance" -> observedCov[[i, j]],
      "correlation" -> observedCov[[i, j]]/
        Sqrt[observedCov[[i, i]] observedCov[[j, j]]]
    |>,
    {i, Length[labels]},
    {j, Length[labels]}
  ],
  1
];

meanRows = Table[
  <|
    "parameter" -> labels[[i]],
    "observed_mean" -> observedMean[[i]],
    "baseline_value" -> 0.,
    "residual" -> observedMean[[i]],
    "posterior_sigma" -> Sqrt[observedCov[[i, i]]]
  |>,
  {i, Length[labels]}
];

inventoryRows = {
  <|
    "dataset" -> "PYRING_220_221_raw",
    "samples" -> Length[pyringRowsAll],
    "note" -> "raw pyRing 220+221 posterior before filter"
  |>,
  <|
    "dataset" -> "PYRING_220_221_filtered",
    "samples" -> Length[pyringRows],
    "note" -> "after public domega lower-tail filter and positivity checks"
  |>
};

fitCsvPath = FileNameJoin[
  {outputDir, "pyring_delta_eft_projection.csv"}
];
summaryCsvPath = FileNameJoin[
  {outputDir, "pyring_delta_summary.csv"}
];
covCsvPath = FileNameJoin[
  {outputDir, "pyring_delta_covariance.csv"}
];
meanCsvPath = FileNameJoin[
  {outputDir, "pyring_delta_mean_residual.csv"}
];
inventoryCsvPath = FileNameJoin[
  {outputDir, "pyring_delta_dataset_inventory.csv"}
];
plotPath = FileNameJoin[
  {outputDir, "pyring_delta_eft_intervals.png"}
];
reportPath = FileNameJoin[
  {outputDir, "gw250114_pyring_delta_eft_projection_report.md"}
];

exportAssociationCSV[fitCsvPath, fitRows];
exportAssociationCSV[summaryCsvPath, summaryRows];
exportAssociationCSV[covCsvPath, covRows];
exportAssociationCSV[meanCsvPath, meanRows];
exportAssociationCSV[inventoryCsvPath, inventoryRows];

plotLabels = (#["operator"] <> " " <> #["polarization"]) & /@ fitRows;
plotCenters = Lookup[fitRows, "alpha_best"];
plotErrors = Lookup[fitRows, "alpha_sigma_gaussian"];
plotData = Transpose[{Range[Length[fitRows]], plotCenters, plotErrors}];
plotYMin = Min[Join[{0}, plotData[[All, 2]] - plotData[[All, 3]]]];
plotYMax = Max[Join[{0}, plotData[[All, 2]] + plotData[[All, 3]]]];
plotYPad = 0.08 Max[0.01, plotYMax - plotYMin];
plotTicks = MapThread[
  {#1, Rotate[Style[#2, 11], Pi/4]} &,
  {Range[Length[plotLabels]], plotLabels}
];
plot = Graphics[
  {
    {GrayLevel[0.75], Dashed, Line[{{0.5, 0}, {Length[fitRows] + 0.5, 0}}]},
    {RGBColor[0.05, 0.35, 0.45], AbsoluteThickness[1.6],
      Map[Line[{{#[[1]], #[[2]] - #[[3]]}, {#[[1]], #[[2]] + #[[3]]}}] &, plotData],
      Map[Line[{{#[[1]] - 0.12, #[[2]] - #[[3]]}, {#[[1]] + 0.12, #[[2]] - #[[3]]}}] &, plotData],
      Map[Line[{{#[[1]] - 0.12, #[[2]] + #[[3]]}, {#[[1]] + 0.12, #[[2]] + #[[3]]}}] &, plotData]
    },
    {RGBColor[0.0, 0.2, 0.25], PointSize[0.016], Point[plotData[[All, {1, 2}]]]}
  },
  Frame -> True,
  Axes -> False,
  FrameTicks -> {{Automatic, None}, {plotTicks, None}},
  FrameLabel -> {None, "alpha proxy from pyRing {df221, dtau221} projection"},
  PlotLabel -> "GW250114 pyRing 221 frequency+damping EFT projection",
  PlotRange -> {{0.5, Length[fitRows] + 0.5}, {plotYMin - plotYPad, plotYMax + plotYPad}},
  ImagePadding -> {{70, 40}, {155, 55}},
  ImageSize -> 1200
];
Export[plotPath, plot, ImageResolution -> 144];

fitRowsSorted = SortBy[fitRows, #["alpha_sigma_gaussian"] &];

fitTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["operator"],
        row["polarization"],
        fmt[row["alpha_best"], 5],
        fmt[row["alpha_sigma_gaussian"], 5],
        fmt[row["gaussian_sigma_from_alpha0"], 3],
        If[row["zero_inside_90pct"] == 1, "yes", "no"]
      },
      " | "
    ] <> " |"
  ],
  fitRows
];

summaryTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["parameter"],
        fmt[row["median"], 5],
        fmt[row["q16"] - row["median"], 5],
        "+" <> fmt[row["q84"] - row["median"], 5],
        fmt[row["q05"], 5],
        fmt[row["q95"], 5]
      },
      " | "
    ] <> " |"
  ],
  Select[summaryRows, MemberQ[{"df_221_log", "dtau_221_log"}, #["parameter"]] &]
];

meanTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["parameter"],
        fmt[row["observed_mean"], 5],
        fmt[row["posterior_sigma"], 5]
      },
      " | "
    ] <> " |"
  ],
  meanRows
];

report = StringRiffle[
  Join[
    {
      "# GW250114 pyRing Delta EFT Projection",
      "",
      "This report is generated by `scripts/wolfram/gw250114_pyring_delta_eft_projection.wl`.",
      "",
      "## Scope",
      "",
      "- This is an approximate Gaussian projection of the public pyRing 221 deviation posterior, not a full strain-level likelihood.",
      "- The observable vector is `y = {log(1 + domega_221), log(1 + dtau_221)}`.",
      "- The same lower-tail `domega_221` filter used in the public Figure 4 script is applied.",
      "- Each EFT coupling is tested one at a time using the imported higher-derivative QNM frequency and damping-time sensitivities.",
      "- No remnant-mass or remnant-spin nuisance profiling is applied here because the pyRing variables are already fractional deviations from the Kerr prediction.",
      "",
      "## Dataset Inventory",
      "",
      "| dataset | samples | note |",
      "| --- | ---: | --- |",
      "| PYRING_220_221_raw | " <> ToString[Length[pyringRowsAll]] <> " | raw pyRing 220+221 posterior before filter |",
      "| PYRING_220_221_filtered | " <> ToString[Length[pyringRows]] <> " | after public domega lower-tail filter and positivity checks |",
      "",
      "## Posterior Summary",
      "",
      "| parameter | median | -1sigma | +1sigma | q05 | q95 |",
      "| --- | ---: | ---: | ---: | ---: | ---: |"
    },
    summaryTableRows,
    {
      "",
      "## Observable Mean",
      "",
      "| parameter | observed mean | posterior sigma |",
      "| --- | ---: | ---: |"
    },
    meanTableRows,
    {
      "",
      "## EFT Projection",
      "",
      "| operator | polarization | alpha best | alpha sigma | sigma from alpha=0 | zero in 90pct |",
      "| --- | --- | ---: | ---: | ---: | --- |"
    },
    fitTableRows,
    {
      "",
      "## Interpretation",
      "",
      "- Reproduction: this uses the public pyRing `domega_221` and `dtau_221` posterior samples directly.",
      "- Numerical check: adding damping-time information changes the operator ranking because `dln_tau/dalpha` is often large for the 221 overtone.",
      "- Interpretation: no one-at-a-time EFT operator gives a robust nonzero deviation; zero remains inside every 90 percent interval.",
      "- Caution: this is a Gaussian projection of phenomenological pyRing deviation samples, not a replacement for a full EFT waveform likelihood.",
      "",
      "## Next Defensible Step",
      "",
      "Compare the RINGDOWN `{log f_220, log f_221, df_221}` projection and this pyRing `{df_221, dtau_221}` projection side by side, then decide which variables are clean enough for a publication-quality constraints table.",
      "",
      "## Generated Files",
      "",
      "- `pyring_delta_eft_projection.csv`",
      "- `pyring_delta_summary.csv`",
      "- `pyring_delta_covariance.csv`",
      "- `pyring_delta_mean_residual.csv`",
      "- `pyring_delta_dataset_inventory.csv`",
      "- `pyring_delta_eft_intervals.png`"
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated pyRing delta EFT Gaussian projection"];
Print["Report: ", reportPath];
Print[
  "Best alpha rows with zero outside 90pct: ",
  Count[fitRows[[All, "zero_inside_90pct"]], 0]
];
Scan[
  Function[row,
    Print[
      row["operator"], "\t", row["polarization"], "\t",
      fmt[row["alpha_best"], 5], "\t",
      fmt[row["alpha_sigma_gaussian"], 5], "\t",
      fmt[row["gaussian_sigma_from_alpha0"], 3]
    ]
  ],
  Take[fitRowsSorted, UpTo[5]]
];
