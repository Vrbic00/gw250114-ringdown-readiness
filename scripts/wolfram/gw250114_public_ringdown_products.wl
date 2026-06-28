(* ::Package:: *)

(* Summarize public GW250114 ringdown/post-merger data products.

   Usage:
     wolframscript -file scripts/wolfram/gw250114_public_ringdown_products.wl
*)

ClearAll["Global`*"];

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[
    {Directory[], "config", "gw250114_public_ringdown_products.wl"}
  ]
];
outputDir = If[
  Length[args] >= 2,
  args[[2]],
  FileNameJoin[{Directory[], "results", "gw250114_public_ringdown_products"}]
];

If[! FileExistsQ[configPath],
  Print["Configuration not found: ", configPath];
  Exit[1];
];

Get[configPath];
If[! AssociationQ[gw250114PublicRingdownConfig],
  Print["Configuration must define gw250114PublicRingdownConfig."];
  Exit[1];
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

config = gw250114PublicRingdownConfig;
releaseRoot = config["release_root"];
pyring220Path = FileNameJoin[{releaseRoot, config["pyring_220_posterior"]}];
pyringDeltaPath = FileNameJoin[
  {releaseRoot, config["pyring_220_221_delta_posterior"]}
];
ringdownDeltaPath = FileNameJoin[
  {releaseRoot, config["ringdown_220_221_delta_posterior"]}
];
eftSensitivityPath = config["eft_sensitivity_path"];
deltaFBound = N[config["pyring_delta_f_plot_bound"]];
alphaProxyMode = config["alpha_proxy_mode"];

requiredFiles = {
  pyring220Path, pyringDeltaPath, ringdownDeltaPath, eftSensitivityPath
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

summaryRowsFor[prefix_String, assocRows_List, columns_List] := Map[
  Function[column,
    Join[
      <|"dataset" -> prefix|>,
      quantileSummary[column, num /@ Lookup[assocRows, column]]
    ]
  ],
  columns
];

summaryRowsForArrays[prefix_String, data_Association] := Map[
  Function[key,
    Join[<|"dataset" -> prefix|>, quantileSummary[key, data[key]]]
  ],
  Keys[data]
];

csvValue[x_] := If[NumericQ[x], N[x], x];

exportAssociationCSV[path_, rows_List] := Module[{fields = Keys[First[rows]]},
  Export[path, Prepend[csvValue /@ Lookup[#, fields] & /@ rows, fields], "CSV"]
];

Print["Reading pyRing 220 posterior"];
pyring220Rows = readDat[pyring220Path];

Print["Reading pyRing 220+221 delta posterior with QNM frequencies"];
pyringDeltaRowsAll = readDat[pyringDeltaPath];
domegaLower = Exp[-deltaFBound] - 1;
pyringDeltaRows = Select[
  pyringDeltaRowsAll,
  num[#["domega_221"]] > domegaLower &&
    1 + num[#["domega_221"]] > 0 &&
    1 + num[#["dtau_221"]] > 0 &
];

pyringDeltaRowsComputed = Map[
  Function[row,
    Join[
      row,
      <|
        "df_221_log" -> Log[1 + num[row["domega_221"]]],
        "dtau_221_log" -> Log[1 + num[row["dtau_221"]]],
        "f221_corrected_over_kerr" ->
          num[row["f_221_corrected"]]/num[row["f_221"]]
      |>
    ]
  ],
  pyringDeltaRows
];

Print["Reading direct RINGDOWN 220+221 delta posterior"];
ringdownData = <|
  "df_221" -> Import[ringdownDeltaPath, {"Datasets", "/df_221"}],
  "f_220" -> Import[ringdownDeltaPath, {"Datasets", "/f_220"}],
  "f_221" -> Import[ringdownDeltaPath, {"Datasets", "/f_221"}]
|>;
ringdownData = Join[
  ringdownData,
  <|"f_221_inferred_kerr" -> ringdownData["f_221"]/Exp[ringdownData["df_221"]]|>
];

eftSensitivityRows = Select[
  importAssociations[eftSensitivityPath],
  txt[#["mode"]] == alphaProxyMode &
];

pyring220Summary = summaryRowsFor[
  "PYRING_220",
  pyring220Rows,
  {"Mf", "af", "cosiota", "logL"}
];

pyringDeltaSummary = summaryRowsFor[
  "PYRING_220_221_filtered",
  pyringDeltaRowsComputed,
  {
    "Mf", "af", "domega_221", "dtau_221", "df_221_log",
    "dtau_221_log", "f_220", "tau_220", "f_221", "tau_221",
    "f_221_corrected", "f221_corrected_over_kerr", "logL"
  }
];

ringdownSummary = summaryRowsForArrays[
  "RINGDOWN_220_221",
  ringdownData
];

dfProxySamples = <|
  "RINGDOWN_220_221" -> ringdownData["df_221"],
  "PYRING_220_221_filtered" -> Lookup[pyringDeltaRowsComputed, "df_221_log"]
|>;

alphaProxyRows = Flatten[
  Table[
    Module[
      {
        sensitivity = num[sens["dln_frequency_dalpha"]],
        alphaSamples, summary
      },
      alphaSamples = dfProxySamples[dataset]/sensitivity;
      summary = quantileSummary["alpha_proxy", alphaSamples];
      <|
        "dataset" -> dataset,
        "operator" -> sens["operator"],
        "polarization" -> sens["polarization"],
        "dln_f221_dalpha" -> sensitivity,
        "n" -> summary["n"],
        "mean" -> summary["mean"],
        "standard_deviation" -> summary["standard_deviation"],
        "q05" -> summary["q05"],
        "q16" -> summary["q16"],
        "median" -> summary["median"],
        "q84" -> summary["q84"],
        "q95" -> summary["q95"],
        "zero_inside_90pct" ->
          Boole[summary["q05"] <= 0 <= summary["q95"]]
      |>
    ],
    {dataset, Keys[dfProxySamples]},
    {sens, eftSensitivityRows}
  ],
  1
];

allSummaryRows = Join[pyring220Summary, pyringDeltaSummary, ringdownSummary];

summaryCsvPath = FileNameJoin[
  {outputDir, "public_ringdown_product_summary.csv"}
];
alphaProxyCsvPath = FileNameJoin[
  {outputDir, "df221_to_eft_alpha_proxy.csv"}
];
datasetCsvPath = FileNameJoin[
  {outputDir, "public_ringdown_dataset_inventory.csv"}
];
plotPath = FileNameJoin[{outputDir, "df221_public_posteriors.png"}];
reportPath = FileNameJoin[
  {outputDir, "gw250114_public_ringdown_products_report.md"}
];

datasetInventoryRows = {
  <|
    "dataset" -> "PYRING_220",
    "source_file" -> config["pyring_220_posterior"],
    "samples" -> Length[pyring220Rows],
    "note" -> "pyRing 220 posterior"
  |>,
  <|
    "dataset" -> "PYRING_220_221_filtered",
    "source_file" -> config["pyring_220_221_delta_posterior"],
    "samples" -> Length[pyringDeltaRows],
    "note" -> "pyRing 220+221 posterior after the published domega_221 lower-tail filter"
  |>,
  <|
    "dataset" -> "PYRING_220_221_raw",
    "source_file" -> config["pyring_220_221_delta_posterior"],
    "samples" -> Length[pyringDeltaRowsAll],
    "note" -> "raw pyRing 220+221 posterior before filter"
  |>,
  <|
    "dataset" -> "RINGDOWN_220_221",
    "source_file" -> config["ringdown_220_221_delta_posterior"],
    "samples" -> Length[ringdownData["df_221"]],
    "note" -> "direct RINGDOWN product with f_220, f_221, and df_221"
  |>
};

exportAssociationCSV[summaryCsvPath, allSummaryRows];
exportAssociationCSV[alphaProxyCsvPath, alphaProxyRows];
exportAssociationCSV[datasetCsvPath, datasetInventoryRows];

hist = Histogram[
  {
    ringdownData["df_221"],
    Lookup[pyringDeltaRowsComputed, "df_221_log"]
  },
  {0.04},
  "PDF",
  ChartLegends -> {"RINGDOWN", "PYRING filtered"},
  Frame -> True,
  Axes -> False,
  FrameLabel -> {"df_221 = log(f_measured/f_Kerr)", "density"},
  PlotLabel -> "GW250114 public ringdown df_221 posteriors",
  ChartStyle -> {Directive[Blue, Opacity[0.45]], Directive[Orange, Opacity[0.45]]},
  PlotRange -> {{-deltaFBound, deltaFBound}, All},
  ImageSize -> 1000
];
Export[plotPath, hist, ImageResolution -> 144];

summaryByDatasetParameter = Association[
  ((#["dataset"] <> ":" <> #["parameter"]) -> #) & /@ allSummaryRows
];

alphaProxyBestRows = TakeSmallestBy[
  Select[alphaProxyRows, #["dataset"] == "RINGDOWN_220_221" &],
  Abs[#["median"]] &,
  UpTo[5]
];

tableLine[row_, cols_] :=
  "| " <> StringRiffle[
    (If[NumericQ[row[#]], fmt[row[#], 5], ToString[row[#]]] & /@ cols),
    " | "
  ] <> " |";

inventoryRows = tableLine[
    #,
    {"dataset", "samples", "source_file", "note"}
  ] & /@ datasetInventoryRows;

dfSummaryRows = Map[
  Function[pair,
    Module[
      {
        dataset = pair[[1]], parameter = pair[[2]], row
      },
      row = summaryByDatasetParameter[pair[[1]] <> ":" <> pair[[2]]];
      "| " <> StringRiffle[
        {
          dataset,
          parameter,
          fmt[row["median"], 5],
          "-" <> fmt[row["median"] - row["q16"], 5],
          "+" <> fmt[row["q84"] - row["median"], 5],
          fmt[row["q05"], 5],
          fmt[row["q95"], 5]
        },
        " | "
      ] <> " |"
    ]
  ],
  {
    {"RINGDOWN_220_221", "df_221"},
    {"RINGDOWN_220_221", "f_220"},
    {"RINGDOWN_220_221", "f_221"},
    {"RINGDOWN_220_221", "f_221_inferred_kerr"},
    {"PYRING_220_221_filtered", "df_221_log"},
    {"PYRING_220_221_filtered", "domega_221"},
    {"PYRING_220_221_filtered", "dtau_221"},
    {"PYRING_220_221_filtered", "f_220"},
    {"PYRING_220_221_filtered", "f_221_corrected"}
  }
];

alphaProxyTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["operator"],
        row["polarization"],
        fmt[row["dln_f221_dalpha"], 4],
        fmt[row["median"], 5],
        "-" <> fmt[row["median"] - row["q16"], 5],
        "+" <> fmt[row["q84"] - row["median"], 5],
        If[row["zero_inside_90pct"] == 1, "yes", "no"]
      },
      " | "
    ] <> " |"
  ],
  alphaProxyBestRows
];

report = StringRiffle[
  Join[
    {
      "# GW250114 Public Ringdown Products",
      "",
      "This report is generated by `scripts/wolfram/gw250114_public_ringdown_products.wl`.",
      "",
      "## Scope",
      "",
      "- These files come from the public GW250114 Zenodo tarball.",
      "- The direct RINGDOWN product contains `f_220`, `f_221`, and `df_221` for 120,000 samples.",
      "- The pyRing delta product contains `domega_221`, `dtau_221`, and converted QNM frequencies.",
      "- Following the published Figure 4 script, the pyRing frequency deviation is `df_221 = log(f_221_corrected/f_221) = log(1 + domega_221)` after removing the low-`domega_221` tail used in that script.",
      "",
      "## Dataset Inventory",
      "",
      "| dataset | samples | source file | note |",
      "| --- | ---: | --- | --- |"
    },
    inventoryRows,
    {
      "",
      "## Main Posterior Summaries",
      "",
      "| dataset | parameter | median | -1sigma | +1sigma | q05 | q95 |",
      "| --- | --- | ---: | ---: | ---: | ---: | ---: |"
    },
    dfSummaryRows,
    {
      "",
      "## One-Dimensional EFT Proxy",
      "",
      "The table maps the direct RINGDOWN `df_221` posterior to a first-order coupling proxy using `df_221 approximately alpha * dln f_221/dalpha`. This ignores damping time, mode correlations, polarization mixing, amplitudes, and first-order EFT validity limits.",
      "",
      "| operator | polarization | dln f221/dalpha | alpha median | -1sigma | +1sigma | zero in 90pct |",
      "| --- | --- | ---: | ---: | ---: | ---: | --- |"
    },
    alphaProxyTableRows,
    {
      "",
      "## Interpretation",
      "",
      "- Reproduction: the public ringdown/post-merger data products are now locally available and summarized.",
      "- Numerical check: the `df_221` variable in the public script is a logarithmic frequency deviation, while pyRing stores the corresponding linear `domega_221`.",
      "- Interpretation: these data are much closer to the actual ringdown test we need than the full IMR posterior used earlier.",
      "- Caution: a physically defensible EFT constraint still needs a likelihood-level map, not only a 1D frequency-deviation proxy.",
      "",
      "## Next Defensible Step",
      "",
      "Build a ringdown-posterior likelihood in the measured variables `f_220`, `f_221`, and `df_221`, then project each EFT fingerprint through the same variables instead of using synthetic `{log f, log tau}` widths.",
      "",
      "## Generated Files",
      "",
      "- `public_ringdown_product_summary.csv`",
      "- `public_ringdown_dataset_inventory.csv`",
      "- `df221_to_eft_alpha_proxy.csv`",
      "- `df221_public_posteriors.png`"
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated public ringdown product summary"];
Print["Report: ", reportPath];
Print["RINGDOWN df_221 median: ",
  fmt[summaryByDatasetParameter["RINGDOWN_220_221:df_221"]["median"], 5]];
Print["PYRING filtered df_221 median: ",
  fmt[summaryByDatasetParameter["PYRING_220_221_filtered:df_221_log"]["median"], 5]];
Print["Samples RINGDOWN/PYRING filtered: ",
  Length[ringdownData["df_221"]], " / ", Length[pyringDeltaRows]];
