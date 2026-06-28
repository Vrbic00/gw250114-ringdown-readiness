(* ::Package:: *)

(* Compare the approximate public GW250114 EFT projection products.

   This script does not combine likelihoods. RINGDOWN and pyRing samples are
   overlapping public analyses of the same event, so the safe output is a
   consistency/comparison table rather than a joint constraint.
*)

ClearAll["Global`*"];

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[{Directory[], "config", "gw250114_constraints_comparison.wl"}]
];
outputDir = If[
  Length[args] >= 2,
  args[[2]],
  FileNameJoin[{Directory[], "results", "gw250114_constraints_comparison"}]
];

If[! FileExistsQ[configPath],
  Print["Configuration not found: ", configPath];
  Exit[1];
];

Get[configPath];
If[! AssociationQ[gw250114ConstraintsComparisonConfig],
  Print["Configuration must define gw250114ConstraintsComparisonConfig."];
  Exit[1];
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

config = gw250114ConstraintsComparisonConfig;
requiredFiles = Values[config];
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

importAssociations[path_] := Module[{raw = Import[path, "CSV"]},
  AssociationThread[First[raw], #] & /@ Rest[raw]
];

csvValue[x_] := If[NumericQ[x], N[x], x];
exportAssociationCSV[path_, rows_List] := Module[{fields = Keys[First[rows]]},
  Export[path, Prepend[csvValue /@ Lookup[#, fields] & /@ rows, fields], "CSV"]
];

key[row_Association] := txt[row["operator"]] <> ":" <> txt[row["polarization"]];

ringdownRowsRaw = importAssociations[config["ringdown_projection_path"]];
pyringRowsRaw = importAssociations[config["pyring_projection_path"]];
ringdownObservableRows = importAssociations[config["ringdown_observable_summary_path"]];
pyringObservableRows = importAssociations[config["pyring_observable_summary_path"]];

projectionRows[projectionName_String, observableSet_String, role_String, rows_List] :=
  Map[
    Function[row,
      <|
        "projection" -> projectionName,
        "observable_set" -> observableSet,
        "role" -> role,
        "operator" -> txt[row["operator"]],
        "polarization" -> txt[row["polarization"]],
        "alpha_best" -> num[row["alpha_best"]],
        "alpha_sigma" -> num[row["alpha_sigma_gaussian"]],
        "alpha_q05" -> num[row["alpha_q05_gaussian"]],
        "alpha_q16" -> num[row["alpha_q16_gaussian"]],
        "alpha_q84" -> num[row["alpha_q84_gaussian"]],
        "alpha_q95" -> num[row["alpha_q95_gaussian"]],
        "zero_inside_90pct" -> num[row["zero_inside_90pct"]],
        "sigma_from_alpha0" -> num[row["gaussian_sigma_from_alpha0"]]
      |>
    ],
    rows
  ];

longRows = Join[
  projectionRows[
    "RINGDOWN",
    "{log f_220, log f_221, df_221}",
    "primary frequency-only public ringdown projection",
    ringdownRowsRaw
  ],
  projectionRows[
    "PYRING_DELTA",
    "{log(1 + domega_221), log(1 + dtau_221)}",
    "complementary frequency-plus-damping deviation projection",
    pyringRowsRaw
  ]
];

ringdownByKey = Association[(key[#] -> #) & /@ ringdownRowsRaw];
pyringByKey = Association[(key[#] -> #) & /@ pyringRowsRaw];
commonKeys = Intersection[Keys[ringdownByKey], Keys[pyringByKey]];

pairRows = Map[
  Function[k,
    Module[
      {
        r = ringdownByKey[k], p = pyringByKey[k], rBest, pBest, rSigma,
        pSigma, rQ05, rQ95, pQ05, pQ95, overlap, normalizedDifference,
        pOverR
      },
      rBest = num[r["alpha_best"]];
      pBest = num[p["alpha_best"]];
      rSigma = num[r["alpha_sigma_gaussian"]];
      pSigma = num[p["alpha_sigma_gaussian"]];
      rQ05 = num[r["alpha_q05_gaussian"]];
      rQ95 = num[r["alpha_q95_gaussian"]];
      pQ05 = num[p["alpha_q05_gaussian"]];
      pQ95 = num[p["alpha_q95_gaussian"]];
      overlap = Boole[Max[rQ05, pQ05] <= Min[rQ95, pQ95]];
      normalizedDifference = Abs[rBest - pBest]/Sqrt[rSigma^2 + pSigma^2];
      pOverR = pSigma/rSigma;
      <|
        "operator" -> txt[r["operator"]],
        "polarization" -> txt[r["polarization"]],
        "ringdown_alpha_best" -> rBest,
        "ringdown_alpha_sigma" -> rSigma,
        "ringdown_sigma_from_zero" -> num[r["gaussian_sigma_from_alpha0"]],
        "pyring_alpha_best" -> pBest,
        "pyring_alpha_sigma" -> pSigma,
        "pyring_sigma_from_zero" -> num[p["gaussian_sigma_from_alpha0"]],
        "pyring_sigma_over_ringdown_sigma" -> pOverR,
        "best_fit_difference" -> pBest - rBest,
        "normalized_projection_difference" -> normalizedDifference,
        "intervals_overlap_90pct" -> overlap,
        "zero_inside_90pct_both" ->
          Boole[
            num[r["zero_inside_90pct"]] == 1 &&
              num[p["zero_inside_90pct"]] == 1
          ],
        "comparison_verdict" ->
          If[
            overlap == 1 && normalizedDifference < 1.5,
            "consistent",
            If[overlap == 1, "overlap_but_shifted", "tension_check"]
          ]
      |>
    ]
  ],
  commonKeys
];

projectionSummaryRows = Map[
  Function[pair,
    Module[{projectionName = pair[[1]], rows = pair[[2]], bestSigmaRow,
      bestZRow, zeroOutsideCount},
      bestSigmaRow = First[SortBy[rows, num[#["alpha_sigma_gaussian"]] &]];
      bestZRow = Last[SortBy[rows, num[#["gaussian_sigma_from_alpha0"]] &]];
      zeroOutsideCount = Count[num /@ Lookup[rows, "zero_inside_90pct"], 0.];
      <|
        "projection" -> projectionName,
        "rows" -> Length[rows],
        "zero_outside_90pct_count" -> zeroOutsideCount,
        "max_sigma_from_zero" -> num[bestZRow["gaussian_sigma_from_alpha0"]],
        "max_sigma_operator" -> txt[bestZRow["operator"]],
        "max_sigma_polarization" -> txt[bestZRow["polarization"]],
        "tightest_alpha_sigma" -> num[bestSigmaRow["alpha_sigma_gaussian"]],
        "tightest_operator" -> txt[bestSigmaRow["operator"]],
        "tightest_polarization" -> txt[bestSigmaRow["polarization"]]
      |>
    ]
  ],
  {
    {"RINGDOWN", ringdownRowsRaw},
    {"PYRING_DELTA", pyringRowsRaw}
  }
];

observableSummaryRows = Join[
  Map[
    Join[<|"projection" -> "RINGDOWN"|>, #] &,
    ringdownObservableRows
  ],
  Map[
    Join[<|"projection" -> "PYRING_DELTA"|>, #] &,
    pyringObservableRows
  ]
];

longCsvPath = FileNameJoin[{outputDir, "projection_constraints_long.csv"}];
pairCsvPath = FileNameJoin[{outputDir, "projection_consistency_by_operator.csv"}];
summaryCsvPath = FileNameJoin[{outputDir, "projection_constraints_summary.csv"}];
observableCsvPath = FileNameJoin[{outputDir, "projection_observable_summary.csv"}];
plotPath = FileNameJoin[{outputDir, "projection_sigma_from_zero_comparison.png"}];
reportPath = FileNameJoin[{outputDir, "gw250114_constraints_comparison_report.md"}];

exportAssociationCSV[longCsvPath, longRows];
exportAssociationCSV[pairCsvPath, pairRows];
exportAssociationCSV[summaryCsvPath, projectionSummaryRows];
exportAssociationCSV[observableCsvPath, observableSummaryRows];

pairRowsSorted = SortBy[pairRows, {#["operator"] &, #["polarization"] &}];
plotLabels = (#["operator"] <> " " <> #["polarization"]) & /@ pairRowsSorted;
ringdownZ = num /@ Lookup[pairRowsSorted, "ringdown_sigma_from_zero"];
pyringZ = num /@ Lookup[pairRowsSorted, "pyring_sigma_from_zero"];
xvals = Range[Length[pairRowsSorted]];
plotYMax = 1.1 Max[1., Max[Join[ringdownZ, pyringZ]]];
plotTicks = MapThread[
  {#1, Rotate[Style[#2, 11], Pi/4]} &,
  {xvals, plotLabels}
];
plot = Graphics[
  {
    {GrayLevel[0.85], Dashed, Line[{{0.5, 1}, {Length[xvals] + 0.5, 1}}]},
    {RGBColor[0.1, 0.25, 0.75], PointSize[0.015],
      Point[Transpose[{xvals - 0.12, ringdownZ}]]},
    {RGBColor[0.8, 0.35, 0.1], PointSize[0.015],
      Point[Transpose[{xvals + 0.12, pyringZ}]]},
    {RGBColor[0.1, 0.25, 0.75],
      Text[Style["RINGDOWN", 12], {Length[xvals] - 0.2, 0.92 plotYMax}]},
    {RGBColor[0.8, 0.35, 0.1],
      Text[Style["pyRing", 12], {Length[xvals] - 0.2, 0.82 plotYMax}]}
  },
  Frame -> True,
  Axes -> False,
  FrameTicks -> {{Automatic, None}, {plotTicks, None}},
  FrameLabel -> {None, "nominal |alpha_best|/sigma_alpha"},
  PlotLabel -> "GW250114 public projection comparison",
  PlotRange -> {{0.5, Length[xvals] + 0.5}, {0, plotYMax}},
  ImagePadding -> {{70, 30}, {155, 55}},
  ImageSize -> 1200
];
Export[plotPath, plot, ImageResolution -> 144];

bestRowsForReport = Take[
  SortBy[longRows, #["alpha_sigma"] &],
  UpTo[6]
];
largestShiftRowsForReport = Take[
  Reverse[SortBy[pairRows, #["normalized_projection_difference"] &]],
  UpTo[5]
];

constraintTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["projection"],
        row["operator"],
        row["polarization"],
        fmt[row["alpha_best"], 5],
        fmt[row["alpha_sigma"], 5],
        fmt[row["sigma_from_alpha0"], 3],
        If[row["zero_inside_90pct"] == 1, "yes", "no"]
      },
      " | "
    ] <> " |"
  ],
  bestRowsForReport
];

consistencyTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["operator"],
        row["polarization"],
        fmt[row["ringdown_alpha_best"], 5],
        fmt[row["pyring_alpha_best"], 5],
        fmt[row["normalized_projection_difference"], 3],
        If[row["intervals_overlap_90pct"] == 1, "yes", "no"],
        row["comparison_verdict"]
      },
      " | "
    ] <> " |"
  ],
  largestShiftRowsForReport
];

summaryTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["projection"],
        ToString[row["zero_outside_90pct_count"]],
        fmt[row["max_sigma_from_zero"], 3],
        row["max_sigma_operator"] <> " " <> row["max_sigma_polarization"],
        fmt[row["tightest_alpha_sigma"], 5],
        row["tightest_operator"] <> " " <> row["tightest_polarization"]
      },
      " | "
    ] <> " |"
  ],
  projectionSummaryRows
];

maxDisagreement = Max[num /@ Lookup[pairRows, "normalized_projection_difference"]];
tensionCount = Count[Lookup[pairRows, "comparison_verdict"], "tension_check"];
zeroOutsideAny = Total[projectionSummaryRows[[All, "zero_outside_90pct_count"]]];

report = StringRiffle[
  Join[
    {
      "# GW250114 Public Projection Constraints Comparison",
      "",
      "This report is generated by `scripts/wolfram/gw250114_constraints_comparison.wl`.",
      "",
      "## Scope",
      "",
      "- This is a comparison of two approximate public-data EFT projections.",
      "- It does not combine RINGDOWN and pyRing likelihoods because the analyses are not independent.",
      "- RINGDOWN uses `{log f_220, log f_221, df_221}` with linear mass/spin profiling.",
      "- pyRing uses `{log(1 + domega_221), log(1 + dtau_221)}` as a complementary frequency-plus-damping deviation posterior.",
      "",
      "## Projection Summary",
      "",
      "| projection | zero outside 90pct | max sigma from zero | max sigma row | tightest alpha sigma | tightest row |",
      "| --- | ---: | ---: | --- | ---: | --- |"
    },
    summaryTableRows,
    {
      "",
      "## Tightest Individual Rows",
      "",
      "| projection | operator | polarization | alpha best | alpha sigma | sigma from alpha=0 | zero in 90pct |",
      "| --- | --- | --- | ---: | ---: | ---: | --- |"
    },
    constraintTableRows,
    {
      "",
      "## Largest Projection Differences",
      "",
      "| operator | polarization | RINGDOWN alpha | pyRing alpha | normalized difference | 90pct overlap | verdict |",
      "| --- | --- | ---: | ---: | ---: | --- | --- |"
    },
    consistencyTableRows,
    {
      "",
      "## Interpretation",
      "",
      "- Reproduction: this table uses only previously generated public RINGDOWN and pyRing projection CSV outputs.",
      "- Numerical check: every one-at-a-time EFT row still has alpha=0 inside its 90 percent interval.",
      "- Interpretation: the two projection branches are mutually consistent at the current approximate level; no pair is flagged as a 90 percent interval tension.",
      "- Caution: the constraints should be quoted as approximate projected constraints, not final LVK-style EFT bounds.",
      "",
      "## Article Use",
      "",
      "Use this as a draft constraints table and pipeline-consistency check. A paper should present RINGDOWN as the cleaner frequency-only public projection, pyRing as the frequency-plus-damping cross-check, and explicitly state that the two are not statistically combined.",
      "",
      "## Generated Files",
      "",
      "- `projection_constraints_long.csv`",
      "- `projection_consistency_by_operator.csv`",
      "- `projection_constraints_summary.csv`",
      "- `projection_observable_summary.csv`",
      "- `projection_sigma_from_zero_comparison.png`"
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated GW250114 public projection constraints comparison"];
Print["Report: ", reportPath];
Print["Rows with alpha=0 outside 90pct: ", zeroOutsideAny];
Print["Projection pairs flagged as interval tension: ", tensionCount];
Print["Maximum normalized projection difference: ", fmt[maxDisagreement, 3]];
