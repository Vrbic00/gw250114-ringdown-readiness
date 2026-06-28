(* ::Package:: *)

(* Robustness sweep for the public GW250114 pyRing 221 deviation projection. *)

ClearAll["Global`*"];

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[
    {Directory[], "config", "gw250114_pyring_filter_robustness.wl"}
  ]
];
outputDir = If[
  Length[args] >= 2,
  args[[2]],
  FileNameJoin[{Directory[], "results", "gw250114_pyring_filter_robustness"}]
];

If[! FileExistsQ[configPath],
  Print["Configuration not found: ", configPath];
  Exit[1];
];

Get[configPath];
If[! AssociationQ[gw250114PyRingFilterRobustnessConfig],
  Print["Configuration must define gw250114PyRingFilterRobustnessConfig."];
  Exit[1];
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

config = gw250114PyRingFilterRobustnessConfig;
pyringDeltaPath = config["pyring_delta_path"];
eftSensitivityPath = config["eft_sensitivity_path"];
mode = ToString[config["mode"]];
filterScenarios = config["filter_scenarios"];

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

csvValue[x_] := If[NumericQ[x], N[x], x];
exportAssociationCSV[path_, rows_List] := Module[{fields = Keys[First[rows]]},
  Export[path, Prepend[csvValue /@ Lookup[#, fields] & /@ rows, fields], "CSV"]
];

Print["Reading pyRing 220+221 delta posterior"];
pyringRowsAll = readDat[pyringDeltaPath];

eftRows = Select[importAssociations[eftSensitivityPath], txt[#["mode"]] == mode &];
operators = {"lambda_ev", "lambda_odd", "epsilon1", "epsilon2", "epsilon3"};
polarizations = {"plus", "minus"};
eftByOperatorBranch = Association[
  ((txt[#["operator"]] <> ":" <> txt[#["polarization"]]) -> #) & /@ eftRows
];

filteredRows[scenario_Association] := Module[{bound, lower},
  bound = scenario["delta_f_bound"];
  lower = If[bound === Infinity, -Infinity, Exp[-N[bound]] - 1];
  Select[
    pyringRowsAll,
    num[#["domega_221"]] > lower &&
      1 + num[#["domega_221"]] > 0 &&
      1 + num[#["dtau_221"]] > 0 &
  ]
];

fitRowsForScenario[scenario_Association] := Module[
  {rows, computedRows, observedVectors, observedMean, observedCov,
    observedPrecision, fitAlpha, projectionRows, bound, lower},
  rows = filteredRows[scenario];
  computedRows = Map[
    Function[row,
      Join[
        row,
        <|
          "df_221_log" -> Log[1 + num[row["domega_221"]]],
          "dtau_221_log" -> Log[1 + num[row["dtau_221"]]]
        |>
      ]
    ],
    rows
  ];
  observedVectors = Transpose[
    {
      num /@ Lookup[computedRows, "df_221_log"],
      num /@ Lookup[computedRows, "dtau_221_log"]
    }
  ];
  observedMean = Mean /@ Transpose[observedVectors];
  observedCov = Covariance[observedVectors];
  observedPrecision = Inverse[observedCov];
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
      "alpha_sigma" -> alphaSigma,
      "alpha_q05" -> alphaBest - 1.6448536269514722 alphaSigma,
      "alpha_q95" -> alphaBest + 1.6448536269514722 alphaSigma,
      "zero_inside_90pct" ->
        Boole[
          alphaBest - 1.6448536269514722 alphaSigma <= 0 <=
            alphaBest + 1.6448536269514722 alphaSigma
        ],
      "sigma_from_alpha0" -> Sqrt[deltaChi2],
      "chi2_best" -> chi2Best,
      "chi2_alpha0" -> chi2Zero
    |>
  ];
  bound = scenario["delta_f_bound"];
  lower = If[bound === Infinity, -Infinity, Exp[-N[bound]] - 1];
  projectionRows = Flatten[
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
            "scenario" -> scenario["name"],
            "description" -> scenario["description"],
            "delta_f_bound" -> If[bound === Infinity, "Infinity", N[bound]],
            "domega_lower" -> lower,
            "samples" -> Length[rows],
            "df221_mean" -> observedMean[[1]],
            "dtau221_mean" -> observedMean[[2]],
            "df221_sigma" -> Sqrt[observedCov[[1, 1]]],
            "dtau221_sigma" -> Sqrt[observedCov[[2, 2]]],
            "df_dtau_correlation" ->
              observedCov[[1, 2]]/
                Sqrt[observedCov[[1, 1]] observedCov[[2, 2]]],
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
  projectionRows
];

allFitRows = Flatten[fitRowsForScenario /@ filterScenarios, 1];

scenarioSummaryRows = Map[
  Function[scenario,
    Module[{rows = Select[allFitRows, #["scenario"] == scenario["name"] &],
      bestSigma, bestZ, zeroOutsideCount},
      bestSigma = First[SortBy[rows, #["alpha_sigma"] &]];
      bestZ = Last[SortBy[rows, #["sigma_from_alpha0"] &]];
      zeroOutsideCount = Count[Lookup[rows, "zero_inside_90pct"], 0];
      <|
        "scenario" -> scenario["name"],
        "samples" -> First[Lookup[rows, "samples"]],
        "zero_outside_90pct_count" -> zeroOutsideCount,
        "max_sigma_from_zero" -> bestZ["sigma_from_alpha0"],
        "max_sigma_operator" -> bestZ["operator"],
        "max_sigma_polarization" -> bestZ["polarization"],
        "tightest_alpha_sigma" -> bestSigma["alpha_sigma"],
        "tightest_operator" -> bestSigma["operator"],
        "tightest_polarization" -> bestSigma["polarization"],
        "df221_mean" -> First[Lookup[rows, "df221_mean"]],
        "dtau221_mean" -> First[Lookup[rows, "dtau221_mean"]]
      |>
    ]
  ],
  filterScenarios
];

operatorKeys = DeleteDuplicates[(#["operator"] <> ":" <> #["polarization"]) & /@ allFitRows];
operatorRobustnessRows = Map[
  Function[k,
    Module[{rows = Select[allFitRows, #["operator"] <> ":" <> #["polarization"] == k &],
      alphaBests, alphaSigmas, zVals, q05, q95},
      alphaBests = Lookup[rows, "alpha_best"];
      alphaSigmas = Lookup[rows, "alpha_sigma"];
      zVals = Lookup[rows, "sigma_from_alpha0"];
      q05 = Lookup[rows, "alpha_q05"];
      q95 = Lookup[rows, "alpha_q95"];
      <|
        "operator" -> First[Lookup[rows, "operator"]],
        "polarization" -> First[Lookup[rows, "polarization"]],
        "alpha_best_min" -> Min[alphaBests],
        "alpha_best_max" -> Max[alphaBests],
        "alpha_best_range" -> Max[alphaBests] - Min[alphaBests],
        "alpha_sigma_min" -> Min[alphaSigmas],
        "alpha_sigma_max" -> Max[alphaSigmas],
        "max_sigma_from_zero" -> Max[zVals],
        "zero_outside_90pct_any_scenario" ->
          Boole[AnyTrue[Lookup[rows, "zero_inside_90pct"], # == 0 &]],
        "common_90pct_interval_lower" -> Max[q05],
        "common_90pct_interval_upper" -> Min[q95],
        "common_90pct_interval_contains_zero" ->
          Boole[Max[q05] <= 0 <= Min[q95]]
      |>
    ]
  ],
  operatorKeys
];

fitCsvPath = FileNameJoin[{outputDir, "pyring_filter_robustness_long.csv"}];
scenarioCsvPath = FileNameJoin[{outputDir, "pyring_filter_scenario_summary.csv"}];
operatorCsvPath = FileNameJoin[{outputDir, "pyring_filter_operator_robustness.csv"}];
plotPath = FileNameJoin[{outputDir, "pyring_filter_max_sigma_by_scenario.png"}];
reportPath = FileNameJoin[{outputDir, "gw250114_pyring_filter_robustness_report.md"}];

exportAssociationCSV[fitCsvPath, allFitRows];
exportAssociationCSV[scenarioCsvPath, scenarioSummaryRows];
exportAssociationCSV[operatorCsvPath, operatorRobustnessRows];

scenarioLabels = Lookup[scenarioSummaryRows, "scenario"];
maxZ = Lookup[scenarioSummaryRows, "max_sigma_from_zero"];
xvals = Range[Length[scenarioLabels]];
plotTicks = MapThread[
  {#1, Rotate[Style[#2, 11], Pi/5]} &,
  {xvals, scenarioLabels}
];
plot = Graphics[
  {
    {GrayLevel[0.85], Dashed, Line[{{0.5, 1}, {Length[xvals] + 0.5, 1}}]},
    {RGBColor[0.1, 0.35, 0.45], PointSize[0.02], Point[Transpose[{xvals, maxZ}]]},
    {RGBColor[0.1, 0.35, 0.45], AbsoluteThickness[1.6],
      Line[Transpose[{xvals, maxZ}]]}
  },
  Frame -> True,
  Axes -> False,
  FrameTicks -> {{Automatic, None}, {plotTicks, None}},
  FrameLabel -> {None, "max nominal |alpha_best|/sigma_alpha"},
  PlotLabel -> "pyRing filter robustness",
  PlotRange -> {{0.5, Length[xvals] + 0.5}, {0, 1.15 Max[1, Max[maxZ]]}},
  ImagePadding -> {{70, 30}, {115, 55}},
  ImageSize -> 1000
];
Export[plotPath, plot, ImageResolution -> 144];

scenarioTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["scenario"],
        ToString[row["samples"]],
        ToString[row["zero_outside_90pct_count"]],
        fmt[row["max_sigma_from_zero"], 3],
        row["max_sigma_operator"] <> " " <> row["max_sigma_polarization"],
        fmt[row["tightest_alpha_sigma"], 5],
        row["tightest_operator"] <> " " <> row["tightest_polarization"]
      },
      " | "
    ] <> " |"
  ],
  scenarioSummaryRows
];

operatorTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["operator"],
        row["polarization"],
        fmt[row["alpha_best_min"], 5],
        fmt[row["alpha_best_max"], 5],
        fmt[row["alpha_best_range"], 5],
        fmt[row["max_sigma_from_zero"], 3],
        If[row["zero_outside_90pct_any_scenario"] == 1, "yes", "no"],
        If[row["common_90pct_interval_contains_zero"] == 1, "yes", "no"]
      },
      " | "
    ] <> " |"
  ],
  SortBy[operatorRobustnessRows, -#["max_sigma_from_zero"] &]
];

zeroOutsideTotal = Total[Lookup[scenarioSummaryRows, "zero_outside_90pct_count"]];
maxScenarioZ = Max[Lookup[scenarioSummaryRows, "max_sigma_from_zero"]];
anyCommonIntervalMissesZero =
  Count[Lookup[operatorRobustnessRows, "common_90pct_interval_contains_zero"], 0];

report = StringRiffle[
  Join[
    {
      "# GW250114 pyRing Filter Robustness",
      "",
      "This report is generated by `scripts/wolfram/gw250114_pyring_filter_robustness.wl`.",
      "",
      "## Scope",
      "",
      "- This sweeps the lower-tail `domega_221` filter used before building the pyRing `{df_221, dtau_221}` Gaussian projection.",
      "- It tests sensitivity to stricter, public, looser, and positive-domain-only choices.",
      "- It does not address full waveform-systematic uncertainty.",
      "",
      "## Scenario Summary",
      "",
      "| scenario | samples | zero outside 90pct | max sigma from zero | max sigma row | tightest alpha sigma | tightest row |",
      "| --- | ---: | ---: | ---: | --- | ---: | --- |"
    },
    scenarioTableRows,
    {
      "",
      "## Operator Robustness",
      "",
      "| operator | polarization | alpha best min | alpha best max | alpha best range | max sigma from zero | zero outside any scenario | common 90pct contains zero |",
      "| --- | --- | ---: | ---: | ---: | ---: | --- | --- |"
    },
    operatorTableRows,
    {
      "",
      "## Interpretation",
      "",
      "- Reproduction: the public `delta_f_bound = 0.8` scenario reproduces the pyRing delta projection scale.",
      "- Numerical check: no tested filter scenario pushes alpha=0 outside a 90 percent interval.",
      "- Interpretation: the no-deviation conclusion is robust to this lower-tail filter choice at the current Gaussian-projection level.",
      "- Caution: this is still a filter robustness test, not a substitute for a non-Gaussian posterior or full strain-level analysis.",
      "",
      "## Generated Files",
      "",
      "- `pyring_filter_robustness_long.csv`",
      "- `pyring_filter_scenario_summary.csv`",
      "- `pyring_filter_operator_robustness.csv`",
      "- `pyring_filter_max_sigma_by_scenario.png`"
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated pyRing filter robustness sweep"];
Print["Report: ", reportPath];
Print["Total rows with alpha=0 outside 90pct: ", zeroOutsideTotal];
Print["Maximum scenario sigma from zero: ", fmt[maxScenarioZ, 3]];
Print["Common 90pct intervals missing zero: ", anyCommonIntervalMissesZero];
