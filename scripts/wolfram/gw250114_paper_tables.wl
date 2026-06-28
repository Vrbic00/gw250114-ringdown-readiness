(* ::Package:: *)

(* Build curated paper-facing tables from generated GW250114 project outputs. *)

ClearAll["Global`*"];

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[{Directory[], "config", "gw250114_paper_tables.wl"}]
];
outputDir = If[
  Length[args] >= 2,
  args[[2]],
  FileNameJoin[{Directory[], "results", "gw250114_paper_tables"}]
];

If[! FileExistsQ[configPath],
  Print["Configuration not found: ", configPath];
  Exit[1];
];

Get[configPath];
If[! AssociationQ[gw250114PaperTablesConfig],
  Print["Configuration must define gw250114PaperTablesConfig."];
  Exit[1];
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

config = gw250114PaperTablesConfig;
requiredFiles = Values[config];
missingFiles = Select[requiredFiles, ! FileExistsQ[#] &];
If[Length[missingFiles] > 0,
  Print["Missing input files: ", missingFiles];
  Exit[1];
];

num[value_] := N[If[NumericQ[value], value, ToExpression[ToString[value]]]];
txt[value_] := ToString[value];
csvValue[x_] := If[NumericQ[x], N[x], x];

importAssociations[path_] := Module[{raw = Import[path, "CSV"]},
  AssociationThread[First[raw], #] & /@ Rest[raw]
];

exportAssociationCSV[path_, rows_List] := Module[{fields = Keys[First[rows]]},
  Export[path, Prepend[csvValue /@ Lookup[#, fields] & /@ rows, fields], "CSV"]
];

fmt[Infinity, digits_Integer : 4] := "inf";
fmt[value_?NumericQ, digits_Integer : 4] := Module[
  {x = N[value]},
  ToString[NumberForm[x, {14, digits}, ExponentFunction -> (Null &)], OutputForm]
];
fmt[value_, digits_Integer : 4] := ToString[value];

constraintsRows = importAssociations[config["constraints_long_path"]];
consistencyRows = importAssociations[config["consistency_path"]];
filterRows = importAssociations[config["filter_summary_path"]];
linearizedRows = importAssociations[config["linearized_summary_path"]];
ringdownSummaryRows = importAssociations[config["public_ringdown_summary_path"]];

mainConstraints = Map[
  Function[row,
    <|
      "projection" -> txt[row["projection"]],
      "observable_set" -> txt[row["observable_set"]],
      "operator" -> txt[row["operator"]],
      "polarization" -> txt[row["polarization"]],
      "alpha_best" -> num[row["alpha_best"]],
      "alpha_sigma" -> num[row["alpha_sigma"]],
      "alpha_90_lower" -> num[row["alpha_q05"]],
      "alpha_90_upper" -> num[row["alpha_q95"]],
      "sigma_from_alpha0" -> num[row["sigma_from_alpha0"]],
      "zero_inside_90pct" -> num[row["zero_inside_90pct"]]
    |>
  ],
  constraintsRows
];

pipelineConsistency = Map[
  Function[row,
    <|
      "operator" -> txt[row["operator"]],
      "polarization" -> txt[row["polarization"]],
      "ringdown_alpha_best" -> num[row["ringdown_alpha_best"]],
      "pyring_alpha_best" -> num[row["pyring_alpha_best"]],
      "normalized_projection_difference" ->
        num[row["normalized_projection_difference"]],
      "intervals_overlap_90pct" -> num[row["intervals_overlap_90pct"]],
      "comparison_verdict" -> txt[row["comparison_verdict"]]
    |>
  ],
  consistencyRows
];

filterRobustness = Map[
  Function[row,
    <|
      "scenario" -> txt[row["scenario"]],
      "samples" -> num[row["samples"]],
      "zero_outside_90pct_count" -> num[row["zero_outside_90pct_count"]],
      "max_sigma_from_zero" -> num[row["max_sigma_from_zero"]],
      "max_sigma_row" ->
        txt[row["max_sigma_operator"]] <> " " <>
          txt[row["max_sigma_polarization"]],
      "tightest_alpha_sigma" -> num[row["tightest_alpha_sigma"]],
      "tightest_row" ->
        txt[row["tightest_operator"]] <> " " <>
          txt[row["tightest_polarization"]]
    |>
  ],
  filterRows
];

linearizedRobustness = Map[
  Function[row,
    <|
      "projection" -> txt[row["projection"]],
      "zero_outside_empirical_90pct_count" ->
        num[row["zero_outside_empirical_90pct_count"]],
      "max_nominal_abs_median_over_sd" ->
        num[row["max_nominal_abs_median_over_sd"]],
      "max_nominal_row" ->
        txt[row["max_nominal_operator"]] <> " " <>
          txt[row["max_nominal_polarization"]],
      "max_abs_interval_asymmetry" ->
        num[row["max_abs_interval_asymmetry"]],
      "max_asymmetry_row" ->
        txt[row["max_asymmetry_operator"]] <> " " <>
          txt[row["max_asymmetry_polarization"]]
    |>
  ],
  linearizedRows
];

observableSelectionQ[row_] := MemberQ[
  {
    "RINGDOWN_220_221:df_221",
    "RINGDOWN_220_221:f_220",
    "RINGDOWN_220_221:f_221",
    "RINGDOWN_220_221:f_221_inferred_kerr",
    "PYRING_220_221_filtered:df_221_log",
    "PYRING_220_221_filtered:dtau_221",
    "PYRING_220_221_filtered:dtau_221_log",
    "PYRING_220_221_filtered:f_220",
    "PYRING_220_221_filtered:f_221_corrected"
  },
  txt[row["dataset"]] <> ":" <> txt[row["parameter"]]
];

publicObservables = Map[
  Function[row,
    <|
      "dataset" -> txt[row["dataset"]],
      "parameter" -> txt[row["parameter"]],
      "samples" -> num[row["n"]],
      "median" -> num[row["median"]],
      "q16" -> num[row["q16"]],
      "q84" -> num[row["q84"]],
      "q05" -> num[row["q05"]],
      "q95" -> num[row["q95"]]
    |>
  ],
  Select[ringdownSummaryRows, observableSelectionQ]
];

mainConstraintsPath = FileNameJoin[
  {outputDir, "table1_main_projected_constraints.csv"}
];
pipelineConsistencyPath = FileNameJoin[
  {outputDir, "table2_pipeline_consistency.csv"}
];
filterPath = FileNameJoin[
  {outputDir, "table3_pyring_filter_robustness.csv"}
];
linearizedPath = FileNameJoin[
  {outputDir, "table4_linearized_posterior_check.csv"}
];
observablesPath = FileNameJoin[
  {outputDir, "table5_public_ringdown_observables.csv"}
];
reportPath = FileNameJoin[
  {outputDir, "gw250114_paper_tables_report.md"}
];

exportAssociationCSV[mainConstraintsPath, mainConstraints];
exportAssociationCSV[pipelineConsistencyPath, pipelineConsistency];
exportAssociationCSV[filterPath, filterRobustness];
exportAssociationCSV[linearizedPath, linearizedRobustness];
exportAssociationCSV[observablesPath, publicObservables];

report = StringRiffle[
  {
    "# GW250114 Paper Tables",
    "",
    "This report is generated by `scripts/wolfram/gw250114_paper_tables.wl`.",
    "",
    "## Generated Tables",
    "",
    "- `table1_main_projected_constraints.csv`: main projected one-at-a-time EFT constraints.",
    "- `table2_pipeline_consistency.csv`: RINGDOWN vs pyRing consistency checks.",
    "- `table3_pyring_filter_robustness.csv`: lower-tail filter robustness summary.",
    "- `table4_linearized_posterior_check.csv`: empirical linearized posterior sanity check.",
    "- `table5_public_ringdown_observables.csv`: selected public ringdown observables.",
    "",
    "## Main Numerical Checks",
    "",
    "- Main constraints rows: " <> ToString[Length[mainConstraints]],
    "- Pipeline consistency rows: " <> ToString[Length[pipelineConsistency]],
    "- Filter robustness scenarios: " <> ToString[Length[filterRobustness]],
    "- Linearized projection summaries: " <> ToString[Length[linearizedRobustness]],
    "- Selected public observable rows: " <> ToString[Length[publicObservables]],
    "",
    "## Publication Note",
    "",
    "These tables are intended for manuscript drafting and machine-readable supplement material. They should be cited as projected public-data constraints, not as a full strain-level EFT likelihood."
  },
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated GW250114 paper-facing tables"];
Print["Report: ", reportPath];
Print["Main constraints rows: ", Length[mainConstraints]];
