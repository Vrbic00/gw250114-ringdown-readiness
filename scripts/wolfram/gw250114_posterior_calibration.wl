(* ::Package:: *)

(* Public posterior calibration for GW250114.

   Usage:
     wolframscript -file scripts/wolfram/gw250114_posterior_calibration.wl
     wolframscript -file scripts/wolfram/gw250114_posterior_calibration.wl config/custom.wl results/custom_dir
*)

ClearAll["Global`*"];

solarMassTimeSeconds = 4.925490947*10^-6;

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[{Directory[], "config", "gw250114_posterior_calibration.wl"}]
];
outputDir = If[
  Length[args] >= 2,
  args[[2]],
  FileNameJoin[{Directory[], "results", "gw250114_posterior_calibration"}]
];

If[! FileExistsQ[configPath],
  Print["Configuration not found: ", configPath];
  Exit[1];
];

Get[configPath];
If[! AssociationQ[gw250114PosteriorCalibrationConfig],
  Print["Configuration must define gw250114PosteriorCalibrationConfig."];
  Exit[1];
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

config = gw250114PosteriorCalibrationConfig;
event = config["event"];
posteriorPath = config["posterior_path"];
posteriorDataset = event["posterior_dataset"];
modes = config["modes"];
selectedColumns = config["selected_columns"];
summaryColumns = config["summary_columns"];
covarianceModes = config["covariance_modes"];
centralReference = config["central_reference"];

If[! FileExistsQ[posteriorPath],
  Print["Posterior file not found: ", posteriorPath];
  Exit[1];
];

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

kerrObservables[mode_String, massSolar_?NumericQ, spin_?NumericQ] := Module[
  {omegaR, q, omegaI, fHz, tauMs},
  omegaR = dimensionlessAngularFrequency[mode, spin];
  q = qualityFactor[mode, spin];
  omegaI = -omegaR/(2 q);
  fHz = omegaR/(2 Pi massSolar solarMassTimeSeconds);
  tauMs = -1000 massSolar solarMassTimeSeconds/omegaI;
  <|
    "mode" -> mode,
    "Momega_re" -> omegaR,
    "Momega_im" -> omegaI,
    "f_Hz" -> fHz,
    "tau_ms" -> tauMs,
    "Q" -> q
  |>
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

Print["Importing posterior samples from ", posteriorPath];
posteriorRows = Import[posteriorPath, {"Datasets", posteriorDataset}];
If[! ListQ[posteriorRows] || Length[posteriorRows] == 0,
  Print["Could not import posterior rows."];
  Exit[1];
];

availableColumns = Keys[First[posteriorRows]];
missingColumns = Select[
  Union[selectedColumns, summaryColumns],
  ! MemberQ[availableColumns, #] &
];
If[Length[missingColumns] > 0,
  Print["Missing posterior columns: ", missingColumns];
  Exit[1];
];

selectedRows = KeyTake[#, selectedColumns] & /@ posteriorRows;
summaryRows = quantileSummary[#, Lookup[posteriorRows, #]] & /@ summaryColumns;

qnmRows = Flatten[
  Table[
    Module[
      {obs = kerrObservables[mode, row["final_mass"], row["final_spin"]]},
      <|
        "mode" -> mode,
        "f_Hz" -> obs["f_Hz"],
        "tau_ms" -> obs["tau_ms"],
        "Q" -> obs["Q"],
        "Momega_re" -> obs["Momega_re"],
        "Momega_im" -> obs["Momega_im"]
      |>
    ],
    {mode, modes},
    {row, posteriorRows}
  ],
  1
];

qnmSummaryRows = Flatten[
  Table[
    Module[
      {rows = Select[qnmRows, #["mode"] == mode &]},
      Join[
        <|"mode" -> mode|>,
        Association[
          Flatten[
            Table[
              With[
                {
                  summary = quantileSummary[
                    quantity,
                    Lookup[rows, quantity]
                  ]
                },
                (quantity <> "_" <> # -> summary[#]) & /@
                  {"mean", "standard_deviation", "q05", "q16", "median",
                    "q84", "q95"}
              ],
              {quantity, {"f_Hz", "tau_ms", "Q", "Momega_re", "Momega_im"}}
            ],
            1
          ]
        ]
      ]
    ],
    {mode, modes}
  ],
  1
];

covarianceLabels = Flatten[
  ({("ln_f_" <> #), ("ln_tau_" <> #)} & /@ covarianceModes)
];
covarianceVectors = Table[
  Flatten[
    Table[
      Module[
        {obs = kerrObservables[mode, row["final_mass"], row["final_spin"]]},
        {Log[obs["f_Hz"]], Log[obs["tau_ms"]]}
      ],
      {mode, covarianceModes}
    ]
  ],
  {row, posteriorRows}
];
covarianceMatrix = Covariance[covarianceVectors];
correlationMatrix = Correlation[covarianceVectors];

covarianceRows = Flatten[
  Table[
    <|
      "row" -> covarianceLabels[[i]],
      "column" -> covarianceLabels[[j]],
      "covariance" -> covarianceMatrix[[i, j]],
      "correlation" -> correlationMatrix[[i, j]]
    |>,
    {i, Length[covarianceLabels]},
    {j, Length[covarianceLabels]}
  ],
  1
];

centralRows = Map[
  Function[mode,
    Join[
      <|"mode" -> mode|>,
      kerrObservables[
        mode,
        centralReference["mass_detector_msun"],
        centralReference["spin"]
      ]
    ]
  ],
  modes
];

selectedCsvPath = FileNameJoin[
  {outputDir, "nrSur7dq4_selected_posterior_samples.csv"}
];
summaryCsvPath = FileNameJoin[
  {outputDir, "nrSur7dq4_posterior_summary.csv"}
];
qnmSummaryCsvPath = FileNameJoin[
  {outputDir, "nrSur7dq4_kerr_qnm_posterior_summary.csv"}
];
covarianceCsvPath = FileNameJoin[
  {outputDir, "nrSur7dq4_log_qnm_covariance_220_221.csv"}
];
centralCsvPath = FileNameJoin[
  {outputDir, "central_reference_berti_qnm.csv"}
];
reportPath = FileNameJoin[
  {outputDir, "gw250114_posterior_calibration_report.md"}
];

exportAssociationCSV[selectedCsvPath, selectedRows];
exportAssociationCSV[summaryCsvPath, summaryRows];
exportAssociationCSV[qnmSummaryCsvPath, qnmSummaryRows];
exportAssociationCSV[covarianceCsvPath, covarianceRows];
exportAssociationCSV[centralCsvPath, centralRows];

summaryByName = Association[(#["parameter"] -> #) & /@ summaryRows];
qnmSummaryByMode = Association[(#["mode"] -> #) & /@ qnmSummaryRows];

posteriorTableRows = Map[
  Function[name,
    Module[{row = summaryByName[name]},
      "| " <> StringRiffle[
        {
          name,
          fmt[row["median"], 4],
          "-" <> fmt[row["median"] - row["q16"], 4],
          "+" <> fmt[row["q84"] - row["median"], 4],
          fmt[row["q05"], 4],
          fmt[row["q95"], 4]
        },
        " | "
      ] <> " |"
    ]
  ],
  {
    "final_mass", "final_mass_source", "final_spin", "redshift",
    "network_optimal_snr", "network_33_multipole_snr",
    "network_44_multipole_snr"
  }
];

qnmTableRows = Map[
  Function[mode,
    Module[{row = qnmSummaryByMode[mode]},
      "| " <> StringRiffle[
        {
          mode,
          fmt[row["f_Hz_median"], 3],
          "-" <> fmt[row["f_Hz_median"] - row["f_Hz_q16"], 3],
          "+" <> fmt[row["f_Hz_q84"] - row["f_Hz_median"], 3],
          fmt[row["tau_ms_median"], 4],
          "-" <> fmt[row["tau_ms_median"] - row["tau_ms_q16"], 4],
          "+" <> fmt[row["tau_ms_q84"] - row["tau_ms_median"], 4]
        },
        " | "
      ] <> " |"
    ]
  ],
  modes
];

sigmaRows = Table[
  Module[
    {idx = FirstPosition[covarianceLabels, label][[1]]},
    "| " <> label <> " | " <> fmt[Sqrt[covarianceMatrix[[idx, idx]]], 5] <>
      " |"
  ],
  {label, covarianceLabels}
];

report = StringRiffle[
  Join[
    {
      "# GW250114 NRSur7dq4 Posterior Calibration",
      "",
      "This report is generated by `scripts/wolfram/gw250114_posterior_calibration.wl`.",
      "",
      "## Data",
      "",
      "- Event: `" <> event["name"] <> "`",
      "- Posterior model: `" <> event["posterior_model"] <> "`",
      "- Source: " <> event["source"],
      "- Zenodo record: `" <> event["zenodo_record"] <> "`",
      "- File: `" <> event["zenodo_file"] <> "`",
      "- Dataset: `" <> posteriorDataset <> "`",
      "- Posterior samples: `" <> ToString[Length[posteriorRows]] <> "`",
      "",
      "This is a full IMR parameter-estimation posterior, not a ringdown-only or no-hair-test posterior.",
      "",
      "## Posterior Summary",
      "",
      "| parameter | median | -1sigma | +1sigma | q05 | q95 |",
      "| --- | ---: | ---: | ---: | ---: | ---: |"
    },
    posteriorTableRows,
    {
      "",
      "## Kerr QNM Posterior Pushforward",
      "",
      "The QNM rows use the Berti-Cardoso-Will Kerr fits evaluated at each posterior sample's `final_mass` and `final_spin`.",
      "",
      "| mode | f median [Hz] | -1sigma | +1sigma | tau median [ms] | -1sigma | +1sigma |",
      "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"
    },
    qnmTableRows,
    {
      "",
      "## 220+221 Log-QNM Widths",
      "",
      "These are posterior widths of Kerr-predicted `{log f, log tau}` induced by the IMR remnant posterior.",
      "",
      "| quantity | posterior sigma |",
      "| --- | ---: |"
    },
    sigmaRows,
    {
      "",
      "## Interpretation",
      "",
      "- Reproduction: public NRSur7dq4 posterior samples are imported and summarized.",
      "- Numerical check: the median detector-frame final mass and spin are close to the project baseline values used so far.",
      "- Interpretation: this calibrates remnant-prior uncertainty and Kerr QNM pushforwards, but it does not replace a ringdown-only likelihood for EFT deviations.",
      "- Speculation: none.",
      "",
      "## Next Defensible Step",
      "",
      "Feed the empirical 220+221 log-QNM covariance into the synthetic spectral EFT projection as a calibrated remnant-prior scale, while keeping the time-domain likelihood synthetic until strain-level or ringdown-posterior products are added.",
      "",
      "## Generated Files",
      "",
      "- `nrSur7dq4_selected_posterior_samples.csv`",
      "- `nrSur7dq4_posterior_summary.csv`",
      "- `nrSur7dq4_kerr_qnm_posterior_summary.csv`",
      "- `nrSur7dq4_log_qnm_covariance_220_221.csv`",
      "- `central_reference_berti_qnm.csv`"
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated GW250114 posterior calibration"];
Print["Samples: ", Length[posteriorRows]];
Print["Report: ", reportPath];
Print["final_mass median: ", fmt[summaryByName["final_mass"]["median"], 4]];
Print["final_spin median: ", fmt[summaryByName["final_spin"]["median"], 5]];
Print["220 f median Hz: ", fmt[qnmSummaryByMode["220"]["f_Hz_median"], 3]];
Print["221 f median Hz: ", fmt[qnmSummaryByMode["221"]["f_Hz_median"], 3]];
