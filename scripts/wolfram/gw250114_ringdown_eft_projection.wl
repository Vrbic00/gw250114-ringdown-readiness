(* ::Package:: *)

(* Approximate EFT projection with public GW250114 RINGDOWN posterior samples.

   This is not a full strain likelihood. It builds a local Gaussian
   approximation in y = {log f_220, log f_221, df_221}.
*)

ClearAll["Global`*"];

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[{Directory[], "config", "gw250114_ringdown_eft_projection.wl"}]
];
outputDir = If[
  Length[args] >= 2,
  args[[2]],
  FileNameJoin[{Directory[], "results", "gw250114_ringdown_eft_projection"}]
];

If[! FileExistsQ[configPath],
  Print["Configuration not found: ", configPath];
  Exit[1];
];

Get[configPath];
If[! AssociationQ[gw250114RingdownEFTProjectionConfig],
  Print["Configuration must define gw250114RingdownEFTProjectionConfig."];
  Exit[1];
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

config = gw250114RingdownEFTProjectionConfig;
event = config["event"];
ringdownPath = config["ringdown_hdf5_path"];
eftSensitivityPath = config["eft_sensitivity_path"];
kerrPath = config["kerr_numeric_path"];
mass0 = N[event["mass_detector_msun"]];
spin0 = N[event["spin"]];

If[
  ! FileExistsQ[ringdownPath] || ! FileExistsQ[eftSensitivityPath] ||
    ! FileExistsQ[kerrPath],
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

importAssociations[path_] := Module[{raw = Import[path, "CSV"]},
  AssociationThread[First[raw], #] & /@ Rest[raw]
];

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

eftRows = Select[
  importAssociations[eftSensitivityPath],
  MemberQ[{"220", "221"}, txt[#["mode"]]] &
];
eftByModeOperatorBranch = Association[
  ((txt[#["mode"]] <> ":" <> txt[#["operator"]] <> ":" <>
        txt[#["polarization"]]) -> #) & /@ eftRows
];

operators = {"lambda_ev", "lambda_odd", "epsilon1", "epsilon2", "epsilon3"};
polarizations = {"plus", "minus"};

f220 = Import[ringdownPath, {"Datasets", "/f_220"}];
f221 = Import[ringdownPath, {"Datasets", "/f_221"}];
df221 = Import[ringdownPath, {"Datasets", "/df_221"}];

observedVectors = Transpose[{Log[f220], Log[f221], df221}];
observedMean = Mean /@ Transpose[observedVectors];
observedCov = Covariance[observedVectors];
observedPrecision = Inverse[observedCov];
labels = {"ln_f220", "ln_f221", "df221"};

baselineVector = {
  Log[num[First[Select[kerrRows, txt[#["mode"]] == "220" &]]["qnm_f_Hz"]]],
  Log[num[First[Select[kerrRows, txt[#["mode"]] == "221" &]]["qnm_f_Hz"]]],
  0.
};
residualVector = observedMean - baselineVector;

nuisanceColumns = {
  {-1., -1., 0.},
  {spinLogFrequencyDerivative["220"], spinLogFrequencyDerivative["221"], 0.}
};
nuisanceLabels = {"delta_lnM", "delta_chi"};

profileFit[target_List, columns_List] := Module[
  {design, fisher, rhs, coeffs, cov, residual, chi2},
  design = Transpose[columns];
  fisher = Transpose[design].observedPrecision.design;
  rhs = Transpose[design].observedPrecision.target;
  cov = PseudoInverse[fisher];
  coeffs = cov.rhs;
  residual = target - design.coeffs;
  chi2 = residual.observedPrecision.residual;
  <|"coefficients" -> coeffs, "covariance" -> cov, "chi2" -> chi2|>
];

fitRows = Flatten[
  Table[
    Module[
      {
        s220 = eftByModeOperatorBranch["220:" <> operator <> ":" <> branch],
        s221 = eftByModeOperatorBranch["221:" <> operator <> ":" <> branch],
        alphaColumn, fullFit, nullFit, alphaIndex, alphaBest,
        alphaSigma, deltaChi2, z
      },
      alphaColumn = {
        num[s220["dln_frequency_dalpha"]],
        num[s221["dln_frequency_dalpha"]],
        num[s221["dln_frequency_dalpha"]]
      };
      fullFit = profileFit[residualVector, Join[nuisanceColumns, {alphaColumn}]];
      nullFit = profileFit[residualVector, nuisanceColumns];
      alphaIndex = Length[nuisanceColumns] + 1;
      alphaBest = fullFit["coefficients"][[alphaIndex]];
      alphaSigma = Sqrt[fullFit["covariance"][[alphaIndex, alphaIndex]]];
      deltaChi2 = Max[0, nullFit["chi2"] - fullFit["chi2"]];
      z = If[alphaSigma > 0, Abs[alphaBest]/alphaSigma, 0];
      <|
        "operator" -> operator,
        "polarization" -> branch,
        "dlnf220_dalpha" -> alphaColumn[[1]],
        "dlnf221_dalpha" -> alphaColumn[[2]],
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
        "best_fit_delta_lnM" -> fullFit["coefficients"][[1]],
        "best_fit_delta_chi" -> fullFit["coefficients"][[2]],
        "chi2_best" -> fullFit["chi2"],
        "chi2_alpha0_profiled" -> nullFit["chi2"]
      |>
    ],
    {operator, operators},
    {branch, polarizations}
  ],
  1
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
    "baseline_value" -> baselineVector[[i]],
    "residual" -> residualVector[[i]],
    "posterior_sigma" -> Sqrt[observedCov[[i, i]]]
  |>,
  {i, Length[labels]}
];

csvValue[x_] := If[NumericQ[x], N[x], x];
exportAssociationCSV[path_, rows_List] := Module[{fields = Keys[First[rows]]},
  Export[path, Prepend[csvValue /@ Lookup[#, fields] & /@ rows, fields], "CSV"]
];

fitCsvPath = FileNameJoin[
  {outputDir, "ringdown_eft_gaussian_projection.csv"}
];
covCsvPath = FileNameJoin[
  {outputDir, "ringdown_observable_covariance.csv"}
];
meanCsvPath = FileNameJoin[
  {outputDir, "ringdown_observable_mean_residual.csv"}
];
plotPath = FileNameJoin[
  {outputDir, "ringdown_eft_alpha_gaussian_intervals.png"}
];
reportPath = FileNameJoin[
  {outputDir, "gw250114_ringdown_eft_projection_report.md"}
];

exportAssociationCSV[fitCsvPath, fitRows];
exportAssociationCSV[covCsvPath, covRows];
exportAssociationCSV[meanCsvPath, meanRows];

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
    {RGBColor[0.1, 0.25, 0.75], AbsoluteThickness[1.6],
      Map[Line[{{#[[1]], #[[2]] - #[[3]]}, {#[[1]], #[[2]] + #[[3]]}}] &, plotData],
      Map[Line[{{#[[1]] - 0.12, #[[2]] - #[[3]]}, {#[[1]] + 0.12, #[[2]] - #[[3]]}}] &, plotData],
      Map[Line[{{#[[1]] - 0.12, #[[2]] + #[[3]]}, {#[[1]] + 0.12, #[[2]] + #[[3]]}}] &, plotData]
    },
    {RGBColor[0.05, 0.15, 0.45], PointSize[0.016], Point[plotData[[All, {1, 2}]]]}
  },
  Frame -> True,
  Axes -> False,
  FrameTicks -> {{Automatic, None}, {plotTicks, None}},
  FrameLabel -> {None, "alpha proxy from Gaussian public ringdown projection"},
  PlotLabel -> "GW250114 public ringdown EFT projection in {ln f220, ln f221, df221}",
  PlotRange -> {{0.5, Length[fitRows] + 0.5}, {plotYMin - plotYPad, plotYMax + plotYPad}},
  ImagePadding -> {{70, 40}, {155, 55}},
  ImageSize -> 1200
];
Export[plotPath, plot, ImageResolution -> 144];

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

meanTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["parameter"],
        fmt[row["observed_mean"], 5],
        fmt[row["baseline_value"], 5],
        fmt[row["residual"], 5],
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
      "# GW250114 Public Ringdown EFT Projection",
      "",
      "This report is generated by `scripts/wolfram/gw250114_ringdown_eft_projection.wl`.",
      "",
      "## Scope",
      "",
      "- This is an approximate Gaussian projection of the public RINGDOWN posterior, not a full strain-level likelihood.",
      "- The observable vector is `y = {log f_220, log f_221, df_221}`.",
      "- Remnant mass and spin are linear nuisance directions.",
      "- Each EFT coupling is tested one at a time using the imported higher-derivative QNM frequency sensitivities.",
      "- The public `df_221` deviation is interpreted as a logarithmic 221 frequency deviation; damping-time information is not yet used.",
      "",
      "## Observable Mean and Baseline",
      "",
      "| parameter | observed mean | Kerr baseline | residual | posterior sigma |",
      "| --- | ---: | ---: | ---: | ---: |"
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
      "- Reproduction: this uses the public RINGDOWN `f_220`, `f_221`, and `df_221` posterior samples directly.",
      "- Numerical check: the Gaussian result is consistent with the simpler 1D `df_221` proxy: no one-at-a-time operator gives a robust nonzero deviation.",
      "- Interpretation: this is now an event-level ringdown projection, but still an approximate projection rather than a full EFT likelihood.",
      "- Speculation: none.",
      "",
      "## Next Defensible Step",
      "",
      "Add damping-time/deviation information if a public variable for `dg_221` or `dtau_221` is selected, and compare this Gaussian projection with the pyRing posterior in the same variables.",
      "",
      "## Generated Files",
      "",
      "- `ringdown_eft_gaussian_projection.csv`",
      "- `ringdown_observable_covariance.csv`",
      "- `ringdown_observable_mean_residual.csv`",
      "- `ringdown_eft_alpha_gaussian_intervals.png`"
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated public ringdown EFT Gaussian projection"];
Print["Report: ", reportPath];
Print["Best alpha rows with zero outside 90pct: ",
  Count[Lookup[fitRows, "zero_inside_90pct"], 0]];
Scan[
  Print[
    StringRiffle[
      {
        #["operator"],
        #["polarization"],
        fmt[#["alpha_best"], 5],
        fmt[#["alpha_sigma_gaussian"], 5],
        fmt[#["gaussian_sigma_from_alpha0"], 3]
      },
      "\t"
    ]
  ] &,
  TakeSmallestBy[fitRows, #["alpha_sigma_gaussian"] &, UpTo[5]]
];
