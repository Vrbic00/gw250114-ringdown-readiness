(* ::Package:: *)

(* Toy time-domain likelihood for synthetic EFT ringdown fingerprints.

   Usage:
     wolframscript -file scripts/wolfram/ringdown_toy_time_domain_likelihood.wl
     wolframscript -file scripts/wolfram/ringdown_toy_time_domain_likelihood.wl config/custom.wl results/custom_dir
*)

ClearAll["Global`*"];

solarMassTimeSeconds = 4.925490947*10^-6;

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[
    {Directory[], "config", "ringdown_toy_time_domain_likelihood.wl"}
  ]
];
outputDir = If[
  Length[args] >= 2,
  args[[2]],
  FileNameJoin[{Directory[], "results", "ringdown_toy_time_domain_likelihood"}]
];

If[! FileExistsQ[configPath],
  Print["Configuration not found: ", configPath];
  Exit[1];
];

Get[configPath];
If[! AssociationQ[toyTimeDomainRingdownConfig],
  Print["Configuration must define association toyTimeDomainRingdownConfig."];
  Exit[1];
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

config = toyTimeDomainRingdownConfig;
event = config["event"];
fitPath = config["fit_path"];
kerrPath = config["kerr_numeric_path"];
modeSets = config["mode_sets"];
amplitudeScenarios = config["amplitude_scenarios"];
operatorOrder = config["operators"];
branchOrder = config["polarizations"];
targetSNR = N[config["target_snr"]];
sampleRate = N[config["sample_rate_Hz"]];
duration = N[config["duration_s"]];
startTime0 = N[config["start_time_s"]];
alphaReference = N[config["alpha_reference"]];
lnMassStep = N[config["ln_mass_derivative_step"]];
spinStep = N[config["spin_derivative_step"]];
startTimeStep = N[config["start_time_derivative_step_s"]];
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

modeParameters[
  mode_,
  massSolar_,
  spin_,
  operator_ : "kerr",
  branch_ : "none",
  alpha_ : 0.
] := Module[
  {
    modeString = txt[mode], operatorString = txt[operator],
    branchString = txt[branch], massValue = num[massSolar],
    spinValue = num[spin], alphaValue = num[alpha],
    momega, shift, omegaR, omegaI
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
  <|
    "frequency_Hz" -> omegaR/(2 Pi massValue solarMassTimeSeconds),
    "tau_s" -> -massValue solarMassTimeSeconds/omegaI
  |>
];

times = N[Range[0, duration, 1/sampleRate]];
sampleCount = Length[times];
frequencies = sampleRate Range[0, sampleCount - 1]/sampleCount;
positiveFrequencyIndices = Range[2, Floor[sampleCount/2] + 1];
df = sampleRate/sampleCount;

noisePSD[f_?NumericQ] := 1 + (70/Max[f, 10])^4 + (f/500)^2;
positivePSD = noisePSD /@ frequencies[[positiveFrequencyIndices]];

coloredInnerProduct[x_List, y_List] := Module[{xf, yf},
  xf = Fourier[N[x], FourierParameters -> {1, -1}];
  yf = Fourier[N[y], FourierParameters -> {1, -1}];
  4 df Total[
    Re[
      Conjugate[xf[[positiveFrequencyIndices]]] *
        yf[[positiveFrequencyIndices]]
    ]/positivePSD
  ]
];

modeBasis[
  mode_,
  massSolar_,
  spin_,
  operator_ : "kerr",
  branch_ : "none",
  alpha_ : 0.,
  t0_ : 0.
] := Module[{pars, shiftedTimes, decay, phase},
  pars = modeParameters[mode, massSolar, spin, operator, branch, alpha];
  If[FailureQ[pars], Return[pars]];
  shiftedTimes = times - t0;
  decay = Exp[-shiftedTimes/pars["tau_s"]];
  phase = 2 Pi pars["frequency_Hz"] shiftedTimes;
  {decay Cos[phase], decay Sin[phase]}
];

basisColumns[
  modes_List,
  massSolar_,
  spin_,
  operator_ : "kerr",
  branch_ : "none",
  alpha_ : 0.,
  t0_ : 0.
] := Flatten[
  modeBasis[#, massSolar, spin, operator, branch, alpha, t0] & /@ modes,
  1
];

coefficientsForModes[modes_List, coefficientAssociation_Association] :=
  Flatten[
    If[KeyExistsQ[coefficientAssociation, txt[#]],
      num /@ coefficientAssociation[txt[#]],
      {0., 0.}
    ] & /@ modes
  ];

waveform[
  modes_List,
  coefficientAssociation_Association,
  massSolar_,
  spin_,
  operator_ : "kerr",
  branch_ : "none",
  alpha_ : 0.,
  t0_ : 0.
] := Module[{columns, coeffs},
  columns = basisColumns[modes, massSolar, spin, operator, branch, alpha, t0];
  coeffs = coefficientsForModes[modes, coefficientAssociation];
  Total[MapThread[#1 #2 &, {coeffs, columns}]]
];

profileInfo[target_List, columns_List] := Module[
  {gram, rhs, coeffs, raw, projected, chi2},
  raw = coloredInnerProduct[target, target];
  If[Length[columns] == 0,
    Return[
      <|"coefficients" -> {}, "chi2" -> raw, "raw" -> raw|>
    ]
  ];
  gram = Table[
    coloredInnerProduct[columns[[i]], columns[[j]]],
    {i, Length[columns]},
    {j, Length[columns]}
  ];
  rhs = coloredInnerProduct[#, target] & /@ columns;
  coeffs = PseudoInverse[gram].rhs;
  projected = rhs.coeffs;
  chi2 = Max[0, raw - projected];
  <|"coefficients" -> coeffs, "chi2" -> chi2, "raw" -> raw|>
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

timeDomainRows = Flatten[
  Table[
    Module[
      {
        modes = modeSet["modes"],
        modeSetName = modeSet["name"],
        coeffAssoc = amplitudeScenario["coefficients"],
        baseWave, baseNorm, scale, baseColumns, amplitudeColumns,
        lnMassDerivative, spinDerivative, startTimeDerivative,
        nuisanceColumns, alphaDerivative, rawInfo, profiled,
        profiledInfo, retained, alpha1, alpha2, alpha3
      },
      baseWave = waveform[modes, coeffAssoc, mass0, spin0, "kerr", "none", 0.,
        startTime0];
      baseNorm = coloredInnerProduct[baseWave, baseWave];
      scale = targetSNR/Sqrt[baseNorm];
      baseColumns = scale basisColumns[modes, mass0, spin0, "kerr", "none", 0.,
          startTime0];
      amplitudeColumns = baseColumns;
      lnMassDerivative = scale (
          waveform[
            modes, coeffAssoc, mass0 Exp[lnMassStep], spin0, "kerr",
            "none", 0., startTime0
          ] -
            waveform[
              modes, coeffAssoc, mass0 Exp[-lnMassStep], spin0, "kerr",
              "none", 0., startTime0
            ]
        )/(2 lnMassStep);
      spinDerivative = scale (
          waveform[
            modes, coeffAssoc, mass0, spin0 + spinStep, "kerr",
            "none", 0., startTime0
          ] -
            waveform[
              modes, coeffAssoc, mass0, spin0 - spinStep, "kerr",
              "none", 0., startTime0
            ]
        )/(2 spinStep);
      startTimeDerivative = scale (
          waveform[
            modes, coeffAssoc, mass0, spin0, "kerr", "none", 0.,
            startTime0 + startTimeStep
          ] -
            waveform[
              modes, coeffAssoc, mass0, spin0, "kerr", "none", 0.,
              startTime0 - startTimeStep
            ]
        )/(2 startTimeStep);
      nuisanceColumns = Join[
        amplitudeColumns,
        {lnMassDerivative, spinDerivative, startTimeDerivative}
      ];
      alphaDerivative = scale (
          waveform[
            modes, coeffAssoc, mass0, spin0, combo["operator"],
            combo["polarization"], alphaStep, startTime0
          ] -
            waveform[
              modes, coeffAssoc, mass0, spin0, combo["operator"],
              combo["polarization"], -alphaStep, startTime0
            ]
        )/(2 alphaStep);
      rawInfo = coloredInnerProduct[alphaDerivative, alphaDerivative];
      profiled = profileInfo[alphaDerivative, nuisanceColumns];
      profiledInfo = profiled["chi2"];
      retained = If[rawInfo <= 0, 0, profiledInfo/rawInfo];
      alpha1 = safeLimit[profiledInfo, 1.];
      alpha2 = safeLimit[profiledInfo, 4.];
      alpha3 = safeLimit[profiledInfo, 9.];
      <|
        "mode_set" -> modeSetName,
        "modes" -> StringRiffle[modes, "+"],
        "amplitude_scenario" -> amplitudeScenario["name"],
        "target_snr" -> targetSNR,
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
          alphaReference Sqrt[profiledInfo]
      |>
    ],
    {modeSet, modeSets},
    {amplitudeScenario, amplitudeScenarios},
    {combo, comboRows}
  ],
  2
];

summaryRows = Flatten[
  Table[
    Module[
      {
        rows = Select[
          timeDomainRows,
          #["mode_set"] == modeSet["name"] &&
            #["amplitude_scenario"] == amplitudeScenario["name"] &
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
        "amplitude_scenario" -> amplitudeScenario["name"],
        "target_snr" -> targetSNR,
        "best_alpha_1sigma_profiled" -> best,
        "median_alpha_1sigma_profiled" -> median,
        "operators_above_3sigma_at_alpha_reference" -> identifiable,
        "operator_count" -> Length[rows],
        "mean_retained_information_fraction" -> meanRetained
      |>
    ],
    {modeSet, modeSets},
    {amplitudeScenario, amplitudeScenarios}
  ],
  1
];

detectabilityCsvPath = FileNameJoin[
  {outputDir, "toy_time_domain_detectability.csv"}
];
summaryCsvPath = FileNameJoin[
  {outputDir, "toy_time_domain_summary.csv"}
];
plotPath = FileNameJoin[
  {outputDir, "toy_time_domain_alpha_1sigma.png"}
];
reportPath = FileNameJoin[
  {outputDir, "toy_time_domain_likelihood_report.md"}
];

exportAssociationCSV[detectabilityCsvPath, timeDomainRows];
exportAssociationCSV[summaryCsvPath, summaryRows];

plotRows = Select[
  timeDomainRows,
  #["mode_set"] == "220_221" &&
    #["amplitude_scenario"] == "moderate_221" &
];
plotLabels = (#["operator"] <> " " <> #["polarization"]) & /@ plotRows;
plotValues = Log10 /@ Lookup[plotRows, "alpha_1sigma_profiled"];

plot = BarChart[
  plotValues,
  ChartLabels -> Placed[Rotate[#, Pi/4] & /@ plotLabels, Below],
  Frame -> True,
  Axes -> False,
  FrameLabel -> {None, "log10 profiled 1-sigma coupling scale"},
  PlotLabel ->
    "Toy time-domain 220+221 identifiability after amplitude, t0, M, chi profiling",
  GridLines -> {None, Automatic},
  ImagePadding -> {{70, 40}, {160, 60}},
  ImageSize -> 1200
];
Export[plotPath, plot, ImageResolution -> 144];

summaryTableRows = Map[
  Function[row,
    "| " <> StringRiffle[
      {
        row["mode_set"],
        row["amplitude_scenario"],
        fmt[row["target_snr"], 1],
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
    timeDomainRows,
    #["mode_set"] == "220_221" &&
      #["amplitude_scenario"] == "moderate_221" &&
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
      "# Toy Time-Domain Ringdown Likelihood",
      "",
      "This report is generated by `scripts/wolfram/ringdown_toy_time_domain_likelihood.wl`.",
      "",
      "## Scope",
      "",
      "- This is still synthetic; no GW250114 strain is used.",
      "- The waveform is a sum of damped cosine/sine basis functions for the selected QNM modes.",
      "- Linear amplitudes/phases, remnant mass, remnant spin, and start time are profiled.",
      "- The inner product uses a simple analytic colored-noise PSD and rescales the Kerr injection to `SNR = " <>
        ToString[targetSNR] <> "`.",
      "- One EFT coupling is enabled at a time through the imported higher-derivative QNM shifts.",
      "",
      "## Summary",
      "",
      "| mode set | amplitude scenario | SNR | best alpha 1sigma | median alpha 1sigma | mean retained information | above 3sigma at alpha_ref |",
      "| --- | --- | ---: | ---: | ---: | ---: | ---: |"
    },
    summaryTableRows,
    {
      "",
      "## Most Detectable 220+221 Fingerprints",
      "",
      "Rows below use the `moderate_221` amplitude scenario.",
      "",
      "| operator | polarization | alpha 1sigma | retained information | Kerr-only sigma at alpha_ref |",
      "| --- | --- | ---: | ---: | ---: |"
    },
    bestTableRows,
    {
      "",
      "## Interpretation",
      "",
      "- This is a stricter test than the spectral likelihood because amplitude, phase, and start-time freedom can absorb some apparent frequency shifts.",
      "- A detectable residual now means the EFT shift changes the waveform shape in a way that is not removable by those nuisance directions.",
      "- The absolute scale still depends on the synthetic SNR, the toy PSD, and the assumed overtone excitation.",
      "- The output is therefore a design study for the real GW250114 likelihood, not an observational constraint.",
      "",
      "## Next Defensible Step",
      "",
      "Calibrate the toy settings against actual GW250114 public products: either use published ringdown posterior samples if available, or build a strain-level likelihood using GWOSC data and detector PSD estimates.",
      "",
      "## Generated Files",
      "",
      "- `toy_time_domain_detectability.csv`",
      "- `toy_time_domain_summary.csv`",
      "- `toy_time_domain_alpha_1sigma.png`"
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated toy time-domain likelihood diagnostics"];
Print["Report: ", reportPath];
Print["Detectability CSV: ", detectabilityCsvPath];
Print["mode_set\tamplitude\tbest_alpha_1sigma\tmedian_alpha_1sigma\tmean_retained_info"];
Scan[
  Print[
    StringRiffle[
      {
        #["mode_set"],
        #["amplitude_scenario"],
        fmt[#["best_alpha_1sigma_profiled"], 5],
        fmt[#["median_alpha_1sigma_profiled"], 5],
        fmt[#["mean_retained_information_fraction"], 4]
      },
      "\t"
    ]
  ] &,
  summaryRows
];
