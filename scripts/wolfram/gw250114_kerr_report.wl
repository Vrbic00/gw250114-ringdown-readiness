(* ::Package:: *)

(* Reproducible Kerr QNM report for GW250114_082203.
   Usage:
     wolframscript -file scripts/wolfram/gw250114_kerr_report.wl
     wolframscript -file scripts/wolfram/gw250114_kerr_report.wl results/custom_dir
*)

ClearAll["Global`*"];

solarMassTimeSeconds = 4.925490947*10^-6;

eventParameters = <|
  "event" -> "GW250114_082203",
  "massDetectorMsun" -> 68.1,
  "massDetectorMinus" -> 0.9,
  "massDetectorPlus" -> 0.8,
  "spin" -> 0.68,
  "spinMinus" -> 0.01,
  "spinPlus" -> 0.01,
  "source" -> "GWOSC NRSur7dq4 PE, O4_Discovery_Papers/GW250114_082203/v1",
  "lvkReferenceTimeMs" -> 0.337
|>;

qnmCoefficients = <|
  "220" -> <|"l" -> 2, "m" -> 2, "n" -> 0,
    "f1" -> 1.5251, "f2" -> -1.1568, "f3" -> 0.1292,
    "q1" -> 0.7000, "q2" -> 1.4187, "q3" -> -0.4990|>,
  "221" -> <|"l" -> 2, "m" -> 2, "n" -> 1,
    "f1" -> 1.3673, "f2" -> -1.0260, "f3" -> 0.1628,
    "q1" -> 0.1000, "q2" -> 0.5436, "q3" -> -0.4731|>,
  "222" -> <|"l" -> 2, "m" -> 2, "n" -> 2,
    "f1" -> 1.3223, "f2" -> -1.0257, "f3" -> 0.1860,
    "q1" -> -0.1000, "q2" -> 0.4206, "q3" -> -0.4256|>,
  "330" -> <|"l" -> 3, "m" -> 3, "n" -> 0,
    "f1" -> 1.8956, "f2" -> -1.3043, "f3" -> 0.1818,
    "q1" -> 0.9000, "q2" -> 2.3430, "q3" -> -0.4810|>,
  "440" -> <|"l" -> 4, "m" -> 4, "n" -> 0,
    "f1" -> 2.3000, "f2" -> -1.5056, "f3" -> 0.2244,
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

kerrQNM[mode_String, massSolar_?NumericQ, spin_?NumericQ] := Module[
  {mSeconds, omegaDimless, omegaSI, fHz, q, tauSeconds},
  mSeconds = massSolar solarMassTimeSeconds;
  omegaDimless = dimensionlessAngularFrequency[mode, spin];
  omegaSI = omegaDimless/mSeconds;
  fHz = omegaSI/(2 Pi);
  q = qualityFactor[mode, spin];
  tauSeconds = 2 q/omegaSI;
  <|
    "mode" -> mode,
    "l" -> qnmCoefficients[mode]["l"],
    "m" -> qnmCoefficients[mode]["m"],
    "n" -> qnmCoefficients[mode]["n"],
    "Momega" -> omegaDimless,
    "f_Hz" -> fHz,
    "Q" -> q,
    "tau_s" -> tauSeconds,
    "tau_ms" -> 1000 tauSeconds
  |>
];

withUncertainty[mode_String, p_Association] := Module[
  {m0, mVals, s0, sVals, samples, central, fVals, tauVals},
  m0 = p["massDetectorMsun"];
  s0 = p["spin"];
  mVals = {m0 - p["massDetectorMinus"], m0, m0 + p["massDetectorPlus"]};
  sVals = {Max[0, s0 - p["spinMinus"]], s0, Min[0.999, s0 + p["spinPlus"]]};
  samples = Flatten[Table[kerrQNM[mode, m, s], {m, mVals}, {s, sVals}], 1];
  central = kerrQNM[mode, m0, s0];
  fVals = Lookup[samples, "f_Hz"];
  tauVals = Lookup[samples, "tau_ms"];
  Join[
    central,
    <|
      "f_Hz_min" -> Min[fVals],
      "f_Hz_max" -> Max[fVals],
      "tau_ms_min" -> Min[tauVals],
      "tau_ms_max" -> Max[tauVals]
    |>
  ]
];

fmt[x_?NumericQ, digits_Integer: 3] := ToString[NumberForm[N[x], {12, digits}], OutputForm];

plusMinusText[value_, min_, max_, digits_Integer: 3] := StringJoin[
  fmt[value, digits],
  " (-",
  fmt[value - min, digits],
  ", +",
  fmt[max - value, digits],
  ")"
];

args = Rest[$ScriptCommandLine];
outputDir = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[{Directory[], "results", "gw250114_kerr_qnm"}]
];

If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

modes = {"220", "221", "222", "330", "440"};
rows = withUncertainty[#, eventParameters] & /@ modes;
tMassMs = 1000 eventParameters["massDetectorMsun"] solarMassTimeSeconds;

csvPath = FileNameJoin[{outputDir, "gw250114_kerr_qnm_berti.csv"}];
plotPath = FileNameJoin[{outputDir, "gw250114_kerr_qnm_f_tau.png"}];
reportPath = FileNameJoin[{outputDir, "gw250114_kerr_qnm_report.md"}];

csvRows = Prepend[
  ({
      #["mode"], #["l"], #["m"], #["n"], #["Momega"], #["f_Hz"],
      #["f_Hz_min"], #["f_Hz_max"], #["Q"], #["tau_ms"],
      #["tau_ms_min"], #["tau_ms_max"]
    } & /@ rows),
  {
    "mode", "l", "m", "n", "Momega", "f_Hz",
    "f_Hz_min", "f_Hz_max", "Q", "tau_ms",
    "tau_ms_min", "tau_ms_max"
  }
];

Export[csvPath, csvRows, "CSV"];

plot = ListPlot[
  (Callout[{#["tau_ms"], #["f_Hz"]}, #["mode"], Above] & /@ rows),
  Frame -> True,
  FrameLabel -> {"damping time tau [ms]", "frequency f [Hz]"},
  PlotLabel -> "GW250114 Kerr QNM baseline",
  PlotRange -> All,
  PlotMarkers -> Automatic,
  GridLines -> Automatic,
  ImageSize -> 900
];

Export[plotPath, plot, ImageResolution -> 144];

markdownTableRows = ("| " <> StringRiffle[
      {
        #["mode"],
        fmt[#["Momega"], 5],
        plusMinusText[#["f_Hz"], #["f_Hz_min"], #["f_Hz_max"], 3],
        fmt[#["Q"], 4],
        plusMinusText[#["tau_ms"], #["tau_ms_min"], #["tau_ms_max"], 3]
      },
      " | "
    ] <> " |") & /@ rows;

report = StringRiffle[
  Join[
    {
      "# GW250114 Kerr QNM Baseline",
      "",
      "This report is generated by `scripts/wolfram/gw250114_kerr_report.wl`.",
      "",
      "Input event parameters:",
      "",
      "- `final_mass_detector = " <>
        ToString[eventParameters["massDetectorMsun"]] <> " -" <>
        ToString[eventParameters["massDetectorMinus"]] <> " +" <>
        ToString[eventParameters["massDetectorPlus"]] <> " Msun`",
      "- `final_spin = " <>
        ToString[eventParameters["spin"]] <> " +/- " <>
        ToString[eventParameters["spinPlus"]] <> "`",
      "- Source: " <> eventParameters["source"],
      "- Computed detector-frame mass timescale: `tM = " <>
        fmt[tMassMs, 4] <> " ms`",
      "",
      "Uncertainty ranges below are conservative corner ranges from the quoted mass and spin intervals.",
      "",
      "| mode | Momega | f_Hz central (-,+) | Q | tau_ms central (-,+) |",
      "| --- | ---: | ---: | ---: | ---: |"
    },
    markdownTableRows,
    {
      "",
      "Generated files:",
      "",
      "- `gw250114_kerr_qnm_berti.csv`",
      "- `gw250114_kerr_qnm_f_tau.png`",
      "",
      "Benchmark comparison:",
      "",
      "- The LVK GW250114 spectroscopy paper uses `tMf = 0.337 ms`; the detector-frame mass used here gives `tM = " <>
        fmt[tMassMs, 4] <> " ms`, consistent at the rounding level.",
      "- The paper reports that post-merger data support at least the 220+221 QNM content and shows 220, 221, and 440 Kerr expectations in the frequency/damping-time plane.",
      "- The present calculation reproduces only the Kerr spectral baseline from remnant mass and spin. It does not yet reproduce the LVK posterior contours or fit strain data.",
      "",
      "Interpretation:",
      "",
      "- Reproduction: Kerr QNM fitting formulas from Berti, Cardoso & Will.",
      "- Numerical check: Wolfram Language calculation with detector-frame remnant mass.",
      "- Physical interpretation: this is a baseline Kerr spectrum, not a beyond-Kerr claim.",
      "- Speculation: none."
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated Kerr QNM report for ", eventParameters["event"]];
Print["CSV: ", csvPath];
Print["Plot: ", plotPath];
Print["Report: ", reportPath];
Print["mode\tf_Hz\tf_min\tf_max\ttau_ms\ttau_min\ttau_max"];
Scan[
  Print[
    StringRiffle[
      {
        #["mode"],
        fmt[#["f_Hz"], 3],
        fmt[#["f_Hz_min"], 3],
        fmt[#["f_Hz_max"], 3],
        fmt[#["tau_ms"], 3],
        fmt[#["tau_ms_min"], 3],
        fmt[#["tau_ms_max"], 3]
      },
      "\t"
    ]
  ] &,
  rows
];
