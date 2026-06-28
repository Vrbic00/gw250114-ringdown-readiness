(* ::Package:: *)

(* Supplied master-potential QNM diagnostics.

   Master equation:
     d^2 psi/dr_*^2 + (omega^2 - V_l(r)) psi == 0

   This script applies first- and third-order WKB estimates around the peak of V_l:
     omega^2 ~= V0 - i (n + 1/2) sqrt(-2 V0'')

   at first order, plus the Iyer-Will third-order correction. All derivatives
   are with respect to tortoise coordinate.

   Usage:
     wolframscript -file scripts/wolfram/master_potential_wkb.wl

   Scope:
     Level 4: user/supplied master potentials.
     Level 5 anchor: Schwarzschild gravitational Regge-Wheeler/Zerilli
     potentials, where the master equations are known in GR.
*)

ClearAll["Global`*"];

solarMassTimeSeconds = 4.925490947*10^-6;
massDetectorMsun = 68.1;
ell = 2;
overtone = 0;

outputDir = FileNameJoin[{Directory[], "results", "master_potential_wkb"}];
If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

fSchwarzschild = 1 - 2/r;
lambdaZ[ll_] := (ll - 1) (ll + 2)/2;

potentialCases = {
  <|
    "name" -> "Schwarzschild_scalar_l2",
    "level" -> "Level 4 supplied/test-field potential",
    "description" -> "Massless scalar test-field potential on Schwarzschild",
    "f" -> fSchwarzschild,
    "V" -> fSchwarzschild (ell (ell + 1)/r^2 + D[fSchwarzschild, r]/r)
  |>,
  <|
    "name" -> "Schwarzschild_EM_l2",
    "level" -> "Level 4 supplied/test-field potential",
    "description" -> "Electromagnetic test-field potential on Schwarzschild",
    "f" -> fSchwarzschild,
    "V" -> fSchwarzschild ell (ell + 1)/r^2
  |>,
  <|
    "name" -> "Schwarzschild_ReggeWheeler_l2",
    "level" -> "Level 5 anchor: known GR gravitational master potential",
    "description" -> "Odd-parity gravitational Regge-Wheeler potential",
    "f" -> fSchwarzschild,
    "V" -> fSchwarzschild (ell (ell + 1)/r^2 - 6/r^3)
  |>,
  <|
    "name" -> "Schwarzschild_Zerilli_l2",
    "level" -> "Level 5 anchor: known GR gravitational master potential",
    "description" -> "Even-parity gravitational Zerilli potential",
    "f" -> fSchwarzschild,
    "V" -> Module[
      {lam = lambdaZ[ell]},
      2 fSchwarzschild (
        lam^2 (lam + 1) r^3 + 3 lam^2 r^2 + 9 lam r + 9
      )/(r^3 (lam r + 3)^2)
    ]
  |>,
  <|
    "name" -> "Toy_modified_RW_eps0p05_l2",
    "level" -> "Level 4 toy supplied potential",
    "description" -> "Toy 5 percent radial deformation of the Regge-Wheeler potential",
    "f" -> fSchwarzschild,
    "V" -> fSchwarzschild (ell (ell + 1)/r^2 - 6/r^3) (1 + 0.05/r^2)
  |>
};

numericPositiveRoots[expr_, variable_] := Module[
  {num, roots},
  num = Numerator[Together[expr]];
  roots = variable /. Quiet[NSolve[num == 0, variable, Reals]];
  Sort[
    DeleteDuplicates[
      Select[N[roots], NumericQ[#] && # > 0 &],
      Abs[#1 - #2] < 10^-8 &
    ]
  ]
];

firstOutside[roots_List, radius_?NumericQ, fExpr_, vExpr_] := Module[
  {candidates},
  candidates = Select[
    roots,
    # > radius + 10^-8 &&
      N[fExpr /. r -> #] > 0 &&
      N[vExpr /. r -> #] > 0 &
  ];
  If[Length[candidates] == 0, Missing["NotFound"], First[candidates]]
];

tortoiseD[expr_, fExpr_] := fExpr D[expr, r];

tortoiseDerivativeExprs[expr_, fExpr_, maxOrder_Integer] := Rest[
  NestList[tortoiseD[#, fExpr] &, expr, maxOrder]
];

safeNumber[x_] := If[MissingQ[x] || ! NumericQ[N[x]], Missing["NotAvailable"], N[x]];

fmt[x_, digits_: 5] := If[
  MissingQ[x],
  "NA",
  ToString[NumberForm[N[x], {12, digits}], OutputForm]
];

wkb3OmegaSquared[v0_, v2_, v3_, v4_, v5_, v6_, n_] := Module[
  {alpha, root, lambda2, omega2},
  If[v2 >= 0, Return[Missing["InvalidWKBPeak"]]];

  alpha = n + 1/2;
  root = Sqrt[-2 v2];

  lambda2 = (1/root) (
    (1/8) (v4/v2) (1/4 + alpha^2) -
    (1/288) (v3/v2)^2 (7 + 60 alpha^2)
  );

  omega2 = (1/(-2 v2)) (
    (5/6912) (v3/v2)^4 (77 + 188 alpha^2) -
    (1/384) ((v3^2 v4)/(v2^3)) (51 + 100 alpha^2) +
    (1/2304) (v4/v2)^2 (67 + 68 alpha^2) +
    (1/288) ((v3 v5)/(v2^2)) (19 + 28 alpha^2) -
    (1/288) (v6/v2) (5 + 4 alpha^2)
  );

  v0 + root lambda2 - I alpha root (1 + omega2)
];

omegaParts[omega_] := <|
  "real" -> If[MissingQ[omega], Missing["NotAvailable"], Re[N[omega]]],
  "imag_abs" -> If[MissingQ[omega], Missing["NotAvailable"], Abs[Im[N[omega]]]]
|>;

qualityFactor[real_, imagAbs_] := If[
  MissingQ[real] || MissingQ[imagAbs] || imagAbs == 0,
  Missing["NotAvailable"],
  real/(2 imagAbs)
];

frequencyHz[omegaReal_] := If[
  MissingQ[omegaReal],
  Missing["NotAvailable"],
  omegaReal/(2 Pi massDetectorMsun solarMassTimeSeconds)
];

dampingMs[omegaImagAbs_] := If[
  MissingQ[omegaImagAbs],
  Missing["NotAvailable"],
  1000 massDetectorMsun solarMassTimeSeconds/omegaImagAbs
];

wkbDiagnose[case_Association] := Module[
  {
    fExpr, vExpr, horizonRoots, horizonOuter, peakRoots, peakRadius,
    derivativeExprs, derivativeValues, v0, vpp, v3, v4, v5, v6,
    omegaSquared1, omega1, parts1, qFactor1, frequencyHz1, dampingMs1,
    omegaSquared3, omega3, parts3, qFactor3, frequencyHz3, dampingMs3,
    status
  },
  fExpr = case["f"];
  vExpr = case["V"];

  horizonRoots = numericPositiveRoots[fExpr, r];
  horizonOuter = If[Length[horizonRoots] > 0, Max[horizonRoots], Missing["NoPositiveHorizon"]];

  peakRoots = numericPositiveRoots[D[vExpr, r], r];
  peakRadius = If[
    MissingQ[horizonOuter],
    Missing["NoHorizon"],
    firstOutside[peakRoots, horizonOuter, fExpr, vExpr]
  ];

  v0 = If[MissingQ[peakRadius], Missing["NotAvailable"], N[vExpr /. r -> peakRadius]];
  derivativeExprs = If[
    MissingQ[peakRadius],
    Missing["NotAvailable"],
    tortoiseDerivativeExprs[vExpr, fExpr, 6]
  ];
  derivativeValues = If[
    MissingQ[derivativeExprs],
    Missing["NotAvailable"],
    N[(# /. r -> peakRadius) & /@ derivativeExprs]
  ];

  vpp = If[MissingQ[derivativeValues], Missing["NotAvailable"], derivativeValues[[2]]];
  v3 = If[MissingQ[derivativeValues], Missing["NotAvailable"], derivativeValues[[3]]];
  v4 = If[MissingQ[derivativeValues], Missing["NotAvailable"], derivativeValues[[4]]];
  v5 = If[MissingQ[derivativeValues], Missing["NotAvailable"], derivativeValues[[5]]];
  v6 = If[MissingQ[derivativeValues], Missing["NotAvailable"], derivativeValues[[6]]];

  omegaSquared1 = If[
    MissingQ[v0] || MissingQ[vpp] || vpp >= 0,
    Missing["InvalidWKBPeak"],
    v0 - I (overtone + 1/2) Sqrt[-2 vpp]
  ];
  omega1 = If[MissingQ[omegaSquared1], Missing["NotAvailable"], Sqrt[omegaSquared1]];
  parts1 = omegaParts[omega1];
  qFactor1 = qualityFactor[parts1["real"], parts1["imag_abs"]];
  frequencyHz1 = frequencyHz[parts1["real"]];
  dampingMs1 = dampingMs[parts1["imag_abs"]];

  omegaSquared3 = If[
    MissingQ[v0] || MissingQ[vpp] || MissingQ[v3] || MissingQ[v4] ||
      MissingQ[v5] || MissingQ[v6] || vpp >= 0,
    Missing["InvalidWKBPeak"],
    wkb3OmegaSquared[v0, vpp, v3, v4, v5, v6, overtone]
  ];
  omega3 = If[MissingQ[omegaSquared3], Missing["NotAvailable"], Sqrt[omegaSquared3]];
  parts3 = omegaParts[omega3];
  qFactor3 = qualityFactor[parts3["real"], parts3["imag_abs"]];
  frequencyHz3 = frequencyHz[parts3["real"]];
  dampingMs3 = dampingMs[parts3["imag_abs"]];

  status = Which[
    MissingQ[horizonOuter], "FAIL:no_positive_horizon",
    MissingQ[peakRadius], "FAIL:no_outer_positive_peak",
    MissingQ[omegaSquared1], "FAIL:invalid_wkb_peak",
    MissingQ[omegaSquared3], "WARN:WKB1_only",
    True, "PASS:WKB3_proxy"
  ];

  <|
    "name" -> case["name"],
    "level" -> case["level"],
    "description" -> case["description"],
    "status" -> status,
    "ell" -> ell,
    "n" -> overtone,
    "horizon_outer_M" -> safeNumber[horizonOuter],
    "peak_radius_M" -> safeNumber[peakRadius],
    "V0_Mminus2" -> safeNumber[v0],
    "Vpp_tortoise_Mminus4" -> safeNumber[vpp],
    "WKB1_Momega_real" -> safeNumber[parts1["real"]],
    "WKB1_Momega_imag_abs" -> safeNumber[parts1["imag_abs"]],
    "WKB1_Q" -> safeNumber[qFactor1],
    "WKB1_f_Hz_M68p1" -> safeNumber[frequencyHz1],
    "WKB1_tau_ms_M68p1" -> safeNumber[dampingMs1],
    "WKB3_Momega_real" -> safeNumber[parts3["real"]],
    "WKB3_Momega_imag_abs" -> safeNumber[parts3["imag_abs"]],
    "WKB3_Q" -> safeNumber[qFactor3],
    "WKB3_f_Hz_M68p1" -> safeNumber[frequencyHz3],
    "WKB3_tau_ms_M68p1" -> safeNumber[dampingMs3],
    "WKB_Momega_real" -> safeNumber[parts1["real"]],
    "WKB_Momega_imag_abs" -> safeNumber[parts1["imag_abs"]],
    "WKB_Q" -> safeNumber[qFactor1],
    "WKB_f_Hz_M68p1" -> safeNumber[frequencyHz1],
    "WKB_tau_ms_M68p1" -> safeNumber[dampingMs1]
  |>
];

rows = wkbDiagnose /@ potentialCases;

csvPath = FileNameJoin[{outputDir, "master_potential_wkb.csv"}];
reportPath = FileNameJoin[{outputDir, "master_potential_wkb_report.md"}];

csvFields = {
  "name", "level", "description", "status", "ell", "n",
  "horizon_outer_M", "peak_radius_M", "V0_Mminus2",
  "Vpp_tortoise_Mminus4",
  "WKB1_Momega_real", "WKB1_Momega_imag_abs", "WKB1_Q",
  "WKB1_f_Hz_M68p1", "WKB1_tau_ms_M68p1",
  "WKB3_Momega_real", "WKB3_Momega_imag_abs", "WKB3_Q",
  "WKB3_f_Hz_M68p1", "WKB3_tau_ms_M68p1",
  "WKB_Momega_real", "WKB_Momega_imag_abs",
  "WKB_Q", "WKB_f_Hz_M68p1", "WKB_tau_ms_M68p1"
};

Export[
  csvPath,
  Prepend[(Lookup[#, csvFields] & /@ rows), csvFields],
  "CSV"
];

tableRows = ("| " <> StringRiffle[
      {
        #["name"],
        #["status"],
        #["level"],
        fmt[#["peak_radius_M"], 5],
        fmt[#["WKB1_Momega_real"], 5],
        fmt[#["WKB1_Momega_imag_abs"], 5],
        fmt[#["WKB3_Momega_real"], 5],
        fmt[#["WKB3_Momega_imag_abs"], 5]
      },
      " | "
    ] <> " |") & /@ rows;

report = StringRiffle[
  Join[
    {
      "# Master-Potential WKB Diagnostics",
      "",
      "Scope: supplied one-dimensional master potentials of the form",
      "",
      "```text",
      "d^2 psi/dr_*^2 + (omega^2 - V_l(r)) psi = 0",
      "```",
      "",
      "The script reports first-order WKB and the Iyer-Will third-order WKB correction at the potential peak.",
      "",
      "First-order WKB:",
      "",
      "```text",
      "omega^2 ~= V0 - i (n + 1/2) sqrt(-2 V0'')",
      "```",
      "",
      "Here all derivatives are with respect to tortoise coordinate. This is still a screening-level diagnostic, not a replacement for continued-fraction, direct-integration, spectral, or time-domain QNM solvers.",
      "",
      "Settings: `l=2`, `n=0`, `M=68.1 Msun` for Hz/ms conversion.",
      "",
      "| potential | status | level | r_peak/M | WKB1 Re | WKB1 Abs(Im) | WKB3 Re | WKB3 Abs(Im) |",
      "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |"
    },
    tableRows,
    {
      "",
      "Interpretation:",
      "",
      "- Level 4 means a user or the code supplies a master potential, so the calculation is stronger than metric-only geodesic diagnostics.",
      "- The scalar and EM rows are test-field examples on Schwarzschild.",
      "- The Regge-Wheeler and Zerilli rows are Level 5 anchors for Schwarzschild in GR because their gravitational master equations are known.",
      "- For arbitrary alternative metrics, Level 5 is only legitimate when the underlying theory supplies linearized perturbation equations or a trusted master potential.",
      "- Third-order WKB is a numerical upgrade over first order, but low `l=2` precision still needs comparison with exact or higher-order references.",
      "",
      "Next numerical upgrade:",
      "",
      "Add a direct-integration, spectral, or time-domain evolution solver for supplied `V_l(r)`. That would turn this from a plausibility filter into a more quantitative QNM tool.",
      ""
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated master-potential WKB diagnostics"];
Print["CSV: ", csvPath];
Print["Report: ", reportPath];
Print["potential\tstatus\tr_peak\tWKB1_Re\tWKB1_ImAbs\tWKB3_Re\tWKB3_ImAbs"];
Scan[
  Print[
    StringRiffle[
      {
        #["name"],
        #["status"],
        fmt[#["peak_radius_M"], 5],
        fmt[#["WKB1_Momega_real"], 5],
        fmt[#["WKB1_Momega_imag_abs"], 5],
        fmt[#["WKB3_Momega_real"], 5],
        fmt[#["WKB3_Momega_imag_abs"], 5]
      },
      "\t"
    ]
  ] &,
  rows
];
