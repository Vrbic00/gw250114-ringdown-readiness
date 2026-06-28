(* ::Package:: *)

(* Static spherical metric diagnostics for one-function metrics

   ds^2 = -f(r) dt^2 + dr^2/f(r) + r^2 dOmega^2

   This is not a full gravitational perturbation solver. It provides
   background sanity checks, geodesic observables, and an eikonal QNM proxy.

   Usage:
     wolframscript -file scripts/wolfram/spherical_metric_diagnostics.wl
*)

ClearAll["Global`*"];

solarMassTimeSeconds = 4.925490947*10^-6;

args = Rest[$ScriptCommandLine];
configPath = If[
  Length[args] >= 1,
  args[[1]],
  FileNameJoin[{Directory[], "config", "spherical_cases.wl"}]
];

If[FileExistsQ[configPath],
  Get[configPath],
  referenceMassMsun = 68.1;
  ell = 2;
  overtone = 0;
  metricCases = {}
];

If[! ValueQ[referenceMassMsun], referenceMassMsun = 68.1];
If[! ValueQ[ell], ell = 2];
If[! ValueQ[overtone], overtone = 0];

massDetectorMsun = referenceMassMsun;

outputDir = FileNameJoin[{Directory[], "results", "spherical_metric_diagnostics"}];
If[! DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]
];

metrics = metricCases;

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

firstOutside[roots_List, horizon_?NumericQ, fExpr_] := Module[
  {candidates},
  candidates = Select[roots, # > horizon + 10^-8 && N[fExpr /. r -> #] > 0 &];
  If[Length[candidates] == 0, Missing["NotFound"], First[candidates]]
];

selectHorizon[roots_List, role_String:"LargestPositive"] := Which[
  Length[roots] == 0, Missing["NoPositiveHorizon"],
  role == "SmallestPositive", Min[roots],
  True, Max[roots]
];

safeNumber[x_] := If[MissingQ[x] || ! NumericQ[N[x]], Missing["NotAvailable"], N[x]];

fmt[x_, digits_: 4] := If[
  MissingQ[x],
  "NA",
  ToString[NumberForm[N[x], {12, digits}], OutputForm]
];

diagnoseMetric[metric_Association] := Module[
  {
    fExpr, horizonRoots, horizonOuter, photonRoots, photonRadius,
    fPhoton, omegaC, lambdaSquared, lambda, l2Expr, iscoRoots, iscoRadius,
    mOmegaReal, mOmegaImagAbs, frequencyHz, dampingMs, status, horizonRole
  },
  fExpr = metric["f"];
  horizonRole = Lookup[metric, "horizonRole", "LargestPositive"];

  horizonRoots = numericPositiveRoots[fExpr, r];
  horizonOuter = selectHorizon[horizonRoots, horizonRole];

  photonRoots = numericPositiveRoots[r D[fExpr, r] - 2 fExpr, r];
  photonRadius = If[
    MissingQ[horizonOuter],
    Missing["NoHorizon"],
    firstOutside[photonRoots, horizonOuter, fExpr]
  ];

  fPhoton = If[MissingQ[photonRadius], Missing["NotAvailable"], N[fExpr /. r -> photonRadius]];
  omegaC = If[MissingQ[photonRadius], Missing["NotAvailable"], Sqrt[fPhoton/photonRadius^2]];
  lambdaSquared = If[
    MissingQ[photonRadius],
    Missing["NotAvailable"],
    N[fPhoton (2 fPhoton - photonRadius^2 (D[fExpr, {r, 2}] /. r -> photonRadius))/(2 photonRadius^2)]
  ];
  lambda = If[MissingQ[lambdaSquared] || lambdaSquared <= 0, Missing["UnstableOrInvalid"], Sqrt[lambdaSquared]];

  l2Expr = Together[r^3 D[fExpr, r]/(2 fExpr - r D[fExpr, r])];
  iscoRoots = numericPositiveRoots[D[l2Expr, r], r];
  iscoRadius = If[
    MissingQ[photonRadius],
    Missing["NotAvailable"],
    firstOutside[iscoRoots, photonRadius, fExpr]
  ];

  mOmegaReal = If[MissingQ[omegaC], Missing["NotAvailable"], ell omegaC];
  mOmegaImagAbs = If[MissingQ[lambda], Missing["NotAvailable"], (overtone + 1/2) lambda];

  frequencyHz = If[
    MissingQ[mOmegaReal],
    Missing["NotAvailable"],
    mOmegaReal/(2 Pi massDetectorMsun solarMassTimeSeconds)
  ];
  dampingMs = If[
    MissingQ[mOmegaImagAbs],
    Missing["NotAvailable"],
    1000 massDetectorMsun solarMassTimeSeconds/mOmegaImagAbs
  ];

  status = Which[
    MissingQ[horizonOuter], "FAIL:no_positive_horizon",
    MissingQ[photonRadius], "FAIL:no_outer_photon_sphere",
    MissingQ[lambda], "FAIL:no_valid_photon_instability",
    MissingQ[iscoRadius], "WARN:no_outer_isco_found",
    True, "PASS:background_proxy"
  ];

  <|
    "name" -> metric["name"],
    "kind" -> metric["kind"],
    "parameters" -> metric["parameters"],
    "horizon_role" -> horizonRole,
    "status" -> status,
    "horizon_outer_M" -> safeNumber[horizonOuter],
    "photon_sphere_M" -> safeNumber[photonRadius],
    "isco_M" -> safeNumber[iscoRadius],
    "Omega_c_M" -> safeNumber[omegaC],
    "lambda_M" -> safeNumber[lambda],
    "eikonal_Momega_real_l2" -> safeNumber[mOmegaReal],
    "eikonal_Momega_imag_abs_n0" -> safeNumber[mOmegaImagAbs],
    "eikonal_f_Hz_M68p1" -> safeNumber[frequencyHz],
    "eikonal_tau_ms_M68p1" -> safeNumber[dampingMs]
  |>
];

rows = diagnoseMetric /@ metrics;

csvPath = FileNameJoin[{outputDir, "spherical_metric_diagnostics.csv"}];
reportPath = FileNameJoin[{outputDir, "spherical_metric_diagnostics_report.md"}];

csvFields = {
  "name", "kind", "parameters", "horizon_role", "status", "horizon_outer_M",
  "photon_sphere_M", "isco_M", "Omega_c_M", "lambda_M",
  "eikonal_Momega_real_l2", "eikonal_Momega_imag_abs_n0",
  "eikonal_f_Hz_M68p1", "eikonal_tau_ms_M68p1"
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
        fmt[#["horizon_outer_M"], 5],
        fmt[#["photon_sphere_M"], 5],
        fmt[#["isco_M"], 5],
        fmt[#["eikonal_Momega_real_l2"], 5],
        fmt[#["eikonal_Momega_imag_abs_n0"], 5],
        fmt[#["eikonal_f_Hz_M68p1"], 3],
        fmt[#["eikonal_tau_ms_M68p1"], 3]
      },
      " | "
    ] <> " |") & /@ rows;

report = StringRiffle[
  Join[
    {
      "# Spherical Metric Diagnostics",
      "",
      "Scope: one-function static spherical metrics",
      "",
      "```text",
      "ds^2 = -f(r) dt^2 + dr^2/f(r) + r^2 dOmega^2",
      "```",
      "",
      "This is a background and eikonal-proxy diagnostic, not a full gravitational perturbation solver.",
      "",
      "Checks included:",
      "",
      "- outer positive horizon from `f(r)=0`,",
      "- outer photon sphere from `r f'(r) - 2 f(r)=0`,",
      "- ISCO proxy from `d L^2/dr = 0` for timelike circular orbits,",
      "- eikonal QNM proxy `M omega ~= l M Omega_c - i (n+1/2) M lambda`.",
      "",
      "For the physical Hz/ms conversion below, the mass scale is `M = 68.1 Msun`, matching the GW250114 detector-frame remnant mass used in the Kerr baseline.",
      "",
      "| metric | status | r_h/M | r_ph/M | r_isco/M | Re(M omega) | Abs(Im(M omega)) | f [Hz] | tau [ms] |",
      "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
    },
    tableRows,
    {
      "",
      "Interpretation:",
      "",
      "- `PASS:background_proxy` means the metric passes only this low-level background/eikonal diagnostic.",
      "- This does not imply that the underlying theory has healthy gravitational perturbations.",
      "- A rigorous gravitational-QNM test requires the theory's perturbation equations or a supplied master potential.",
      "- The Schwarzschild row is the validation anchor: `r_h=2M`, `r_ph=3M`, `r_isco=6M`, and `M Omega_c = M lambda = 1/(3 sqrt(3))`.",
      "",
      "Recommended next level:",
      "",
      "Allow a user-supplied potential `V_l(r)` and compute test-field or supplied-master-equation QNMs. That would be a physically stronger online diagnostic than using the metric alone.",
      ""
    }
  ],
  "\n"
];

Export[reportPath, report, "Text"];

Print["Generated spherical metric diagnostics"];
Print["CSV: ", csvPath];
Print["Report: ", reportPath];
Print["metric\tstatus\tr_h\tr_ph\tr_isco\tMomega_Re\tMomega_ImAbs"];
Scan[
  Print[
    StringRiffle[
      {
        #["name"],
        #["status"],
        fmt[#["horizon_outer_M"], 5],
        fmt[#["photon_sphere_M"], 5],
        fmt[#["isco_M"], 5],
        fmt[#["eikonal_Momega_real_l2"], 5],
        fmt[#["eikonal_Momega_imag_abs_n0"], 5]
      },
      "\t"
    ]
  ] &,
  rows
];
