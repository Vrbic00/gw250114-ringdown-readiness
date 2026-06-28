(* ::Package:: *)

(* Copy this file to config/custom_cases.wl and edit the candidate rows.

   Run with:
     wolframscript -file scripts/wolfram/deviation_detectability_report.wl config/custom_cases.wl
*)

referenceMassMsun = 68.1;
ell = 2;
overtone = 0;
precisionTargetsPct = {1, 3, 5, 10};

fSchwarzschild = 1 - 2/r;

metricCases = {
  <|
    "name" -> "Schwarzschild",
    "kind" -> "GR vacuum benchmark",
    "parameters" -> "M=1",
    "reference" -> "Schwarzschild",
    "horizonRole" -> "LargestPositive",
    "f" -> fSchwarzschild
  |>,
  <|
    "name" -> "MyMetric_eps0p05",
    "kind" -> "candidate one-function metric",
    "parameters" -> "M=1, epsilon=0.05",
    "reference" -> "Schwarzschild",
    "horizonRole" -> "LargestPositive",
    "f" -> 1 - 2/r + 0.05/r^3
  |>
};

potentialCases = {
  <|
    "name" -> "Schwarzschild_ReggeWheeler_l2",
    "sector" -> "gravitational_odd",
    "level" -> "Level 5 anchor: known GR gravitational master potential",
    "description" -> "Odd-parity gravitational Regge-Wheeler potential",
    "reference" -> "Schwarzschild_ReggeWheeler_l2",
    "f" -> fSchwarzschild,
    "V" -> fSchwarzschild (ell (ell + 1)/r^2 - 6/r^3)
  |>,
  <|
    "name" -> "MyPotential_RW_deformation_l2",
    "sector" -> "gravitational_candidate",
    "level" -> "Level 4 supplied potential",
    "description" -> "Example deformation of the Regge-Wheeler potential",
    "reference" -> "Schwarzschild_ReggeWheeler_l2",
    "f" -> fSchwarzschild,
    "V" -> fSchwarzschild (ell (ell + 1)/r^2 - 6/r^3) (1 + 0.05/r^2)
  |>
};
