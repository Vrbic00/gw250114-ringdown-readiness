(* ::Package:: *)

(* Editable case registry for static spherical metric diagnostics.

   Conventions:
   - r is dimensionless radius in units of M.
   - one-function metrics use ds^2 = -f(r) dt^2 + dr^2/f(r) + r^2 dOmega^2.
   - supplied potentials use the same tortoise coordinate convention
     dr_*/dr = 1/f(r).
*)

referenceMassMsun = 68.1;
ell = 2;
overtone = 0;
precisionTargetsPct = {1, 3, 5, 10};

fSchwarzschild = 1 - 2/r;
lambdaZ[ll_] := (ll - 1) (ll + 2)/2;

decimalTag[x_] := StringReplace[ToString[NumberForm[x, {4, 3}, NumberPoint -> "p"]], " " -> ""];

metricCases = Join[
  {
  <|
    "name" -> "Schwarzschild",
    "kind" -> "GR vacuum benchmark",
    "parameters" -> "M=1",
    "reference" -> "Schwarzschild",
    "horizonRole" -> "LargestPositive",
    "f" -> fSchwarzschild
  |>
  },
  Table[
    <|
      "name" -> "ReissnerNordstrom_q" <> decimalTag[q],
      "kind" -> "GR electrovac validation case",
      "parameters" -> "M=1, Q/M=" <> ToString[q],
      "reference" -> "Schwarzschild",
      "horizonRole" -> "LargestPositive",
      "f" -> 1 - 2/r + q^2/r^2
    |>,
    {q, {0.2, 0.5, 0.8, 0.95}}
  ],
  Table[
    <|
      "name" -> "Bardeen_g" <> decimalTag[g],
      "kind" -> "regular black-hole candidate metric",
      "parameters" -> "M=1, g/M=" <> ToString[g],
      "reference" -> "Schwarzschild",
      "horizonRole" -> "LargestPositive",
      "f" -> 1 - 2 r^2/(r^2 + g^2)^(3/2)
    |>,
    {g, {0.2, 0.4, 0.6}}
  ],
  Table[
    <|
      "name" -> "Hayward_l" <> decimalTag[l],
      "kind" -> "regular black-hole candidate metric",
      "parameters" -> "M=1, l/M=" <> ToString[l],
      "reference" -> "Schwarzschild",
      "horizonRole" -> "LargestPositive",
      "f" -> 1 - 2 r^2/(r^3 + 2 l^2)
    |>,
    {l, {0.2, 0.4, 0.6}}
  ],
  Table[
    With[
      {lam = pair[[1]], tag = pair[[2]]},
      <|
        "name" -> "SchwarzschildDeSitter_lambda" <> tag,
        "kind" -> "GR Lambda validation case",
        "parameters" -> "M=1, Lambda M^2=" <> ToString[N[lam]],
        "reference" -> "Schwarzschild",
        "horizonRole" -> "SmallestPositive",
        "f" -> 1 - 2/r - lam r^2/3
      |>
    ],
    {pair, {{1/100000, "1em5"}, {1/10000, "1em4"}, {1/2000, "5em4"}}}
  ],
  {
    <|
      "name" -> "ToyDeformation_eps0p05",
      "kind" -> "toy Schwarzschild deformation",
      "parameters" -> "M=1, epsilon=0.05, f=1-2/r+epsilon/r^3",
      "reference" -> "Schwarzschild",
      "horizonRole" -> "LargestPositive",
      "f" -> 1 - 2/r + 0.05/r^3
    |>
  }
];

potentialCases = {
  <|
    "name" -> "Schwarzschild_scalar_l2",
    "sector" -> "scalar",
    "level" -> "Level 4 supplied/test-field potential",
    "description" -> "Massless scalar test-field potential on Schwarzschild",
    "reference" -> "Schwarzschild_scalar_l2",
    "f" -> fSchwarzschild,
    "V" -> fSchwarzschild (ell (ell + 1)/r^2 + D[fSchwarzschild, r]/r)
  |>,
  <|
    "name" -> "Schwarzschild_EM_l2",
    "sector" -> "electromagnetic",
    "level" -> "Level 4 supplied/test-field potential",
    "description" -> "Electromagnetic test-field potential on Schwarzschild",
    "reference" -> "Schwarzschild_EM_l2",
    "f" -> fSchwarzschild,
    "V" -> fSchwarzschild ell (ell + 1)/r^2
  |>,
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
    "name" -> "Schwarzschild_Zerilli_l2",
    "sector" -> "gravitational_even",
    "level" -> "Level 5 anchor: known GR gravitational master potential",
    "description" -> "Even-parity gravitational Zerilli potential",
    "reference" -> "Schwarzschild_ReggeWheeler_l2",
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
    "sector" -> "gravitational_toy",
    "level" -> "Level 4 toy supplied potential",
    "description" -> "Toy 5 percent radial deformation of the Regge-Wheeler potential",
    "reference" -> "Schwarzschild_ReggeWheeler_l2",
    "f" -> fSchwarzschild,
    "V" -> fSchwarzschild (ell (ell + 1)/r^2 - 6/r^3) (1 + 0.05/r^2)
  |>
};
