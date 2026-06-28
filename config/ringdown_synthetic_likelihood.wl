(* Synthetic ringdown identifiability configuration.

   This is not an observational GW250114 likelihood. It assigns illustrative
   measurement widths to log f and log tau, then asks how much of a one-at-a-
   time EFT QNM shift survives profiling over remnant mass and spin.
*)

syntheticRingdownLikelihoodConfig = <|
  "event" -> <|
    "name" -> "GW250114_082203",
    "mass_detector_msun" -> 68.1,
    "spin" -> 0.68,
    "source" -> "GWOSC NRSur7dq4 PE, O4_Discovery_Papers/GW250114_082203/v1"
  |>,
  "fit_path" -> FileNameJoin[
    {Directory[], "data", "beyond_kerr_qnm_selected_fits.csv"}
  ],
  "kerr_numeric_path" -> FileNameJoin[
    {
      Directory[],
      "results",
      "gw250114_kerr_qnm",
      "qnm_solver_crosscheck.csv"
    }
  ],
  "mode_sets" -> {
    <|"name" -> "220_only", "modes" -> {"220"}|>,
    <|"name" -> "220_221", "modes" -> {"220", "221"}|>
  },
  "operators" -> {
    "lambda_ev", "lambda_odd", "epsilon1", "epsilon2", "epsilon3"
  },
  "polarizations" -> {"plus", "minus"},
  "uncertainty_scenarios" -> {
    <|
      "name" -> "optimistic_1pct_f_5pct_tau",
      "sigma_lnf" -> 0.01,
      "sigma_lntau" -> 0.05
    |>,
    <|
      "name" -> "moderate_2pct_f_10pct_tau",
      "sigma_lnf" -> 0.02,
      "sigma_lntau" -> 0.10
    |>,
    <|
      "name" -> "conservative_5pct_f_20pct_tau",
      "sigma_lnf" -> 0.05,
      "sigma_lntau" -> 0.20
    |>
  },
  "alpha_reference" -> 0.001,
  "spin_derivative_step" -> 0.0005,
  "alpha_derivative_step" -> 0.000001,
  "confusion_mode_set" -> "220_221",
  "confusion_scenario" -> "moderate_2pct_f_10pct_tau"
|>;
