(* Posterior-informed synthetic EFT spectral projection.

   This is still not an EFT observation. It compares free remnant profiling
   with Gaussian remnant priors inferred from the public NRSur7dq4 full-IMR
   posterior.
*)

posteriorInformedProjectionConfig = <|
  "event" -> <|
    "name" -> "GW250114_082203",
    "mass_detector_msun" -> 68.1,
    "spin" -> 0.68,
    "source" -> "GWOSC O4_Discovery_Papers/GW250114_082203/v1"
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
  "posterior_samples_path" -> FileNameJoin[
    {
      Directory[],
      "results",
      "gw250114_posterior_calibration",
      "nrSur7dq4_selected_posterior_samples.csv"
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
  "remnant_prior_scenarios" -> {
    <|"name" -> "free_remnant", "scale" -> Infinity|>,
    <|"name" -> "nrSur7dq4_imr_prior", "scale" -> 1.0|>,
    <|"name" -> "loose_3x_imr_prior", "scale" -> 3.0|>
  },
  "alpha_reference" -> 0.001,
  "spin_derivative_step" -> 0.0005
|>;
