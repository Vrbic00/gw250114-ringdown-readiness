(* Toy time-domain ringdown likelihood configuration.

   This remains synthetic. It asks how one-at-a-time EFT QNM shifts survive
   after profiling over linear mode amplitudes/phases, remnant mass, remnant
   spin, and a small start-time shift.
*)

toyTimeDomainRingdownConfig = <|
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
  "amplitude_scenarios" -> {
    <|
      "name" -> "weak_221",
      "coefficients" -> <|
        "220" -> {1.0, 0.0},
        "221" -> {0.088, -0.234}
      |>
    |>,
    <|
      "name" -> "moderate_221",
      "coefficients" -> <|
        "220" -> {1.0, 0.0},
        "221" -> {0.159, -0.313}
      |>
    |>
  },
  "operators" -> {
    "lambda_ev", "lambda_odd", "epsilon1", "epsilon2", "epsilon3"
  },
  "polarizations" -> {"plus", "minus"},
  "target_snr" -> 25.0,
  "sample_rate_Hz" -> 8192.0,
  "duration_s" -> 0.030,
  "start_time_s" -> 0.0,
  "alpha_reference" -> 0.001,
  "ln_mass_derivative_step" -> 0.001,
  "spin_derivative_step" -> 0.0005,
  "start_time_derivative_step_s" -> 0.00001,
  "alpha_derivative_step" -> 0.000001
|>;
