(* Cano et al. higher-derivative rotating-QNM bridge.

   The configured couplings are the dimensionless EFT quantities

     alpha_ev  = (ell/M)^4 lambda_ev
     alpha_odd = (ell/M)^4 lambda_odd
     alpha_i   = (ell/M)^6 epsilon_i

   from arXiv:2307.07431. The 0.01 examples are sensitivity illustrations,
   not observational constraints.
*)

eftQNMConfig = <|
  "event" -> <|
    "name" -> "GW250114_082203",
    "mass_detector_msun" -> 68.1,
    "spin" -> 0.68,
    "source" -> "GWOSC NRSur7dq4 PE, O4_Discovery_Papers/GW250114_082203/v1"
  |>,
  "coefficient_path" -> FileNameJoin[
    {Directory[], "data", "qnm_higher_derivative_coefficients.csv"}
  ],
  "validation_path" -> FileNameJoin[
    {Directory[], "data", "qnm_higher_derivative_validation_chi0p7.csv"}
  ],
  "kerr_numeric_path" -> FileNameJoin[
    {
      Directory[],
      "results",
      "gw250114_kerr_qnm",
      "qnm_solver_crosscheck.csv"
    }
  ],
  "example_alpha" -> 0.01,
  "scenarios" -> {
    <|
      "name" -> "lambda_ev_only",
      "couplings" -> <|"lambda_ev" -> 0.01|>
    |>,
    <|
      "name" -> "lambda_odd_only",
      "couplings" -> <|"lambda_odd" -> 0.01|>
    |>,
    <|
      "name" -> "epsilon1_only",
      "couplings" -> <|"epsilon1" -> 0.01|>
    |>,
    <|
      "name" -> "epsilon2_only",
      "couplings" -> <|"epsilon2" -> 0.01|>
    |>,
    <|
      "name" -> "epsilon3_only",
      "couplings" -> <|"epsilon3" -> 0.01|>
    |>
  }
|>;
