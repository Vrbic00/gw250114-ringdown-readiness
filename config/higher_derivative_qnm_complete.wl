(* Complete selected-mode EFT QNM spectrum from Cano et al. (2024).

   The example coupling is deliberately small because overtones can be much
   more sensitive than fundamental modes. It is not an observational bound.
*)

completeEFTQNMConfig = <|
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
  "modes" -> {"220", "221", "222", "330", "440"},
  "example_alpha" -> 0.001
|>;
