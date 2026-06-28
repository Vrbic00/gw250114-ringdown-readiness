(* Approximate EFT projection using public GW250114 pyRing 221 deviation samples. *)

gw250114PyRingDeltaEFTProjectionConfig = <|
  "event" -> <|
    "name" -> "GW250114_082203",
    "mass_detector_msun" -> 68.1,
    "spin" -> 0.68
  |>,
  "pyring_delta_path" -> FileNameJoin[
    {
      Directory[],
      "data", "raw", "GW250114_data_release", "data",
      "posterior_with_qnm_frequencies.dat"
    }
  ],
  "eft_sensitivity_path" -> FileNameJoin[
    {
      Directory[],
      "results",
      "higher_derivative_qnm_complete",
      "gw250114_complete_eft_sensitivities.csv"
    }
  ],
  "pyring_delta_f_plot_bound" -> 0.8,
  "mode" -> "221"
|>;
