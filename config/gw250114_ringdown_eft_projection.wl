(* Approximate EFT projection using public GW250114 RINGDOWN posterior samples. *)

gw250114RingdownEFTProjectionConfig = <|
  "event" -> <|
    "name" -> "GW250114_082203",
    "mass_detector_msun" -> 68.1,
    "spin" -> 0.68
  |>,
  "ringdown_hdf5_path" -> FileNameJoin[
    {
      Directory[],
      "data", "raw", "GW250114_data_release", "data",
      "220+221+df221+dg221_6M_f220meas_f221meas_df221meas_120Ksamps.hdf5"
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
  "kerr_numeric_path" -> FileNameJoin[
    {
      Directory[],
      "results",
      "gw250114_kerr_qnm",
      "qnm_solver_crosscheck.csv"
    }
  ],
  "modes" -> {"220", "221"}
|>;
