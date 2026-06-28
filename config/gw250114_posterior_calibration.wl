(* GW250114 public posterior calibration.

   The posterior file is the preferred NRSur7dq4 PE posterior from the GWOSC
   event API / Zenodo data release. It is not a ringdown-only posterior.
*)

gw250114PosteriorCalibrationConfig = <|
  "event" -> <|
    "name" -> "GW250114_082203",
    "posterior_model" -> "NRSur7dq4",
    "source" -> "GWOSC O4_Discovery_Papers/GW250114_082203/v1",
    "zenodo_record" -> "10.5281/zenodo.16877101",
    "zenodo_file" -> "posterior_samples_NRSur7dq4.h5",
    "posterior_dataset" ->
      "/bilby-NRSur7dq4_prod-reweighted/posterior_samples"
  |>,
  "posterior_path" -> FileNameJoin[
    {Directory[], "data", "raw", "posterior_samples_NRSur7dq4.h5"}
  ],
  "modes" -> {"220", "221", "222", "330", "440"},
  "selected_columns" -> {
    "log_likelihood",
    "final_mass",
    "final_mass_source",
    "final_spin",
    "redshift",
    "radiated_energy",
    "mass_1",
    "mass_2",
    "mass_ratio",
    "chi_eff",
    "chi_p",
    "network_optimal_snr",
    "network_matched_filter_snr",
    "network_21_multipole_snr",
    "network_33_multipole_snr",
    "network_44_multipole_snr"
  },
  "summary_columns" -> {
    "final_mass",
    "final_mass_source",
    "final_spin",
    "redshift",
    "radiated_energy",
    "mass_1",
    "mass_2",
    "mass_ratio",
    "chi_eff",
    "chi_p",
    "network_optimal_snr",
    "network_matched_filter_snr",
    "network_21_multipole_snr",
    "network_33_multipole_snr",
    "network_44_multipole_snr"
  },
  "covariance_modes" -> {"220", "221"},
  "central_reference" -> <|
    "mass_detector_msun" -> 68.1,
    "spin" -> 0.68
  |>
|>;
