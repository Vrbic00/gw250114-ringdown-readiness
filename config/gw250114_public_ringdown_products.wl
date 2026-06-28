(* Public GW250114 ringdown data products from the Zenodo release tarball. *)

gw250114PublicRingdownConfig = <|
  "release_root" -> FileNameJoin[
    {Directory[], "data", "raw", "GW250114_data_release", "data"}
  ],
  "pyring_220_posterior" -> "pyring_220_posterior.dat",
  "pyring_220_221_delta_posterior" -> "posterior_with_qnm_frequencies.dat",
  "ringdown_220_221_delta_posterior" ->
    "220+221+df221+dg221_6M_f220meas_f221meas_df221meas_120Ksamps.hdf5",
  "eft_sensitivity_path" -> FileNameJoin[
    {
      Directory[],
      "results",
      "higher_derivative_qnm_complete",
      "gw250114_complete_eft_sensitivities.csv"
    }
  ],
  "pyring_delta_f_plot_bound" -> 0.8,
  "alpha_proxy_mode" -> "221"
|>;
