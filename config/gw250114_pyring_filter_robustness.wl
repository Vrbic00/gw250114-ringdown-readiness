(* Robustness sweep for public GW250114 pyRing 221 deviation projection. *)

gw250114PyRingFilterRobustnessConfig = <|
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
  "mode" -> "221",
  "filter_scenarios" -> {
    <|
      "name" -> "strict_df_bound_0p5",
      "description" -> "stricter lower-tail cut than the public Figure 4 script",
      "delta_f_bound" -> 0.5
    |>,
    <|
      "name" -> "public_df_bound_0p8",
      "description" -> "same lower-tail cut as the public Figure 4 script",
      "delta_f_bound" -> 0.8
    |>,
    <|
      "name" -> "loose_df_bound_1p2",
      "description" -> "looser lower-tail cut than the public Figure 4 script",
      "delta_f_bound" -> 1.2
    |>,
    <|
      "name" -> "positive_domain_only",
      "description" -> "only require positive log-domain variables",
      "delta_f_bound" -> Infinity
    |>
  }
|>;
