(* Side-by-side comparison of public GW250114 EFT projection products. *)

gw250114ConstraintsComparisonConfig = <|
  "ringdown_projection_path" -> FileNameJoin[
    {
      Directory[],
      "results", "gw250114_ringdown_eft_projection",
      "ringdown_eft_gaussian_projection.csv"
    }
  ],
  "pyring_projection_path" -> FileNameJoin[
    {
      Directory[],
      "results", "gw250114_pyring_delta_eft_projection",
      "pyring_delta_eft_projection.csv"
    }
  ],
  "ringdown_observable_summary_path" -> FileNameJoin[
    {
      Directory[],
      "results", "gw250114_ringdown_eft_projection",
      "ringdown_observable_mean_residual.csv"
    }
  ],
  "pyring_observable_summary_path" -> FileNameJoin[
    {
      Directory[],
      "results", "gw250114_pyring_delta_eft_projection",
      "pyring_delta_mean_residual.csv"
    }
  ]
|>;
