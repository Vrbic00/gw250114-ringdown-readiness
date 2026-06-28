(* Curated paper-facing table build for the GW250114 public projection project. *)

gw250114PaperTablesConfig = <|
  "constraints_long_path" -> FileNameJoin[
    {
      Directory[], "results", "gw250114_constraints_comparison",
      "projection_constraints_long.csv"
    }
  ],
  "consistency_path" -> FileNameJoin[
    {
      Directory[], "results", "gw250114_constraints_comparison",
      "projection_consistency_by_operator.csv"
    }
  ],
  "filter_summary_path" -> FileNameJoin[
    {
      Directory[], "results", "gw250114_pyring_filter_robustness",
      "pyring_filter_scenario_summary.csv"
    }
  ],
  "linearized_summary_path" -> FileNameJoin[
    {
      Directory[], "results", "gw250114_linearized_posterior_projection",
      "linearized_projection_summary.csv"
    }
  ],
  "public_ringdown_summary_path" -> FileNameJoin[
    {
      Directory[], "results", "gw250114_public_ringdown_products",
      "public_ringdown_product_summary.csv"
    }
  ]
|>;
