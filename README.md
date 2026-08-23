# Repository Structure

```text
.
├── literature_review/                # Some literature review materials
│
├── simulation/                       # Monte Carlo simulation study
│   ├── meta-Cox_model_-_Sim.R        # Data generation, SCAD/MCP hierarchical fits,
│   │                                 # competitor methods, and full simulation study (R)
│   ├── sim_results_full.csv          # Per-replication results for the full scenario grid
│   │                                 # (nsim = 500)
│   ├── sim_results_full.rds          # Same results, saved as an R object
│   ├── sim_summary_full.csv          # Replication-averaged summary statistics by
│   │                                 # scenario and method
│   │
│   ├── plot_MCC.png                  # MCC by scenario (small demo simulation)
│   ├── plot_bias.png                 # Alpha bias by scenario (small demo simulation)
│   ├── plot_FDR.png                  # FDR by scenario (small demo simulation)
│   ├── plot_MSE.png                  # MSE(theta) by scenario (small demo simulation)
│   ├── plot_coef_SCAD.png            # SCAD (Hierarchical): estimated vs true theta_kj
│   │                                 # by study
│   ├── plot_coef_MCP.png             # MCP (Hierarchical): estimated vs true theta_kj
│   │                                 # by study
│   ├── plot_het_SCAD.png             # SCAD: study-specific deviations, true vs estimated
│   ├── plot_het_MCP.png              # MCP: study-specific deviations, true vs estimated
│   │
│   ├── plot_full_MCC.png             # MCC across all scenarios (full simulation)
│   ├── plot_full_TPR.png             # TPR across all scenarios (full simulation)
│   ├── plot_full_FDR.png             # FDR across all scenarios (full simulation)
│   ├── plot_full_F1.png              # F1 across all scenarios (full simulation)
│   ├── plot_full_alpha_bias.png      # Alpha bias across all scenarios (full simulation)
│   ├── plot_full_fp_bias.png         # fp_bias across all scenarios (full simulation)
│   ├── plot_full_MSE_theta.png       # MSE(theta) across all scenarios (full simulation)
│   └── plot_full_exact_match.png     # Exact active-set match rate across all scenarios
│                                     # (full simulation)
│
└── application/                      # Application to curatedOvarianData, a multi-cohort
    │                                 # ovarian cancer survival dataset
    ├── meta-Cox_model_-_Anl.R        # Loads and preprocesses curatedOvarianData, fits
    │                                 # hierarchical and competitor methods, and produces
    │                                 # all tables and figures below (R)
    ├── Table7_1_cohort_chars.csv     # Cohort-level characteristics (N, events, censoring, etc.)
    ├── Table7_2_biomarkers.csv       # Selected biomarkers (genes) and estimated effects
    ├── Table7_2_biomarkers_CI.csv    # Selected biomarkers with confidence intervals
    ├── Table7_3_method_summary.csv   # Selected gene-set sizes and summary metrics by method
    │
    ├── Fig7_1_cohort_bar.png         # Cohort composition: events vs censored patients
    │                                 # per study
    ├── Fig7_2_kepm.png              # Kaplan-Meier survival curves by cohort
    ├── Fig7_3_forest.png             # Global vs study-specific log-hazard ratios
    │                                 # (forest plot)
    ├── Fig7_4_deviation_tile.png     # Study-specific deviations from global effects
    │                                 # (heatmap)
    ├── Fig7_5_mcp_vs_scad.png        # MCP vs SCAD hierarchical study-specific estimates
    └── Fig7_6_selection_sizes.png    # Selected gene-set sizes by method
