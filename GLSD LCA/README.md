# Respond Surface Methodology (AnD feed)

Slimmed bundle of the files needed to reproduce the anaerobic-digester-feed response-surface surrogate and its diagnostics. The Julia environment in this folder is self contained.

## Included
- Inputs: `AD_database.db`, digester temperature workbooks (`Jan-Dec 2021...xlsx`, `Jan-July 2022...xlsx`), cleaned inventory/LCA tables (`cleaned_data_per_cubic_meter_AnD_feed_in.csv`, `lcia_results_per_cubic_meter_AnD_feed_in.csv`, `raw_inputs_with_feed.csv`, `panel_FU_digester_feed.csv`).
- Scripts: `LCI_Calculation_FU_PerCubicMeter_AD_feedin V2.jl`, `RSM_FU_digester_feed.jl`, `RSM_Layer1_V7 AD Feed In.jl`, `Raw_to_Impact_Surrogate_2_fixed.jl`, `Raw_Ratio_To_Impact_2_fixed.jl`, `Sensitivity_Analysis_Beta_SSO_AnD_feed_in.jl`, `Transformer_Model.jl`.
- Outputs/results: model summaries, VIF/sensitivity tables (`model_full_summary_per_cubic_meter_AnD_feed_in.csv`, `model_summary_per_cubic_meter_AnD_feed_in.csv`, `vif_result_per_cubic_meter_AnD_feed_in.csv`, `beta_SSO_sensitivity_per_cubic_meter_AnD_feed_in.csv`, `SigmaE_cov_test.csv`, `SigmaE_cor_test.csv`, prediction CSVs), plus figures in `plots_per_cubic_meter_AnD_feed_in/` and `plots_per_AD_feed_in/`.

## How to run
From this folder:
1) Install deps once: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
2) (Optional) Rebuild LCI predictors from raw SCADA/DB: `julia --project=. "LCI_Calculation_FU_PerCubicMeter_AD_feedin V2.jl"`
3) Build the panel and basic OLS checks: `julia --project=. RSM_FU_digester_feed.jl`
4) Train the multi-response surrogate + plots/tables: `julia --project=. "RSM_Layer1_V7 AD Feed In.jl"`
5) Apply surrogate to raw inventories (OpenLCA residual checks): `julia --project=. Raw_to_Impact_Surrogate_2_fixed.jl`
6) Optional alternate raw-ratio surrogate: `julia --project=. Raw_Ratio_To_Impact_2_fixed.jl`
7) Optional SSO sensitivity table: `julia --project=. Sensitivity_Analysis_Beta_SSO_AnD_feed_in.jl`

The scripts read/write files relative to this directory; plots land in `plots_per_cubic_meter_AnD_feed_in/` and `plots_per_AD_feed_in/`. The provided CSVs/figures already reflect a clean run; re-running will overwrite them.
