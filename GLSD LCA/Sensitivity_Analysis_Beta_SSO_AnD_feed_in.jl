using CSV
using DataFrames
using Statistics

# ---------------------------
# 1. Filenames (per your current workflow)
# ---------------------------
model_file = "model_full_summary_per_cubic_meter_AnD_feed_in.csv"
lcia_file  = "lcia_results_per_cubic_meter_AnD_feed_in.csv"
# This is your LCI / predictor table and must include the SSO_LCI column
x_file     = "cleaned_data_dedup_per_cubic_meter_AnD_feed_in.csv"

# Output sensitivity results
out_file   = "beta_SSO_sensitivity_per_cubic_meter_AnD_feed_in.csv"

# ---------------------------
# 2. Read data
# ---------------------------
model_df = CSV.read(model_file, DataFrame)
x_df     = CSV.read(x_file, DataFrame)
lcia_df  = CSV.read(lcia_file, DataFrame)

@assert nrow(x_df) == nrow(lcia_df) "LCI and LCIA row counts differ; please check alignment"

# SSO predictor (already the per-m³ AnD feed LCI)
@assert :SSO_LCI in propertynames(x_df)
x_sso = x_df.SSO_LCI

# All impact indicator names, read from the regression summary
indicators = unique(model_df.Indicator)

# ---------------------------
# 3. Prepare result table
# ---------------------------
results = DataFrame(
    Indicator             = String[],
    beta_SSO              = Float64[],
    SE                    = Float64[],
    CI_low                = Float64[],
    CI_high               = Float64[],
    StdCoeff              = Float64[],
    Elasticity            = Float64[],
    DeltaY_10pct          = Float64[],
    RelDeltaY_10pct_pct   = Float64[]
)

# ---------------------------
# 4. Compute sensitivity for each indicator
# ---------------------------
for ind in indicators
    # 4.1 Grab the coefficient and standard error of SSO_LCI from the regression summary
    row = filter(r -> r.Indicator == ind && r.Variable == "SSO_LCI", model_df)

    if nrow(row) == 0
        @warn "Indicator $ind does not have an SSO_LCI coefficient; skipping this indicator"
        continue
    end

    β  = row.Coef[1]
    se = row.StdError[1]

    # 95% confidence interval (normal approximation)
    ci_low  = β - 1.96 * se
    ci_high = β + 1.96 * se

    # 4.2 Pull the LCIA time series corresponding to this indicator
    # Column names in lcia_results_per_cubic_meter_AnD_feed_in.csv
    # match model_df.Indicator, e.g., "Global Climate Change Potential - GCCP"
    @assert ind in names(lcia_df) "LCIA data does not contain column $ind"
    y = lcia_df[!, ind]

    # 4.3 Drop possible missings to keep x / y aligned
    mask = .!(ismissing.(x_sso) .| ismissing.(y))
    x_clean = x_sso[mask]
    y_clean = y[mask]

    μx = mean(x_clean)
    μy = mean(y_clean)
    σx = std(x_clean)
    σy = std(y_clean)

    # Standardized coefficient: slope after scaling x and y
    β_std = β * (σx / σy)

    # Elasticity: β * (mean(x) / mean(y))
    elasticity = β * (μx / μy)

    # Absolute and relative (% ) change when SSO increases by 10%
    Δx_10 = 0.10 * μx
    Δy_10 = β * Δx_10
    rel_Δy_pct = 100 * Δy_10 / μy

    push!(results, (
        String(ind),
        β,
        se,
        ci_low,
        ci_high,
        β_std,
        elasticity,
        Δy_10,
        rel_Δy_pct
    ))
end

# ---------------------------
# 5. Write results
# ---------------------------
CSV.write(out_file, results)
println("SSO sensitivity results written to: $out_file")
