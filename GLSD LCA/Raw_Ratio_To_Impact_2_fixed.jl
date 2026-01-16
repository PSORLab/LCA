import Pkg
Pkg.activate(@__DIR__)

function ensure_deps(deps::Vector{String})
    missing = [d for d in deps if isnothing(Base.find_package(d))]
    if !isempty(missing)
        @info "Instantiating environment to install missing deps" missing
        Pkg.instantiate()
    end
end

ensure_deps(["CSV", "DataFrames", "Random", "Statistics", "Plots", "StatsPlots", "Measures"])

using CSV
using DataFrames
using Random
using Statistics
using Plots, StatsPlots, Measures

# ------------------------------------------------------------------
# Predict environmental impacts from raw SCADA ratios (raw / digester feed)
# using the coefficient matrix exported by Raw_to_Impact_Surrogate.jl
#
# Inputs:
#   - raw_inputs_with_feed.csv
#   - raw_over_feed_coeff_matrix.csv
# Optional (for residual diagnostics):
#   - lcia_results_per_cubic_meter_AnD_feed_in.csv
# Outputs:
#   - prediction_raw_ratio_to_impact_surrogate.csv
#   - residual heatmaps + boxplots in plots_per_AD_feed_in/
# ------------------------------------------------------------------

default(fontfamily = "Arial")

# ----------------------
# File paths
# ----------------------
raw_path  = "raw_inputs_with_feed.csv"
coef_path = "raw_over_feed_coeff_matrix.csv"
out_pred  = "prediction_raw_ratio_to_impact_surrogate.csv"

# Used only if present
lcia_results_path = "lcia_results_per_cubic_meter_AnD_feed_in.csv"
PANEL_PATH = "panel_FU_digester_feed.csv"  # created by RSM_FU_digester_feed.jl
MODEL_SUMMARY_PATH = "model_full_summary_per_cubic_meter_AnD_feed_in.csv"  # for intercepts

# Plot output folder (matches your Overleaf folder naming)
output_dir = "plots_per_AD_feed_in"
mkpath(output_dir)

# ----------------------
# Load raw inputs
# ----------------------
raw_inputs = CSV.read(raw_path, DataFrame)
rename!(raw_inputs, Symbol.(names(raw_inputs)))
@info "Loaded raw inputs" nrow=nrow(raw_inputs) ncol=ncol(raw_inputs)

function ensure_time_column!(df::DataFrame)
    for n in names(df)
        if lowercase(strip(String(n))) == "time"
            if n != :time
                rename!(df, n => :time)
            end
            return true
        end
    end
    return false
end

# Ensure time exists (needed for alignment with LCIA results)
has_time = ensure_time_column!(raw_inputs)
if !has_time
    @warn "raw_inputs_with_feed.csv has no 'time' column; residual diagnostics will be skipped."
end

# ----------------------
# Load coefficient matrix
# ----------------------
coef_df = CSV.read(coef_path, DataFrame)
rename!(coef_df, Symbol.(names(coef_df)))
coef_col_syms = Symbol.(names(coef_df))
predictor_col = if :Predictor in coef_col_syms
    :Predictor
elseif :Variable in coef_col_syms
    :Variable
else
    error("Coefficient matrix must have a 'Predictor' or 'Variable' column. Found: $(names(coef_df))")
end
coef_df = select(coef_df, predictor_col, Not(predictor_col))

impact_cols = Symbol.(names(coef_df)[2:end])
@info "Impact categories" impact_cols

# Coef matrix: rows=predictors, cols=impacts
coef_mat = Matrix(coef_df[:, 2:end])
predictor_names = string.(coef_df[!, predictor_col])

# Convenience lookup
pred_idx = Dict(predictor_names[i] => i for i in eachindex(predictor_names))

# ----------------------
# Compute per-feed predictors (match LCI workflow)
# ----------------------
function compute_predictors(df::DataFrame)
    # Column names expected in raw_inputs_with_feed.csv
    function get_col(candidates::Vector{Symbol})
        for c in candidates
            if c in propertynames(df)
                return df[!, c]
            end
        end
        error("Missing required column. Expected one of: $(candidates)")
    end

    feed_kgal = get_col([:dig_tfeed])
    sso_mass  = get_col([:SSO_mass_kg, :SSO_mass])
    ae_raw    = get_col([:AE_raw_MJ, :AE_MJ])
    ah_raw    = get_col([:AH_raw_MJ, :AH_MJ])
    ng_raw    = get_col([:ng_energy_MJ, :ng_energy])
    steel_raw = get_col([:steel_raw_const, :steel])

    # Denominators consistent with LCI_Calculation_FU_PerCubicMeter_AD_feedin V2.jl
    flow_m3_per_day   = feed_kgal .* 3.78541
    flow_Mgal_per_day = feed_kgal .* 1e-3

    safe_div(num, den) = map(num, den) do n, d
        if ismissing(n) || ismissing(d) || d <= 0
            missing
        else
            n / d
        end
    end

    DataFrame(
        time   = hasproperty(df, :time) ? df[!, :time] : collect(1:nrow(df)),
        SSO_LCI= safe_div(sso_mass, flow_m3_per_day),
        AE     = safe_div(ae_raw,   flow_m3_per_day),
        AH     = safe_div(ah_raw,   flow_m3_per_day),
        ng_LCI = safe_div(ng_raw,   flow_m3_per_day),
        steel  = safe_div(steel_raw, flow_Mgal_per_day)
    )
end

predictors_df = compute_predictors(raw_inputs)

# Map predictor columns → model predictor names
# (these predictor names must match the first column of raw_over_feed_coeff_matrix.csv)
predictor_syms = Symbol.(predictor_names)
missing_predictors = setdiff(predictor_syms, propertynames(predictors_df))
if !isempty(missing_predictors)
    error("Missing predictors in input: $(string.(missing_predictors))")
end

# Design matrix aligned to coef row order
Xmat = Matrix(predictors_df[:, predictor_syms])

# ----------------------
# Intercepts (from model summary, if available)
# ----------------------
intercepts = zeros(Float64, length(impact_cols))
if isfile(MODEL_SUMMARY_PATH)
    model_df = CSV.read(MODEL_SUMMARY_PATH, DataFrame)
    rename!(model_df, Symbol.(names(model_df)))
    if all(x -> x in propertynames(model_df), [:Indicator, :Variable, :Coef])
        intercept_map = Dict(String(row.Indicator) => Float64(row.Coef)
                             for row in eachrow(model_df) if row.Variable == "Intercept")
        missing_impacts = setdiff(String.(impact_cols), collect(keys(intercept_map)))
        if !isempty(missing_impacts)
            @warn "Missing intercepts for some impacts; defaulting to 0" missing_impacts
        end
        intercepts = [get(intercept_map, String(c), 0.0) for c in impact_cols]
    else
        @warn "Model summary missing required columns for intercept extraction" cols=names(model_df)
    end
else
    @warn "Model summary not found; assuming zero intercepts" MODEL_SUMMARY_PATH
end

# ----------------------
# Predict impacts
# ----------------------
Yhat = Xmat * coef_mat
Yhat .+= intercepts'

# Save predictions
pred_df = DataFrame(time = predictors_df.time)
for (j, c) in enumerate(impact_cols)
    pred_df[!, c] = Yhat[:, j]
end
CSV.write(out_pred, pred_df)
@info "Saved predictions" out_pred

# ----------------------
# Helper: abbreviate impact names
# ----------------------
function short_label(name::AbstractString)
    # Extract trailing abbreviation after " - ", if present
    if occursin(" - ", name)
        return strip(split(name, " - ")[end])
    elseif occursin("Water Use", name)
        return "WU"
    else
        return name
    end
end

# ----------------------
# Residual diagnostics (if LCIA results exist)
# ----------------------
if has_time && isfile(lcia_results_path)
    lcia_df = CSV.read(lcia_results_path, DataFrame)
    rename!(lcia_df, Symbol.(names(lcia_df)))
    if !ensure_time_column!(lcia_df)
        error("$lcia_results_path must contain a 'time' column")
    end

    # Keep only impacts that exist in BOTH predicted and LCIA tables
    pred_impacts = setdiff(propertynames(pred_df), [:time])
    lcia_impacts = setdiff(propertynames(lcia_df), [:time])
    common_impacts = sort(intersect(pred_impacts, lcia_impacts))

    if isempty(common_impacts)
        error("No common LCIA indicator columns found between predictions and LCIA results.")
    end

    # Rename LCIA columns to *_y to avoid collisions
    lcia_sub = select(lcia_df, [:time; common_impacts])
    rename!(lcia_sub, Dict(c => Symbol("$(c)_y") for c in common_impacts))

    # Merge predictions with LCIA truth
    df_merged = leftjoin(pred_df, lcia_sub, on = :time)
    df_merged = dropmissing(df_merged)

    @info "Merged prediction + LCIA tables" nrow=nrow(df_merged)

    # ----------------------
    # Held-out test set (prefer reading the partition labels saved by the inventory-based run)
    # ----------------------
    test_set = nothing
    if isfile(PANEL_PATH)
        panel = CSV.read(PANEL_PATH, DataFrame)
        if !(:time in propertynames(panel))
            if :date in propertynames(panel)
                rename!(panel, :date => :time)
            else
                error("Panel file $(PANEL_PATH) must contain a 'time' (or 'date') column.")
            end
        end
        if !(:partition in propertynames(panel))
            error("Panel file $(PANEL_PATH) must contain a 'partition' column.")
        end
        test_set = Set(string.(panel.time[panel.partition .== "test"]))
    else
        @warn "($(PANEL_PATH)) not found; falling back to an unstratified 80/10/10 split with seed=2025. For paper figures, generate panel_FU_digester_feed.csv by running RSM_FU_digester_feed.jl."
        n_total = nrow(raw_inputs)
        idx = shuffle(MersenneTwister(2025), 1:n_total)
        n_train = floor(Int, 0.8 * n_total)
        n_val   = floor(Int, 0.1 * n_total)
        test_idx = idx[(n_train + n_val + 1):end]
        test_times = raw_inputs[test_idx, :time]
        test_set = Set(string.(test_times))
    end

    df_test = df_merged[in.(string.(df_merged.time), Ref(test_set)), :]
    df_test = sort(df_test, :time)

    if nrow(df_test) == 0
        error("After aligning by time, the test-set intersection is empty. Check that raw_inputs_with_feed.csv and lcia_results_per_cubic_meter_AnD_feed_in.csv use the same time keys.")
    end

    @info "Test-set residual diagnostics" nrow=nrow(df_test)

    # ----------------------
    # Compute residuals on the test set
    # ----------------------
    residuals_raw = Dict{String, Vector{Float64}}()
    residuals_rel = Dict{String, Vector{Float64}}()   # signed, percent
    abs_rel_error = Dict{String, Vector{Float64}}()   # unsigned, percent

    for c_sym in common_impacts
        c = String(c_sym)
        yhat = Vector{Float64}(df_test[!, c_sym])
        ytrue = Vector{Float64}(df_test[!, Symbol("$(c)_y")])

        r = yhat .- ytrue

        # Robust denominator: avoid exploding % errors when the true value is very small
        eps = 1e-6 * max(1.0, median(abs.(ytrue)))
        den = clamp.(abs.(ytrue), eps, Inf)

        rel = 100 .* r ./ den

        residuals_raw[c] = r
        residuals_rel[c] = rel
        abs_rel_error[c] = abs.(rel)
    end

    # ----------------------
    # Performance summary (test)
    # ----------------------
    perf = DataFrame(Indicator = String[], R2_test = Float64[], RMSE_test = Float64[])
    for c_sym in common_impacts
        c = String(c_sym)
        yhat = Vector{Float64}(df_test[!, c_sym])
        ytrue = Vector{Float64}(df_test[!, Symbol("$(c)_y")])

        r = yhat .- ytrue
        rmse = sqrt(mean(r .^ 2))

        sst = sum((ytrue .- mean(ytrue)).^2)
        r2 = sst > 0 ? 1.0 - sum(r.^2) / sst : NaN

        push!(perf, (c, r2, rmse))
    end

    @info "Test performance summary (raw-ratio surrogate)" perf

    # ----------------------
    # Heatmaps (raw residuals; relative residuals)
    # ----------------------
    indicator_order = [String(c) for c in common_impacts]
    ylabels_abbr = short_label.(indicator_order)

    function heatplot(dict::Dict{String, Vector{Float64}}; title::String, cbar::String)
        mat = reduce(hcat, [dict[k] for k in indicator_order])
        matT = permutedims(mat)
        maxabs = maximum(abs, matT)
        
        ntest = size(matT, 2)
        xtick_pos = range(1, ntest; length = min(10, ntest))
        xtick_pos = unique(round.(Int, collect(xtick_pos)))

        heatmap(
            1:ntest, 1:length(indicator_order), matT;
            xticks = (xtick_pos, string.(xtick_pos)),
            yticks = (1:length(indicator_order), ylabels_abbr),
            xlabel = "Test sample index",
            ylabel = "Impact category",
            title = title,
            color = :coolwarm,
            clims = (-maxabs, maxabs),
            size = (1400, 650),
            left_margin = 10mm,
            bottom_margin = 10mm,
            top_margin = 10mm,
            right_margin = 10mm,
            xguidefont = font("Arial", 14),
            yguidefont = font("Arial", 14),
            tickfont = font("Arial", 10),
            titlefont = font("Arial", 14),
            legend = false,
            colorbar_title = cbar,
            grid = false
        )
    end

    p_raw = heatplot(residuals_raw;
        title = "Residual heatmap (raw, signed) — raw-ratio surrogate (test set)",
        cbar  = "Residual"
    )
    savefig(p_raw, joinpath(output_dir, "residuals_raw_heatmap_raw_ratio.svg"))

    p_rel = heatplot(residuals_rel;
        title = "Residual heatmap (relative, signed) — raw-ratio surrogate (test set)",
        cbar  = "Residual (%)"
    )
    savefig(p_rel, joinpath(output_dir, "residuals_rel_heatmap_raw_ratio.svg"))

    # ----------------------
    # Boxplots
    # ----------------------
    df_rel = DataFrame(
        Indicator = repeat(ylabels_abbr, inner = nrow(df_test)),
        RelResidualPct = vcat([residuals_rel[k] for k in indicator_order]...),
        AbsRelErrorPct = vcat([abs_rel_error[k] for k in indicator_order]...)
    )

    bp_signed = @df df_rel boxplot(
        :Indicator, :RelResidualPct;
        xlabel = "Impact category",
        ylabel = "Relative residual (%)",
        title  = "Relative residuals (signed) — raw-ratio surrogate (test set)",
        legend = false,
        xrotation = 45,
        size = (1400, 600),
        left_margin = 10mm,
        bottom_margin = 15mm,
        top_margin = 10mm,
        right_margin = 10mm,
        outliers = true
    )
    hline!(bp_signed, [0.0]; color=:black, lw=1.0, alpha=0.6)
    savefig(bp_signed, joinpath(output_dir, "boxplot_rel_residual_signed_raw_ratio.svg"))

    bp_abs = @df df_rel boxplot(
        :Indicator, :AbsRelErrorPct;
        xlabel = "Impact category",
        ylabel = "Absolute relative error (%)",
        title  = "Absolute relative errors — raw-ratio surrogate (test set)",
        legend = false,
        xrotation = 45,
        size = (1400, 600),
        left_margin = 10mm,
        bottom_margin = 15mm,
        top_margin = 10mm,
        right_margin = 10mm,
        outliers = true
    )
    savefig(bp_abs, joinpath(output_dir, "boxplot_abs_rel_error_raw_ratio.svg"))

    @info "Saved plots" output_dir
else
    @warn "LCIA results file not found — skipping residual diagnostics" lcia_results_path
end
