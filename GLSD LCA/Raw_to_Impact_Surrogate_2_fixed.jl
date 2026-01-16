import Pkg
Pkg.activate(@__DIR__)

function ensure_deps(deps::Vector{String})
    missing = [d for d in deps if isnothing(Base.find_package(d))]
    if !isempty(missing)
        @info "Instantiating environment to install missing deps" missing
        Pkg.instantiate()
    end
end

ensure_deps(["CSV", "DataFrames", "LinearAlgebra", "Statistics", "Random", "Plots", "StatsPlots", "Measures"])

using CSV, DataFrames, LinearAlgebra, Statistics, Random
using Plots, StatsPlots, Measures

# ------------------------------------------------------------
# Raw data → LCI predictors → LCIA surrogate (no retraining)
# ------------------------------------------------------------
# This script:
#   1) loads the per-FU predictor table (constructed from SCADA/engineering logic),
#   2) applies the fitted multi-response linear coefficients,
#   3) exports predictions, coefficient matrices, and residual diagnostics plots.

gr()
default(fontfamily = "Arial")

# -------------------------
# Paths
# -------------------------
INPUT_DATA_PATH         = "cleaned_data_per_cubic_meter_AnD_feed_in.csv"
MODEL_COEF_PATH         = "model_full_summary_per_cubic_meter_AnD_feed_in.csv"
LCIA_RESULTS_PATH       = "lcia_results_per_cubic_meter_AnD_feed_in.csv"
OUTPUT_PRED_PATH        = "prediction_raw_to_impact_surrogate.csv"
OUTPUT_COEF_MATRIX_PATH = "raw_over_feed_coeff_matrix.csv"
PANEL_PATH              = "panel_FU_digester_feed.csv"  # created by RSM_FU_digester_feed.jl

# Plot output directory used by the paper
PLOT_DIR = "plots_per_AD_feed_in"
mkpath(PLOT_DIR)

# -------------------------
# Include predictor construction logic (used upstream in the workflow)
# -------------------------
include("LCI_Calculation_FU_PerCubicMeter_AD_feedin V2.jl")

# -------------------------
# Load predictor table
# -------------------------
function ensure_time_column!(df::DataFrame)
    for n in propertynames(df)
        if lowercase(string(n)) == "time"
            n == :time || rename!(df, n => :time)
            return true
        end
    end
    return false
end

df_cleaned = CSV.read(INPUT_DATA_PATH, DataFrame)
rename!(df_cleaned, collect(propertynames(df_cleaned)) .|> Symbol)
@assert ensure_time_column!(df_cleaned) "Expected a :time column in $(INPUT_DATA_PATH)."
time_col = df_cleaned.time

# Keep all numeric predictors (exclude time + any bookkeeping index)
X_lci_full = df_cleaned[:, Not([:time, :index])]

# -------------------------
# Load fitted coefficients (long format)
# -------------------------
coef_df = CSV.read(MODEL_COEF_PATH, DataFrame)
rename!(coef_df, collect(propertynames(coef_df)) .|> Symbol)

# Identify LCIA categories and predictor names
impacts    = unique(String.(coef_df.Indicator))
predictors = unique(String.(filter(!=("Intercept"), coef_df.Variable)))
predictor_syms = Symbol.(predictors)

# Construct coefficient matrices
β0 = Dict(imp => coef_df[(coef_df.Indicator .== imp) .& (coef_df.Variable .== "Intercept"), :Coef][1]
         for imp in impacts)

β = Dict((pred, imp) => coef_df[(coef_df.Indicator .== imp) .& (coef_df.Variable .== pred), :Coef][1]
         for pred in predictors, imp in impacts)

# Subset and order X
X_dedup = X_lci_full[:, predictor_syms]

# -------------------------
# Predict LCIA indicators
# -------------------------
Y_pred = DataFrame(time = time_col)
for imp in impacts
    yhat = fill(β0[imp], nrow(X_dedup))
    for (pred_sym, pred_str) in zip(predictor_syms, predictors)
        yhat .+= β[(pred_str, imp)] .* X_dedup[!, pred_sym]
    end
    Y_pred[!, Symbol(imp)] = yhat
end

CSV.write(OUTPUT_PRED_PATH, Y_pred)
println("\n✅ Saved predictions: $(OUTPUT_PRED_PATH)")

# -------------------------
# Export coefficient matrix used by the controller-ready implementation
# (rows = predictors, cols = impacts)
# -------------------------
coef_mat = Array{Float64}(undef, length(predictors), length(impacts))
for (i, pred) in enumerate(predictors)
    for (j, imp) in enumerate(impacts)
        coef_mat[i, j] = β[(pred, imp)]
    end
end

coef_mat_df = DataFrame(coef_mat, Symbol.(impacts))
coef_mat_df[!, :Predictor] = predictors
select!(coef_mat_df, :Predictor, Symbol.(impacts)...)  # Predictor first
CSV.write(OUTPUT_COEF_MATRIX_PATH, coef_mat_df)
println("✅ Saved coefficient matrix: $(OUTPUT_COEF_MATRIX_PATH)")

# -------------------------
# Residual diagnostics (requires OpenLCA outputs)
# -------------------------
if isfile(LCIA_RESULTS_PATH)
    lcia_df = CSV.read(LCIA_RESULTS_PATH, DataFrame)
    rename!(lcia_df, collect(propertynames(lcia_df)) .|> Symbol)
    @assert ensure_time_column!(lcia_df) "Expected a :time column in $(LCIA_RESULTS_PATH)."

    # Determine common impact columns in LCIA table
    impacts_in_lcia = [c for c in impacts if Symbol(c) in propertynames(lcia_df)]
    @assert !isempty(impacts_in_lcia) "No matching LCIA indicator columns found in $(LCIA_RESULTS_PATH)."

    # Inner join on :time to guarantee day-level alignment between predictors and OpenLCA outputs
    lcia_cols = vcat([:time], Symbol.(impacts_in_lcia))
    lcia_keep = select(lcia_df, lcia_cols)
    df_join = innerjoin(df_cleaned, lcia_keep, on=:time, makeunique=true)
    df_join = sort(df_join, :time)
    df_join = dropmissing(df_join, Symbol.(impacts_in_lcia))

    # Rebuild prediction matrix in the joined order
    X_sorted = df_join[:, predictor_syms]
    times_sorted = df_join.time

    Y_pred_sorted = Matrix{Float64}(undef, nrow(X_sorted), length(impacts_in_lcia))
    for (j, imp) in enumerate(impacts_in_lcia)
        yhat = fill(β0[imp], nrow(X_sorted))
        for (pred_sym, pred_str) in zip(predictor_syms, predictors)
            yhat .+= β[(pred_str, imp)] .* X_sorted[!, pred_sym]
        end
        Y_pred_sorted[:, j] = yhat
    end

    Y_true_sorted = Matrix(df_join[:, Symbol.(impacts_in_lcia)])

    # -------------------------
    # Held-out test set (prefer reading the partition labels saved by the inventory-based run)
    # -------------------------
    if isfile(PANEL_PATH)
        panel = CSV.read(PANEL_PATH, DataFrame)
        rename!(panel, collect(propertynames(panel)) .|> Symbol)
        if !ensure_time_column!(panel)
            if :date in propertynames(panel)
                rename!(panel, :date => :time)
            else
                error("$(PANEL_PATH) must contain a :time (or :date) column.")
            end
        end
        @assert :partition in propertynames(panel) "Expected a :partition column in $(PANEL_PATH)."
        test_set = Set(string.(panel.time[panel.partition .== "test"]))
        test_idx = findall(in.(string.(times_sorted), Ref(test_set)))
        @assert !isempty(test_idx) "No overlapping test days between $(PANEL_PATH) and $(LCIA_RESULTS_PATH)."
    else
        @warn "($(PANEL_PATH)) not found; falling back to an unstratified 80/10/10 split with seed=2025."
        N = nrow(df_join)
        idx = collect(1:N)
        Random.seed!(2025)
        shuffle!(idx)
        n_train = round(Int, 0.80*N)
        n_val   = round(Int, 0.10*N)
        test_idx = idx[n_train+n_val+1:end]
    end

    Y_test_true = Y_true_sorted[test_idx, :]
    Y_test_pred = Y_pred_sorted[test_idx, :]

    # -------------------------
    # Residuals: define residual = surrogate - OpenLCA
    # -------------------------
    residuals_raw = Dict{String, Vector{Float64}}()
    residuals_rel = Dict{String, Vector{Float64}}()  # signed %, normalized by |truth|

    for (j, imp) in enumerate(impacts_in_lcia)
        yv = Y_test_true[:, j]
        pv = Y_test_pred[:, j]
        residual_vec = pv .- yv

        # robust normalization to avoid division blow-ups near zero
        eps = 1e-6 * max(1.0, median(abs.(yv)))
        denom = clamp.(abs.(yv), eps, Inf)

        residuals_raw[imp] = residual_vec
        residuals_rel[imp] = 100.0 .* residual_vec ./ denom
    end

    # -------------------------
    # Plot helpers
    # -------------------------
    short_label(s::AbstractString) = begin
        if occursin("Water Use", s)
            return "WU"
        end
        m = match(r"-\s*([A-Za-z0-9]+)\s*$", s)
        return m === nothing ? s : m.captures[1]
    end

    indicator_order = impacts_in_lcia
    ylabels_abbr = short_label.(indicator_order)

    function heatplot(mat_dict::Dict{String, Vector{Float64}};
                      title::String,
                      file::String,
                      ylabel::String)
        test_len = length(first(values(mat_dict)))
        mat = reduce(hcat, [mat_dict[imp] for imp in indicator_order])
        matT = permutedims(mat)  # rows = categories

        # symmetric color limits centered at zero
        maxabs = maximum(abs, matT)
        maxabs = max(maxabs, eps())

        xtick_pos = range(1, test_len, length=10) |> collect
        xtick_pos = unique(round.(Int, xtick_pos))

        p = heatmap(
            1:test_len, 1:length(indicator_order), matT;
            yticks = (1:length(indicator_order), ylabels_abbr),
            xticks = (xtick_pos, string.(xtick_pos)),
            xlabel = "Test Sample Index",
            ylabel = ylabel,
            title  = title,
            size   = (1200, 600),
            left_margin = 10mm,
            bottom_margin = 10mm,
            top_margin = 10mm,
            color = :coolwarm,
            clims = (-maxabs, maxabs),
            colorbar = true,
            legend = false,
            xgrid = false,
            ygrid = false,
            grid = false,
            tickfont = font("Arial", 10),
            guidefont = font("Arial", 14),
            titlefont = font("Arial", 14)
        )

        # Draw only an outer border (avoid heavy cell-by-cell grids)
        plot!(p; xlims = (0.5, test_len + 0.5), ylims = (0.5, length(indicator_order) + 0.5))
        vline!(p, [0.5, test_len + 0.5]; color = :black, lw = 1.0)
        hline!(p, [0.5, length(indicator_order) + 0.5]; color = :black, lw = 1.0)

        savefig(p, joinpath(PLOT_DIR, file))
        println("✅ Saved: $(joinpath(PLOT_DIR, file))")
        return p
    end

    # Heatmaps
    heatplot(residuals_raw;
        title = "Residual Heatmap (Raw, Test Set): Surrogate - OpenLCA",
        file  = "residuals_raw_heatmap_raw_to_impact.svg",
        ylabel = "Impact Category")

    heatplot(residuals_rel;
        title = "Residual Heatmap (Relative %, Test Set): Surrogate - OpenLCA",
        file  = "residuals_rel_heatmap_raw_to_impact.svg",
        ylabel = "Impact Category")

    # -------------------------
    # Boxplots (test set)
    # -------------------------
    # Signed relative residuals (%): distribution by impact category
    rel_df = DataFrame(
        Category = repeat(ylabels_abbr, inner = length(residuals_rel[indicator_order[1]])),
        RelResidualPct = vcat([residuals_rel[imp] for imp in indicator_order]...)
    )

    p_box_rel = boxplot(
        rel_df.Category,
        rel_df.RelResidualPct;
        legend = false,
        xlabel = "Impact Category",
        ylabel = "Relative residual (%, surrogate - OpenLCA)",
        title  = "Relative residual distribution (test set)",
        xrotation = 45,
        size = (1200, 500),
        left_margin = 10mm,
        bottom_margin = 15mm,
        tickfont = font("Arial", 10),
        guidefont = font("Arial", 14),
        titlefont = font("Arial", 14)
    )
    hline!(p_box_rel, [0.0]; color = :black, lw = 1.0)
    savefig(p_box_rel, joinpath(PLOT_DIR, "boxplot_rel_residual_signed_raw_to_impact.svg"))
    println("✅ Saved: $(joinpath(PLOT_DIR, "boxplot_rel_residual_signed_raw_to_impact.svg"))")

    # Absolute relative error (%)
    abs_rel_df = DataFrame(
        Category = rel_df.Category,
        AbsRelErrorPct = abs.(rel_df.RelResidualPct)
    )

    p_box_abs = boxplot(
        abs_rel_df.Category,
        abs_rel_df.AbsRelErrorPct;
        legend = false,
        xlabel = "Impact Category",
        ylabel = "Absolute relative error (%)",
        title  = "Absolute relative error distribution (test set)",
        xrotation = 45,
        size = (1200, 500),
        left_margin = 10mm,
        bottom_margin = 15mm,
        tickfont = font("Arial", 10),
        guidefont = font("Arial", 14),
        titlefont = font("Arial", 14)
    )
    savefig(p_box_abs, joinpath(PLOT_DIR, "boxplot_rel_error_absolute_raw_to_impact.svg"))
    println("✅ Saved: $(joinpath(PLOT_DIR, "boxplot_rel_error_absolute_raw_to_impact.svg"))")

else
    println("\n⚠️ Skipped residual diagnostics: $(LCIA_RESULTS_PATH) not found.")
end
