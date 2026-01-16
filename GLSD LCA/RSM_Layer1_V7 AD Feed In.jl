using CSV, DataFrames, GLM, StatsModels
using Plots
using Statistics
using StatsBase
using Random
using Distributions
using Measures
using LinearAlgebra
gr()  # set GR as the plotting backend

Random.seed!(42)  # fix the split

# Step 1: Read and preprocess predictors
X_all = CSV.read("cleaned_data_per_cubic_meter_AnD_feed_in.csv", DataFrame)
X = select(X_all, Not(["index", "time", "A_SSO_WTE", "A_SSO_LF", "building", "concrete", "earth", "gravel","elec"]))  # keep the actual predictors

function compute_vif(X::DataFrame)
    vif_table = DataFrame(Variable=String[], VIF=Float64[])
    for col in names(X)
        y = X[!, col]
        others = select(X, Not(col))
        rhs_terms = reduce(+, Term.(Symbol.(names(others))))
        formula = Term(Symbol(col)) ~ rhs_terms
        model = lm(formula, hcat(others, DataFrame(Symbol(col) => y)))
        r2 = GLM.r2(model)
        vif = 1 / max(1e-12, (1 - r2))
        push!(vif_table, (string(col), vif))
    end
    return vif_table
end

println("🔍 VIF collinearity check results:")
vif_table = compute_vif(X)
display(vif_table)
CSV.write("vif_result_per_cubic_meter_AnD_feed_in.csv", vif_table)

# Set threshold
threshold = 0.95
# Step 1: Handle highly correlated variables
function remove_highly_correlated(X::DataFrame, threshold::Float64)
    cor_matrix = cor(Matrix(X))
    drop_set = Set{String}()
    var_names = names(X) .|> string  # variable names (strings)

    for i in 1:length(var_names)-1
        for j in i+1:length(var_names)
            corr_val = cor_matrix[i, j]
            if abs(corr_val) ≥ threshold
                println("⚠️ Highly correlated pair: $(var_names[i]) & $(var_names[j]) → corr = $(round(corr_val, digits=4))")
                push!(drop_set, var_names[j])  # keep the first, drop the second
            end
        end
    end
    return Symbol.(setdiff(var_names, drop_set))
end

kept_vars = remove_highly_correlated(X, threshold)
println("\n✅ Variables kept:")
display(kept_vars)

# Export X after correlation pruning
X_dedup = select(X, kept_vars)
CSV.write("cleaned_data_dedup_per_cubic_meter_AnD_feed_in.csv", X_dedup)

# Step 2: Read and preprocess dependent variables
Y_all = CSV.read("lcia_results_per_cubic_meter_AnD_feed_in.csv", DataFrame)

# Remove rows with missing dependent variables and sync features/time
mask = completecases(Y_all)
missing_cnt = count(!, mask)
if missing_cnt > 0
    println("⚠️ Detected $(missing_cnt) rows with missing LCIA entries; removed corresponding samples.")
end
Y_all = Y_all[mask, :]
X_dedup = X_dedup[mask, :]

Y = select(Y_all, Not(["index", "time"]))  # LCIA indicators
time_col = Y_all.time                      # time column for output

# Step 3: Prepare output containers
full_summary_table = DataFrame(
    Indicator = String[],
    Variable = String[],
    Coef = Float64[],
    StdError = Float64[],
    t = Float64[],
    P = Float64[],
    R2 = Float64[],
    RMSE = Float64[]
)

summary_table = DataFrame(Indicator = String[], R2 = Float64[], RMSE = Float64[])
mkpath("plots_per_cubic_meter_AnD_feed_in")

# —— Aggregated residual plot data (store test only) ——
combined_abs = Dict{String, Vector{Float64}}()
combined_rel = Dict{String, Vector{Float64}}()

# Step 4: Multi-input, multi-output modeling (strict 80/10/10)
println("\n==============================")
println("Multi-input multi-output linear surrogate (fit all indicators at once)")

# 80/10/10 split (shared across all indicators)
N = nrow(X_dedup)
idx = shuffle(1:N)
ntr = floor(Int, 0.80N)
nva = floor(Int, 0.10N)
train_idx = idx[1:ntr]
val_idx   = idx[ntr+1 : ntr+nva] 
test_idx  = idx[ntr+nva+1 : end]

# Matrix-ify and add an intercept column
function with_intercept(mat::AbstractMatrix)
    return hcat(ones(size(mat, 1)), mat)
end

X_train = with_intercept(Matrix(X_dedup[train_idx, :]))
X_val   = with_intercept(Matrix(X_dedup[val_idx, :]))
X_test  = with_intercept(Matrix(X_dedup[test_idx, :]))
X_trainval = with_intercept(Matrix(vcat(X_dedup[train_idx, :], X_dedup[val_idx, :])))

Y_train = Matrix(Y[train_idx, :])
Y_val   = Matrix(Y[val_idx, :])
Y_test  = Matrix(Y[test_idx, :])
Y_trainval = Matrix(vcat(Y[train_idx, :], Y[val_idx, :]))

# Fit multi-output linear regression (least squares)
function fit_multioutput(X::AbstractMatrix, Y::AbstractMatrix)
    return X \ Y  # shape: (p, k)
end

# Train (train used for val; train+val used for final)
beta_train = fit_multioutput(X_train, Y_train)
Y_val_pred = X_val * beta_train

# Validation evaluation (per indicator)
for (j, y_name) in enumerate(names(Y))
    y_val_vec = Y_val[:, j]
    y_pred_vec = Y_val_pred[:, j]
    r2_val = 1 - sum((y_val_vec .- y_pred_vec).^2) / sum((y_val_vec .- mean(y_val_vec)).^2)
    rmse_val = sqrt(mean((y_val_vec .- y_pred_vec).^2))
    println("(Val) $y_name → R² = $(round(r2_val, digits=6)), RMSE = $(round(rmse_val, digits=6))")
end

# Final model (train+val)
beta_final = fit_multioutput(X_trainval, Y_trainval)

# Test evaluation
Y_test_pred = X_test * beta_final

# ============================================================
# Residual covariance/correlation across LCIA categories (TEST)
# Fixes A9: ensures we use residuals E (not Y) and prevents GCCP/SFP duplication
# ============================================================

using LinearAlgebra
using DataFrames
using CSV

# E_test: signed residual matrix on the held-out test set (N_test × C)
E_test = Y_test .- Y_test_pred

N_test = size(E_test, 1)
D = size(X_test, 2)          # X_test already includes intercept because you used with_intercept(...)
dof_test = N_test - D
@assert dof_test > 0 "DoF for Sigma_E is non-positive. Check N_test and the number of regressors."

# Residual covariance (Sigma_E) and correlation (R_E)
SigmaE = (E_test' * E_test) / dof_test
sd = sqrt.(diag(SigmaE))
R_E = SigmaE ./ (sd * sd')

# ---- Build a short-name label set (optional but matches your SI tables) ----
impact_full = collect(names(Y))

label_map = Dict(
    "Eutrophication Potential - EP" => "EP",
    "Cumulative Energy Demand - CED" => "CED",
    "Particulate Matter Formation Potential - PMFP" => "PMFP",
    "Fossil Depletion Potential - FDP" => "FDP",
    "Water Use" => "WU",
    "Global Climate Change Potential - GCCP" => "GCCP",
    "Acidification Potential - AP" => "AP",
    "Smog Formation Potential - SFP" => "SFP",
)

impact_short = [get(label_map, s, s) for s in impact_full]

# ---- A9 hard guard: GCCP and SFP must be distinct columns ----
i_gccp = findfirst(==("GCCP"), impact_short)
i_sfp  = findfirst(==("SFP"),  impact_short)

@assert i_gccp !== nothing && i_sfp !== nothing "Cannot find GCCP and/or SFP in impact_short. Check label_map / LCIA column names."
@assert i_gccp != i_sfp "GCCP and SFP resolve to the same column index. This indicates a label/indexing bug."
@assert !all(E_test[:, i_gccp] .== E_test[:, i_sfp]) "A9 triggered: GCCP and SFP residual vectors are identical. You are likely reading GCCP twice."

# ---- Export to CSV so you can paste/update SI tables safely ----
SigmaE_df = DataFrame(SigmaE, Symbol.(impact_short))
insertcols!(SigmaE_df, 1, :Indicator => impact_short)

R_E_df = DataFrame(R_E, Symbol.(impact_short))
insertcols!(R_E_df, 1, :Indicator => impact_short)

CSV.write("SigmaE_cov_test.csv", SigmaE_df)
CSV.write("SigmaE_cor_test.csv", R_E_df)

println("✅ Wrote residual covariance/correlation tables:")
println("   - SigmaE_cov_test.csv")
println("   - SigmaE_cor_test.csv")

r2_dict = Dict{String, Float64}()
rmse_dict = Dict{String, Float64}()

# Evaluation and residual storage
for (j, y_name) in enumerate(names(Y))
    y_test_vec = Y_test[:, j]
    y_pred_vec = Y_test_pred[:, j]
    residuals = y_test_vec .- y_pred_vec
    r2_test = 1 - sum((y_test_vec .- y_pred_vec).^2) / sum((y_test_vec .- mean(y_test_vec)).^2)
    rmse_test = sqrt(mean(residuals .^ 2))
    push!(summary_table, (y_name, r2_test, rmse_test))
    r2_dict[y_name] = r2_test
    rmse_dict[y_name] = rmse_test

    # —— Aggregate plot data (store test only) ——
    eps = 1e-6 * max(1.0, median(abs.(y_test_vec)))
    rel_resid = 100 .* residuals ./ clamp.(abs.(y_test_vec), eps, Inf)
    combined_abs[y_name] = residuals
    combined_rel[y_name] = rel_resid
end

# Coefficients and significance (per indicator output)
var_names = ["Intercept"; string.(names(X_dedup))]
XtX_inv = inv(X_trainval' * X_trainval)
dof = size(X_trainval, 1) - size(X_trainval, 2)
tdist = TDist(dof)

for (j, y_name) in enumerate(names(Y))
    y_pred_trainval = X_trainval * beta_final[:, j]
    resid = Y_trainval[:, j] .- y_pred_trainval
    sigma2 = sum(resid .^ 2) / dof
    stderr_vals = sqrt.(diag(XtX_inv) .* sigma2)
    coef_vals = beta_final[:, j]
    tvals = coef_vals ./ stderr_vals
    pvals = [2 * (1 - cdf(tdist, abs(t))) for t in tvals]

    r2_test = r2_dict[y_name]
    rmse_test = rmse_dict[y_name]

    println("\nModel coefficients (final_model, multi-output): $y_name")
    show(IOContext(stdout, :compact => true, :limit => false),
         DataFrame(Variable = var_names,
                   Coef = coef_vals,
                   StdError = stderr_vals,
                   t = tvals,
                   P = pvals))
    println()

    for i in 1:length(coef_vals)
        push!(full_summary_table, (
            y_name,
            var_names[i],
            coef_vals[i],
            stderr_vals[i],
            tvals[i],
            pvals[i],
            r2_test,
            rmse_test
        ))
    end
    println("(Test) $y_name → R² = $(round(r2_test, digits=6)), RMSE = $(round(rmse_test, digits=6))")
end

# Save test predictions (single file, multi-output)
pred_df = DataFrame(time = time_col[test_idx])
for (j, y_name) in enumerate(names(Y))
    pred_df[!, Symbol("prediction_" * replace(y_name, r"[ \-]+" => "_"))] = Y_test_pred[:, j]
end
CSV.write("prediction_per_cubic_meter_AnD_feed_in_multi.csv", pred_df)

# Save final model summary table (test results)
CSV.write("model_summary_per_cubic_meter_AnD_feed_in.csv", summary_table)
println("\n✅ All models finished (80/10/10; final metrics based on 10% test). Outputs:")
println("📁 Multi-indicator predictions (test): prediction_per_cubic_meter_AnD_feed_in_multi.csv")
println("🖼️  Fit and residual plots (single indicator, test): plots_per_cubic_meter_AnD_feed_in/*.png")
println("📊 Model evaluation summary (test): model_summary_per_cubic_meter_AnD_feed_in.csv")
println("\n📋 Full coefficient summary (final_model on train+val):")
@show full_summary_table
CSV.write("model_full_summary_per_cubic_meter_AnD_feed_in.csv", full_summary_table)
println("📘 Full coefficient summary saved: model_full_summary_per_cubic_meter_AnD_feed_in.csv")

# —— Combined residual overview (test set only) ——
using ColorSchemes
using Plots
gr()

colors = distinguishable_colors(length(names(Y)))

println("\n📊 Drawing combined residual plots (test 10%)...")

# —— Initialize subplots: legend on top plot, none on bottom ——
p_abs = plot(title = "Absolute Residuals (Test)",
             xlabel = "Index", ylabel = "Residual",
             legend = :bottomright, grid = false, titlefont = font(14))

p_rel = plot(title = "Relative Residuals (Test)",
             xlabel = "Index", ylabel = "Relative Residual (%)",
             legend = false, grid = false, titlefont = font(14))   # keep lower plot legend disabled

for (i, y_name) in enumerate(names(Y))
    res  = get(combined_abs, y_name, Float64[])
    rrel = get(combined_rel, y_name, Float64[])
    if !isempty(res)
        scatter!(p_abs, 1:length(res), res;
                 label = y_name, markersize = 3, color = colors[i], markerstrokewidth = 0)
    end
    if !isempty(rrel)
        scatter!(p_rel, 1:length(rrel), rrel;
                 label = nothing, markersize = 3, color = colors[i], markerstrokewidth = 0)  # no labels on the lower plot
    end
end

# Force once more to prevent later overrides
plot!(p_abs; legend = :bottomleft)
plot!(p_rel; legend = false)

final = plot(p_abs, p_rel,
             layout = @layout([a; b]),
             size = (1000, 800),
             guidefont = font(12), tickfont = font(10))

savefig(final, "plots_per_cubic_meter_AnD_feed_in/residuals_combined.svg")
println("✅ Combined plot saved to plots_per_cubic_meter_AnD_feed_in/residuals_combined.svg")

# —— Heatmaps (test residuals, absolute and relative) ——
indicator_order = collect(names(Y))
test_len = length(test_idx)
abs_mat = reduce(hcat, [combined_abs[y] for y in indicator_order])
rel_mat = reduce(hcat, [combined_rel[y] for y in indicator_order])

default(fontfamily = "Arial")
# Wrap indicator names every two words to shorten row height
fix_names = replace.(indicator_order, "Particulate Matter Formation Potentia" => "Particulate Matter Formation Potential")
pair_lines(s) = begin
    w = split(s, ' ')
    pairs = [join(w[i:min(i+1, end)], " ") for i in 1:2:length(w)]
    join(pairs, "\n")
end
ylabels_full = pair_lines.(fix_names)
ypos_res = 1:length(indicator_order)
xtick_pos = 1:5:test_len
xtick_lbl = string.(xtick_pos)
# Absolute residual heatmap (same style as coefficient plot)
abs_heat = heatmap(1:test_len, ypos_res, permutedims(abs_mat);
                   xticks = (xtick_pos, xtick_lbl),
                   yticks = (ypos_res, ylabels_full),
                   xlabel = "Test Sample Index", ylabel = "Environmental Impact",
                   title = "Residual Heatmap of Multi-Response Linear Surrogate (Absolute, Test Set)\nColumns = Test Samples (n=$(test_len)); Rows = Environmental Impacts; Color = Absolute Residual",
                   size = (1400, 900),
                   left_margin = 10mm, bottom_margin = 12mm, top_margin = 10mm,
                   xguidefont = font("Arial", 16), yguidefont = font("Arial", 14),
                   tickfont = font("Arial", 10), titlefont = font("Arial", 18),
                   color = :coolwarm, legend = false, colorbar = true,
                   xgrid = true, ygrid = true, grid = true,
                   gridalpha = 1.0, gridcolor = :black, gridlinewidth = 1,
                   foreground_color_grid = :black, gridstyle = :solid)
nx = test_len; ny = length(indicator_order)
plot!(abs_heat; xlims = (0.5, nx + 0.5), ylims = (0.5, ny + 0.5))
for x in 0.5:1:(nx + 0.5)
    vline!(abs_heat, [x]; color = :black, lw = 1.0, alpha = 1.0)
end
for y in 0.5:1:(ny + 0.5)
    hline!(abs_heat, [y]; color = :black, lw = 1.0, alpha = 1.0)
end
savefig(abs_heat, "plots_per_cubic_meter_AnD_feed_in/residuals_abs_heatmap.svg")

# —— Residual heatmap (Test set, raw residuals) ——

indicator_order = collect(names(Y))
test_len = length(test_idx)

# Use the raw residual matrix rather than absolute values
# If your raw residual dictionary has a different name, swap combined_raw accordingly
res_mat = reduce(hcat, [combined_abs[y] for y in indicator_order])

default(fontfamily = "Arial")

# Split indicator names into two lines to reduce row height
fix_names = replace.(indicator_order,
    "Particulate Matter Formation Potentia" =>
    "Particulate Matter Formation Potential")

pair_lines(s) = begin
    w = split(s, ' ')
    pairs = [join(w[i:min(i+1, end)], " ") for i in 1:2:length(w)]
    join(pairs, "\n")
end
ylabels_full = pair_lines.(fix_names)

# Keep only the abbreviation after the dash for ticks (e.g., "... - SFP" -> "SFP")
short_label(s) = begin
    m = match(r"-\s*(\S+)\s*$", s)
    m === nothing ? strip(s) : m.captures[1]
end
ylabels_abbr = short_label.(indicator_order)

ypos_res  = 1:length(indicator_order)
xtick_pos = 1:5:test_len
xtick_lbl = string.(xtick_pos)

# Symmetric color scale centered at 0 to show positive/negative residuals
maxabs_res = maximum(abs, res_mat)

res_heat = heatmap(1:test_len, ypos_res, permutedims(res_mat);
    xticks = (xtick_pos, xtick_lbl),
    yticks = (ypos_res, ylabels_abbr),
    xlabel = "Test sample index",
    ylabel = "Environmental impact",
    title  = "Raw Residual Heatmap of Multi-Response Linear Surrogate (Test Set)",
    size   = (1400, 900),
    left_margin   = 10mm,
    bottom_margin = 12mm,
    top_margin    = 10mm,
    xguidefont = font("Arial", 18),
    yguidefont = font("Arial", 18),
    tickfont   = font("Arial", 18),
    titlefont  = font("Arial", 18),
    color   = :coolwarm,
    clims   = (-maxabs_res, maxabs_res),   # critical: symmetric around zero
    legend  = false,
    colorbar = true,
    xgrid = false, ygrid = false, grid = false)

# If grid lines feel too noisy, keep them off (already grid=false).
# If you still want thin borders, add an outer frame without per-cell grids:

nx = test_len
ny = length(indicator_order)
plot!(res_heat; xlims = (0.5, nx + 0.5), ylims = (0.5, ny + 0.5))
# Draw only the outer border
vline!(res_heat, [0.5, nx + 0.5]; color = :black, lw = 1.0)
hline!(res_heat, [0.5, ny + 0.5]; color = :black, lw = 1.0)

savefig(res_heat, "plots_per_cubic_meter_AnD_feed_in/residuals_raw_heatmap.svg")
println("✅ Residual heatmap saved: plots_per_cubic_meter_AnD_feed_in/residuals_raw_heatmap.svg")

# Relative residual heatmap
max_rel = maximum(abs, rel_mat)
rel_heat = heatmap(1:test_len, ypos_res, permutedims(rel_mat);
                   xticks = (xtick_pos, xtick_lbl),
                   yticks = (ypos_res, ylabels_abbr),
                   xlabel = "Test Sample Index", ylabel = "Environmental Impact",
                   title = "Relative Residual Heatmap of Multi-Response Linear Surrogate (Test Set)",
                   size = (1400, 900),
                   left_margin = 10mm, bottom_margin = 12mm, top_margin = 10mm,
                   xguidefont = font("Arial", 18), yguidefont = font("Arial", 18),
                   tickfont = font("Arial", 18), titlefont = font("Arial", 18),
                   color = :coolwarm, clims = (-max_rel, max_rel),
                   legend = false, colorbar = true,
                   xgrid = false, ygrid = false, grid = false)
plot!(rel_heat; xlims = (0.5, nx + 0.5), ylims = (0.5, ny + 0.5))
# Keep only the outer border, no internal grid
vline!(rel_heat, [0.5, nx + 0.5]; color = :black, lw = 1.0)
hline!(rel_heat, [0.5, ny + 0.5]; color = :black, lw = 1.0)
savefig(rel_heat, "plots_per_cubic_meter_AnD_feed_in/residuals_rel_heatmap.svg")
println("✅ Heatmaps saved: abs/rel heatmap to plots_per_cubic_meter_AnD_feed_in/")

# —— R² bar chart (test 10%) ——
r2_labels = replace.(indicator_order, " " => "\n")
r2_labels = replace.(r2_labels, "Particulate\nMatter\nFormation\nPotentia" => "Particulate\nMatter\nFormation\nPotential")
r2_vals = [r2_dict[ind] for ind in indicator_order]
r2_pos = 1:length(indicator_order)
r2_bar = bar(r2_pos, r2_vals;
             xticks = (r2_pos, r2_labels),
             xlabel = "Environmental Impact", ylabel = "R² (Test)",
             title = "R² of Multi-Response Linear Surrogate (Test Set)\n80/10/10 split; metrics computed on held-out 10% test",
             size = (1400, 600),
             left_margin = 10mm, bottom_margin = 12mm, top_margin = 10mm,
             xguidefont = font("Arial", 16), yguidefont = font("Arial", 16),
             tickfont = font("Arial", 10), titlefont = font("Arial", 18),
             legend = false, color = :steelblue,
             ylims = (0, 1.05),
             xgrid = false, ygrid = true, grid = true,
             gridalpha = 1.0, gridcolor = :black, gridlinewidth = 1,
             foreground_color_grid = :black, gridstyle = :solid)
# Annotate exact R² values on the bars
for i in eachindex(r2_pos)
    ylab = max(r2_vals[i] - 0.04, 0.02)  # place inside the bar near the top to avoid tick overlap
    annotate!(r2_bar, r2_pos[i], ylab,
              text(string(round(r2_vals[i], digits = 6)), 9, :black, :center))
end
savefig(r2_bar, "plots_per_cubic_meter_AnD_feed_in/r2_bar_test.svg")
println("✅ R² bar chart saved: plots_per_cubic_meter_AnD_feed_in/r2_bar_test.svg")

# —— Coefficient matrix heatmap (multi-response linear surrogate B matrix) —— #

coef_rows = string.(names(X_dedup))        # predictors (no intercept)

desired_cols = [
    "Acidification Potential - AP",
    "Eutrophication Potential - EP",
    "Fossil Depletion Potential - FDP",
    "Global Climate Change Potential - GCCP",
    "Particulate Matter Formation Potentia - PMFP",
    "Cumulative Energy Demand - CED",
    "Smog Formation Potential - SFP",
    "Water Use"
]

# Keep all indicators but order them in a sensible way
col_order = vcat([c for c in desired_cols if c in names(Y)],
                 [c for c in names(Y) if c ∉ desired_cols])

col_idx = [findfirst(==(c), names(Y)) for c in col_order]

# Raw coefficient matrix (drop intercept row)
coef_mat_raw = beta_final[2:end, col_idx]

# Log10 compression: sign(β)*log10(1 + |β|)
coef_mat = sign.(coef_mat_raw) .* log10.(1 .+ abs.(coef_mat_raw))
maxabs = maximum(abs, coef_mat)

# Colorbar ticks (symmetric around 0, labeled as 10^k)
tick_max = floor(Int, maxabs)
colorbar_ticks_vals   = tick_max == 0 ? [0] : collect(-tick_max:tick_max)
colorbar_tick_labels  = [v == 0 ? "0" :
                         v > 0  ? "10^$(v)" :
                                  "-10^$(abs(v))" for v in colorbar_ticks_vals]

# Short, multi-line x-labels: AP, EP, etc.
col_label_map = Dict(
    "Acidification Potential - AP"              => "AP",
    "Eutrophication Potential - EP"            => "EP",
    "Fossil Depletion Potential - FDP"         => "FDP",
    "Global Climate Change Potential - GCCP"   => "GCCP",
    "Particulate Matter Formation Potentia - PMFP" =>
        "PMFP",
    "Cumulative Energy Demand - CED"           => "CED",
    "Smog Formation Potential - SFP"           => "SFP",
    "Water Use"                                => "WU"
)
xlabels = [get(col_label_map, c, c) for c in col_order]
xpos    = 1:length(col_order)

# Predictor display names
pred_map = Dict(
    "steel"   => "Steel (CAPEX)",
    "ng_LCI"  => "Natural gas consumption",
    "AH"      => "Avoided heat",
    "AE"      => "Avoided electricity",
    "SSO_LCI" => "SSO feed to digester"
)
coef_rows_disp = [get(pred_map, r, r) for r in coef_rows]
pair_two_words(s) = begin
    w = split(s, ' ')
    parts = [join(w[i:min(i+1, end)], " ") for i in 1:2:length(w)]
    join(parts, "\n")
end
coef_rows_disp_wrap = pair_two_words.(coef_rows_disp)
# Force avoided heat/electricity onto separate lines
coef_rows_disp_wrap = replace.(coef_rows_disp_wrap,
    "Avoided heat" => "Avoided\nheat",
    "Avoided electricity" => "Avoided\nelectricity",
    "Steel (CAPEX)" => "Steel\n(CAPEX)")
yposeq         = 1:length(coef_rows_disp)

default(fontfamily = "Arial")

heat_coef = heatmap(
    xpos, yposeq, coef_mat;
    xticks          = (xpos, xlabels),
    yticks          = (yposeq, coef_rows_disp_wrap),
    xrotation       = 0,
    size            = (1000, 650),
    left_margin     = 10mm,
    bottom_margin   = 8mm,
    top_margin      = 5mm,
    guidefont       = font("Arial", 14),   # axis labels
    tickfont        = font("Arial", 11),   # tick labels
    titlefont       = font("Arial", 16),
    xlabel          = "LCIA category",
    ylabel          = "Predictor (per-FU LCI flow)",
    title           = "Coefficients of multi-response linear surrogate",
    color           = :coolwarm,
    clims           = (-maxabs, maxabs),
    legend          = false,
    colorbar        = true,
    colorbar_title  = "sign * log10(1 + |beta|)",
    colorbar_titlefont = font("Arial", 12),
    colorbar_tickfont  = font("Arial", 10),
    colorbar_ticks  = (colorbar_ticks_vals, colorbar_tick_labels),
    xgrid           = false,
    ygrid           = false,
    grid            = false
)

# Set plot limits to align with centers of cells
nx = length(col_order)
ny = length(coef_rows_disp)
plot!(heat_coef; xlims = (0.5, nx + 0.5), ylims = (0.5, ny + 0.5))

# Optional: very light white grid just to separate cells (not heavy black boxes)
for x in 1:nx-1
    vline!(heat_coef, [x + 0.5]; color = RGBA(1,1,1,0.4), lw = 0.5)
end
for y in 1:ny-1
    hline!(heat_coef, [y + 0.5]; color = RGBA(1,1,1,0.4), lw = 0.5)
end

# Annotate only "large" coefficients to avoid clutter
annot_threshold = 0.05          # adjust to 0.1 or 1 if you want fewer/more labels
for i in 1:ny, j in 1:nx
    val = coef_mat_raw[i, j]
    if abs(val) > annot_threshold
        annotate!(
            heat_coef,
            j, i,
            text(string(round(val, digits = 2)), 8, :black, :center)
        )
    end
end

savefig(heat_coef,
        "plots_per_cubic_meter_AnD_feed_in/coefficients_heatmap.svg")
println("✅ Saved coefficient heatmap: plots_per_cubic_meter_AnD_feed_in/coefficients_heatmap.svg")

# --- Export split indices for reproducible diagnostics across scripts ---
split_df = DataFrame(
    index = 1:N,
    time  = time_col,
    split = fill("train", N)
)
split_df.split[val_idx] .= "val"
split_df.split[test_idx] .= "test"

CSV.write("train_val_test_split.csv", split_df)
println("✅ Saved split file: train_val_test_split.csv")
