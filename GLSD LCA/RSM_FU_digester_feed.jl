using CSV, DataFrames, GLM, StatsModels, Statistics, StatsBase, Random, LinearAlgebra

# --------- Configuration: input filenames (match your existing pipeline names) ---------
const PREDICTOR_FILE = "cleaned_data_per_cubic_AD_feed.csv"           # from LCI_Calculation_FU_PerCubicMeter_AD_feedin V2.jl
const LCIA_FILE      = "lcia_results_per_cubic_meter_AnD_feed_in.csv" # from Doing LCA Auto Per Cubic Meter AnD feed in.py
const RAW_FILE       = "raw_inputs_with_feed.csv"                     # from LCI_Calculation_FU_PerCubicMeter_AD_feedin V2.jl

# =============== Utilities: VIF and correlation pruning ===============

"Compute VIF for each DataFrame column (regress that column on the others; VIF = 1 / (1 - R²))"
function calc_vif(df::DataFrame)
    X = Matrix{Float64}(df)
    n, p = size(X)
    res = DataFrame(variable = String[], vif = Float64[])
    for j in 1:p
        y = X[:, j]
        others = setdiff(1:p, [j])
        Xj = X[:, others]
        Xdesign = hcat(ones(n), Xj)
        β = Xdesign \ y
        ŷ = Xdesign * β
        ss_resid = sum((y .- ŷ).^2)
        ss_tot   = sum((y .- mean(y)).^2)
        R2  = 1 - ss_resid / ss_tot
        vif = (1 - R2) ≈ 0 ? Inf : 1 / (1 - R2)
        push!(res, (string(names(df)[j]), vif))
    end
    rename!(res, [:variable, :vif])
    return res
end

"Simple correlation pruning: drop the latter column when |ρ| ≥ threshold"
function drop_high_corr(df::DataFrame; threshold = 0.99)
    X = Matrix{Float64}(df)
    p = size(X, 2)
    C = cor(X)
    nm = names(df)
    drop = Set{String}()
    for i in 1:(p-1), j in (i+1):p
        if abs(C[i, j]) ≥ threshold
            push!(drop, nm[j])  # keep the first occurrence, drop later ones
        end
    end
    return drop
end

# =============== Step 0: Load data & align complete cases ===============

println("Reading data...")
to_symbol_names!(df::DataFrame) = (rename!(df, Symbol.(names(df))); df)

X_raw  = CSV.read(PREDICTOR_FILE, DataFrame) |> to_symbol_names!  # Predictors: SSO_LCI, AE, AH, CH4, NG_LCI, FeCl3, CHP, steel, ...
Y_all  = CSV.read(LCIA_FILE, DataFrame)        |> to_symbol_names! # Dependents: index, time, LCIA indicators
raw_in = CSV.read(RAW_FILE, DataFrame)         |> to_symbol_names! # Raw quantities: dig_tfeed, etc.

# Use only LCIA indicator columns for completecases (drop index, time)
Y_ind = select(Y_all, Not([:index, :time]); copycols = false)
mask  = completecases(Y_ind)

X_raw  = X_raw[mask, :]
Y_all  = Y_all[mask, :]
raw_in = raw_in[mask, :]

N = nrow(X_raw)
@assert N == nrow(Y_all) == nrow(raw_in) "Row counts after masking do not match."

index_col = X_raw.index
time_col  = X_raw.time

# S_t: digester feed volume, m³/day (dig_tfeed is kgal/day)
kgal_to_m3 = 3.78541
@assert "dig_tfeed" in names(raw_in) "raw_inputs_with_feed.csv needs a dig_tfeed (kgal/day) column."
S_t = raw_in.dig_tfeed .* kgal_to_m3

# =============== Step 1: Build candidate X, compute VIF, and prune correlations ===============

# Drop index, time, and two optional credit columns if present
drop_cols = String[]
for nm in names(X_raw)
    if nm in ("index", "time", "A_SSO_WTE", "A_SSO_LF")
        push!(drop_cols, nm)
    end
end
X_candidates = select(X_raw, Not(drop_cols))

println("Computing VIF on raw candidate predictors...")
vif_table = calc_vif(X_candidates)
CSV.write("vif_prescreen_FU.csv", vif_table)
println("  -> wrote vif_prescreen_FU.csv")

println("Pruning highly correlated predictors (|ρ| ≥ 0.99)...")
drop_set = drop_high_corr(X_candidates; threshold = 0.99)
X_dedup  = select(X_candidates, Not(collect(drop_set)))
CSV.write("cleaned_data_dedup_FU_digester_feed.csv", X_dedup)
println("  -> wrote cleaned_data_dedup_FU_digester_feed.csv")

# =============== Step 2: LCIA Y matrix (keep indicator columns only) ===============
# Re-select indicator columns from filtered Y_all to avoid using an unfiltered view
Y = select(Y_all, Not(["index", "time"]); copycols = false)
indicator_names = names(Y)    # LCIA indicator names; regress each later

# =============== Step 3: 80/10/10 train/val/test split ===============

Random.seed!(2025)            # Align with the reproducibility manifest in the paper
idx_all = collect(1:N)
shuffle!(idx_all)

n_train = round(Int, 0.80 * N)
n_val   = round(Int, 0.10 * N)

train_idx = idx_all[1:n_train]
val_idx   = idx_all[n_train+1 : n_train+n_val]
test_idx  = idx_all[n_train+n_val+1 : end]

partition = fill("train", N)
partition[val_idx]  .= "val"
partition[test_idx] .= "test"

# =============== Step 4: Build panel_FU_digester_feed.csv ===============

panel = DataFrame(
    index     = index_col,
    date      = time_col,
    S_t       = S_t,          # Volume for 1 m³ AnD feed FU, m³/day
    partition = partition,
)

# Add predictors used for modeling (all columns in X_dedup)
for nm in names(X_dedup)
    panel[!, nm] = X_dedup[!, nm]
end

# Add LCIA indicators with a unified _FU suffix (e.g., GCCP_FU, CED_FU, ...)
for nm in indicator_names
    newname = string(nm) * "_FU"
    panel[!, newname] = Y[!, nm]
end

CSV.write("panel_FU_digester_feed.csv", panel)
println("  -> wrote panel_FU_digester_feed.csv")

# =============== Step 5: OLS main-effect multiple linear regression; output mlr_results_FU.csv ===============

X_train = X_dedup[train_idx, :]
X_val   = X_dedup[val_idx, :]
X_test  = X_dedup[test_idx, :]

Y_train = Y[train_idx, :]
Y_val   = Y[val_idx, :]
Y_test  = Y[test_idx, :]

summary_table = DataFrame(
    Indicator = String[],
    R2_test   = Float64[],
    RMSE_test = Float64[]
)

full_summary_table = DataFrame(
    Indicator = String[],
    Variable  = String[],
    Coef      = Float64[],
    StdError  = Float64[],
    t         = Float64[],
    P         = Float64[],
    R2_test   = Float64[],
    RMSE_test = Float64[],
    N_test    = Int[]
)

predictor_terms = Term.(Symbol.(names(X_dedup)))  # e.g. SSO_LCI, AE, AH, CH4, NG_LCI, FeCl3, CHP, steel, ...

for nm in indicator_names
    y_train = Y_train[!, nm]
    y_test  = Y_test[!, nm]

    df_train = hcat(DataFrame(y = y_train), X_train)
    df_test  = hcat(DataFrame(y = y_test),  X_test)

    form = Term(:y) ~ sum(predictor_terms)
    model = lm(form, df_train)

    # --- test-set metrics ---
    yhat_test = predict(model, df_test)
    rss  = sum((y_test .- yhat_test).^2)
    tss  = sum((y_test .- mean(y_test)).^2)
    r2_test   = 1 - rss / tss
    rmse_test = sqrt(mean((y_test .- yhat_test).^2))
    n_test    = length(y_test)

    push!(summary_table, (string(nm), r2_test, rmse_test))

    ct = coeftable(model)
    coef_vals   = ct.cols[1]
    stderr_vals = ct.cols[2]
    tvals       = ct.cols[3]
    pvals       = ct.cols[4]
    var_names   = ct.rownms

    for i in eachindex(var_names)
        push!(full_summary_table, (
            string(nm),          # Indicator
            string(var_names[i]),# Variable (includes (Intercept))
            coef_vals[i],
            stderr_vals[i],
            tvals[i],
            pvals[i],
            r2_test,
            rmse_test,
            n_test
        ))
    end
end

CSV.write("model_summary_FU.csv", summary_table)   # Short version: one row per indicator
CSV.write("mlr_results_FU.csv", full_summary_table)
println("  -> wrote model_summary_FU.csv and mlr_results_FU.csv")

# =============== Step 6: split_summary_FU.csv (representativeness table) ===============

@assert "SSO_LCI" in names(panel) "panel needs an SSO_LCI column."

split_summary = DataFrame(
    partition    = String[],
    N_days       = Int[],
    SSO_LCI_mean = Float64[],
    SSO_LCI_sd   = Float64[],
    S_t_mean     = Float64[],
    S_t_sd       = Float64[]
)

for (name, idxs) in [("train", train_idx), ("val", val_idx), ("test", test_idx)]
    push!(split_summary, (
        name,
        length(idxs),
        mean(panel.SSO_LCI[idxs]),
        std(panel.SSO_LCI[idxs]),
        mean(panel.S_t[idxs]),
        std(panel.S_t[idxs])
    ))
end

CSV.write("split_summary_FU.csv", split_summary)
println("  -> wrote split_summary_FU.csv")

println("Done.")
