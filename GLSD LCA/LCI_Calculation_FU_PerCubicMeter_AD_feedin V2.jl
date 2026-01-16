##############################
# Preamble & data loading
##############################

using SQLite, DataFrames, Statistics, Plots, CSV, Dates, Missings, XLSX, IfElse

include("./Transformer_Model.jl")

# Load database (fallback to ADM1 code directory if not found in current folder)
function load_ad_database()
    candidates = [
        joinpath(@__DIR__, "AD_database.db"),
        joinpath(@__DIR__, "ADM1 code", "AD_database.db")
    ]

    last_err = nothing
    for path in candidates
        if !isfile(path) || filesize(path) == 0
            continue
        end

        db = SQLite.DB(path)
        try
            full_table_info = SQLite.tables(db)
            data  = DBInterface.execute(db, "SELECT * FROM ad_data")   |> DataFrame
            names = DBInterface.execute(db, "SELECT * FROM title_data") |> DataFrame
            return full_table_info, data, names, path
        catch err
            last_err = err
        finally
            DBInterface.close!(db)
        end
    end

    error("Could not find a database containing ad_data/title_data tables. Last error: $(last_err)")
end

full_table_info, data, names, db_path = load_ad_database()
println("Using database: ", db_path)

##############################
# Raw variable preparation
##############################

dates = copy(data.Date) # Date in String form
dig1_pH  = copy(data.v7001) # pH units
dig2_pH  = copy(data.v7031) # pH units
dig3_pH  = copy(data.v7061) # pH units
dig4_ph  = copy(data.v7062) # pH units
dig1_TS  = copy(data.v7007) # Total solids; mg/L
dig2_TS  = copy(data.v7037) # Total solids; mg/L
dig3_TS  = copy(data.v7067) # Total solids; mg/L
dig1_pTS = copy(data.v7006) # Percent total solids; %
dig2_pTS = copy(data.v7036) # Percent total solids; %
dig3_pTS = copy(data.v7066) # Percent total solids; %
dig1_pVS = copy(data.v7011) # Percent volatile solids; %
dig2_pVS = copy(data.v7041) # Percent volatile solids; %
dig3_pVS = copy(data.v7071) # Percent volatile solids; %
dig4_PVS = copy(data.v7072) # Percent volatile solids; %
dig1_alk = copy(data.v7016) # Alkalinity; mg/L
dig2_alk = copy(data.v7046) # Alkalinity; mg/L
dig3_alk = copy(data.v7076) # Alkalinity; mg/L
dig4_alk = copy(data.v7077) # Alkalinity; mg/L
dig1_amm = copy(data.v7019) # Ammonia; mg/kg dry
dig2_amm = copy(data.v7054) # Ammonia; mg/kg dry
dig3_amm = copy(data.v7070) # Ammonia; mg/kg dry
dig1_VA  = copy(data.v7021) # Volatile acids; mg/L
dig2_VA  = copy(data.v7051) # Volatile acids; mg/L
dig3_VA  = copy(data.v7081) # Volatile acids; mg/L
dig4_VA  = copy(data.v7082) # Volatile acids; mg/L
dig1_VAalk = copy(data.v7026) # VA/Alk ratio
dig2_VAalk = copy(data.v7056) # VA/Alk ratio
dig3_VAalk = copy(data.v7086) # VA/Alk ratio
dig4_VAalk = copy(data.v7087) # VA/Alk ratio
dig1_CODt = copy(data.v7017) # Total COD; mg/L
dig1_CODs = copy(data.v7018) # Soluble COD; mg/L
dig2_CODt = copy(data.v7052) # Total COD; mg/L
dig2_CODs = copy(data.v7053) # Soluble COD; mg/L
dig3_CODt = copy(data.v7068) # Total COD; mg/L
dig3_CODs = copy(data.v7069) # Soluble COD; mg/L
dig1_HRT = copy(data.v7141) # Hydraulic retention time; days
dig2_HRT = copy(data.v7146) # Hydraulic retention time; days
dig3_HRT = copy(data.v7151) # Hydraulic retention time; days
dig1_feed = copy(data.v7106) # Digester feed; kgal
dig2_feed = copy(data.v7111) # Digester feed; kgal
dig3_feed = copy(data.v7116) # Digester feed; kgal
dig_tfeed = copy(data.v7121) # Total feed; kgal (overall digester inflow)
dig_PSSO = copy(data.v15932) # PSSO to digesters; (early records mix gal / kgal)
out_gas1 = copy(data.v7231)  # Total biogas production; k ft^3
out_gas2 = copy(data.v16342) # Total biogas production; k ft^3 (alternate)
GBT_cake = copy(data.v11223) # Cake % to gravity belt thickener; mg/L
GBT1 = copy(data.v9126) # Gravity belt thickening #1; %
GBT2 = copy(data.v9131) # Gravity belt thickening #2; %
food_pH = copy(data.v13027) # Food waste tank; pH units
food_pTS = copy(data.v13022) # Percent total solids; %
food_pVS = copy(data.v13024) # Percent volatile solids; %
food_VA = copy(data.v13025) # Volatile acids; mg/L
food_amm = copy(data.v13030) # Ammonia; mg/kg dry
food_pho = copy(data.v13041) # Phosphorus; mg/kg dry
food_CODt = copy(data.v13028) # Total COD; mg/L
food_CODs = copy(data.v13029) # Soluble COD; mg/L
blend_pH = copy(data.v12637) # Blend tank pH; pH units
blend_pTS = copy(data.v12632) # Percent total solids; %
blend_VA = copy(data.v12635) # Volatile acids; mg/L
blend_pVS = copy(data.v12634) # Percent volatile solids; %
blend_alk = copy(data.v12636) # Alkalinity; mg/L
sludge_pVS = copy(data.v9076) # Percent volatile solids in primary sludge; %
centrifuge_feed = copy(data.v9196) # Centrifuge feed; kgal (?)
centrifuge_cakesolids = copy(data.v9011) # Centrifuge cake solids; (units messy, fixed later)
unknown = copy(data.v9211) # "C1-GLSD Weight Scale WT [wet ton]" (??)
centrate_TS = copy(data.v9122) # Centrate total solids; mg/L
centrifuge1_solids = copy(data.v9111) # "Solids centrifuge feed #1" (?)
centrifuge2_solids = copy(data.v9112) # "Solids centrifuge feed #2" (?)
centrifuge_cakesludge = copy(data.v9108) # Centrifuge cake sludge; % (?)
secondary_solids = copy(data.v9161) # Solids secondary flow; kgal
PSSO_flow = copy(data.v13084) # PSSO yesterday's flow; gal
eng1_ch4 = copy(data.v13991) # Methane percent to engine 1; %
eng2_ch4 = copy(data.v14011) # Methane percent to engine 2; %

##############################
# Gas, PSSO, and cake cleaning
##############################

# Fill out_gas1 with out_gas2 where available; prefer out_gas2
out_gas = copy(out_gas2)
out_gas[ismissing.(out_gas)] = out_gas1[ismissing.(out_gas)]

# Set outliers to missing: 0 or >= 1e4 kft³
for i in eachindex(out_gas)
    if !ismissing(out_gas[i]) && (out_gas[i] == 0.0 || out_gas[i] >= 1e4)
        out_gas[i] = missing
    end
end

# Early PSSO mixes gal/kgal: treat >=100 as gal and divide by 1000 → kgal
for i in eachindex(dig_PSSO)
    if !ismissing(dig_PSSO[i]) && dig_PSSO[i] >= 100
        dig_PSSO[i] = dig_PSSO[i] / 1000
    end
end

# Apply the same simple unit correction to cake solids (>=100 treated as ×1000 typo)
for i in eachindex(centrifuge_cakesolids)
    if !ismissing(centrifuge_cakesolids[i]) && centrifuge_cakesolids[i] >= 100
        centrifuge_cakesolids[i] = centrifuge_cakesolids[i] / 1000
    end
end

##############################
# Time and flow (EPA data expansion)
##############################

times = Date.(dates, dateformat"mm/dd/yyyy")

"""
    expand_monthly_to_daily(monthly::Vector{<:Real}, is_leap_year::Bool)

Expand 12 monthly totals into a daily vector with uniform distribution.
"""
function expand_monthly_to_daily(monthly::Vector{<:Real}, is_leap_year::Bool)
    @assert length(monthly) == 12 "monthly vector must have 12 elements"

    days_in_month = [31, is_leap_year ? 29 : 28, 31, 30, 31, 30,
                     31, 31, 30, 31, 30, 31]

    daily = Float64[]

    for (month_total, days) in zip(monthly, days_in_month)
        append!(daily, fill(month_total / days, days))
    end

    return daily
end

monthly_data_2017 = [
    907.19, 790.55, 923.3, 1367.92, 1006.58, 925.22,
    666.58, 565.02, 610.77, 599.34, 688.19, 660.83
]

monthly_data_2018 = [
    837.84, 900.93, 1233.86, 1191.6, 852.02, 642.45,
    638.3, 776.04, 804.82, 878.51, 1799.55, 1195.08
]

monthly_data_2019 = [
    1016.84, 921.07, 1075.18, 1120.05, 1068.94, 855.6,
    779.12, 664.48, 520.14, 714.25, 718.98, 1089.93
]

monthly_data_2020 = [
    931.98, 831.57, 892.49, 1162.12, 825.71, 552.41,
    570.11, 526.03, 493.34, 565.93, 615.44, 991.41
]

monthly_data_2021 = [
    853.86, 663.53, 866.70, 886.42, 937.90, 736.89,
    1295.92, 888.19, 968.60, 914.39, 969.74, 841.90
]

monthly_data_2022 = [
    814.99, 967.43, 964.04, 892.44, 713.31, 584.42,
    528.18, 492.61, 559.53, 599.49, 585.05, 824.85
]

all_daily_influent_flowrate = vcat(
    fill(23, 365*6+2), # early years: assume 23 MGD
    expand_monthly_to_daily(monthly_data_2017, false),
    expand_monthly_to_daily(monthly_data_2018, false),
    expand_monthly_to_daily(monthly_data_2019, false),
    expand_monthly_to_daily(monthly_data_2020, true),
    expand_monthly_to_daily(monthly_data_2021, false),
    expand_monthly_to_daily(monthly_data_2022, false)
)

# GLSD-estimated daily influent (kgal → Mgal)
GLSD_estimate_daily_influent_flowrate = dig_tfeed .* 1e-3  # Mgal/day

# Guard against 0 or negative values to avoid divide-by-zero later
flow_safe = similar(GLSD_estimate_daily_influent_flowrate, Union{Missing, Float64})

for i in eachindex(GLSD_estimate_daily_influent_flowrate)
    v = GLSD_estimate_daily_influent_flowrate[i]
    if !ismissing(v) && v > 0
        flow_safe[i] = v
    else
        flow_safe[i] = missing
    end
end

# Convert to baseline volume m³/day
den_safe = flow_safe .* 3785.41  # m³/day

##############################
# CAPEX & SSO LCI
##############################

# Time series of total feed + PSSO
a1 = plot(times, dig_tfeed, ylabel = "(kgal)",
          label = "Total Feed",
          legend = :topleft,
          title = "Digester PSSO & Total Feed Over Time")
plot!(a1, times, dig_PSSO, label = "PSSO")

# Assume SSO density 3.99 kg/gal; compute SSO LCI (kg/m³ wastewater)
SSO_LCI = dig_PSSO .* 1000 .* 3.99 ./ den_safe
A_SSO_WTE = SSO_LCI .* 0.68
A_SSO_LF  = SSO_LCI .* 0.32

# CAPEX allocation (per m³)
building = 3.67e-6 * 23 ./ flow_safe
concrete = 7.37e-7 * 23 ./ flow_safe
earth    = 4.44e-6 * 23 ./ flow_safe
gravel   = 3.60e-4 * 23 ./ flow_safe
steel    = 4.91e-5 * 23 ./ flow_safe
elec     = 0.065    * 23 ./ flow_safe

a11 = plot(times, SSO_LCI, label = "SSO",
           ylabel = "(SSO kg per m^3 wastewater treated)",
           title = "SSO Over Time")
hline!(a11, [2.06], label="EPA Partial Capacity")
hline!(a11, [4.13], label="EPA Full Capacity")

##############################
# Biogas / CHP / AE / AH
##############################

a2 = plot(times, out_gas .* 1000,
          ylabel = "(ft^3/day)",
          label = "GLSD Total Gas Production",
          title = "Biogas Production Per Day Over Time")
hline!(a2, [7.98e5], label="EPA Partial Low")
hline!(a2, [1.28e6], label="EPA Full Low")
hline!(a2, [1.78e6], label="EPA Full Base")
hline!(a2, [1.11e6], label="EPA Partial Base")
hline!(a2, [3.93e5], label="EPA Base")

a3 = plot(times,
          out_gas .* 0.0283168 .* 1000 ./ den_safe,
          label ="GLSD Biogas Data",
          ylabel = "(m^3 biogas per m^3 wastewater treated)")
hline!(a3, [0.428], label="Full Low")
hline!(a3, [0.598], label="Full Base")
hline!(a3, [0.371], label="Partial Base")
hline!(a3, [0.267], label="Partial Low")

a4 = plot(times, eng1_ch4,
          label ="CHP Engine 1 Methane Percent",
          ylabel = "(%)", xlabel = "Time",
          legend = :topleft)
plot!(a4, times, eng2_ch4,
      label ="CHP Engine 2 Methane Percent")

a5 = plot(times, dig_PSSO .* 1000,
          label ="Total gallons of PSSO to digesters",
          ylabel = "gal", xlabel = "Time",
          legend = :topleft)
plot!(a5, times, PSSO_flow,
      label ="PSSO yesterday's flow")

a6 = plot(times, centrifuge_cakesolids,
          label = "Centrifuge Cake Solids",
          ylabel = "(Short Dry Tons/day)",
          xlabel = "Time", legend = :topleft)

# Gas and energy required for the pelletization line
gas_to_pelletization   = centrifuge_cakesolids .* 8462 ./ (550/948) # ft^3/day
energy_to_pelletization = centrifuge_cakesolids .* 8462            # MJ/day

a7 = plot(times, gas_to_pelletization,
          label = "Gas to Pelletization",
          ylabel = "(ft^3/day)", xlabel = "Time")

a71 = plot(times,
           gas_to_pelletization ./ 35.5 ./ den_safe,
           label = "GLSD Data",
           ylabel = "(m^3/m^3 wastewater)",
           xlabel = "Time",
           legend = :topleft,
           title = "Gas to Pelletization Per m^3 Wastewater Treated")
hline!(a71, [0.105], label="EPA Report Partial")
hline!(a71, [0.142], label="EPA Report Full")

a8 = plot(times, out_gas .* 1000 .- gas_to_pelletization,
          label = "Gas to CHP",
          ylabel = "(ft^3/day)", xlabel = "Time",
          legend = :topleft, title = "Gas to CHP Over Time")

a81 = plot(times, out_gas .* 1000,
           label = "Total Gas Production",
           ylabel = "(ft^3/day)", xlabel = "Time",
           legend = :topleft, title = "Total Gas Production Over Time")
plot!(a81, times, gas_to_pelletization,
      label = "Gas to Pelletization",
      ylabel = "(ft^3/day)", xlabel = "Time")

# CHP heating value: ft³/day → MJ/day
gas_energy_to_CHP = (out_gas .* 1000 .- gas_to_pelletization) .* 550 ./ 948

CHP = (out_gas .* 1000 .- gas_to_pelletization) ./ 35.3 ./ den_safe

a82 = plot(times, CHP, label = "GLSD Data",
           ylabel = "(m^3/m^3 wastewater)",
           xlabel = "Time", legend = :topleft,
           title = "Gas to CHP Per m^3 Wastewater Treated")
hline!(a82, [0.098], label="EPA Report Partial Low")
hline!(a82, [0.184], label="EPA Report Full Low")
hline!(a82, [0.369], label="EPA Report Full Base")
hline!(a82, [0.213], label="EPA Report Partial Base")
hline!(a82, [0],     label="EPA Report Base")

# Avoided heat / electricity
AH = (energy_to_pelletization .+ gas_energy_to_CHP .* 0.39) ./ den_safe
a9 = plot(times, AH, label = "GLSD Data",
          ylabel = "(MJ/m^3 wastewater)",
          xlabel = "Time", legend = :topleft,
          title = "Avoided Heat Over Time")
hline!(a9, [3.14], label="EPA Partial Low")
hline!(a9, [4.74], label="EPA Full Low")
hline!(a9, [6.58], label="EPA Full Base")
hline!(a9, [4.29], label="EPA Partial Base")
hline!(a9, [2.1],  label="EPA Base")

AE = (gas_energy_to_CHP .* 0.40) ./ den_safe
a10 = plot(times, AE, label = "GLSD Data",
           ylabel = "(MJ/m^3 wastewater)",
           xlabel = "Time", legend = :topleft,
           title = "Avoided Electricity Over Time")
hline!(a10, [0.226], label="EPA Partial Low")
hline!(a10, [0.424], label="EPA Full Low")
hline!(a10, [0.848], label="EPA Full Base")
hline!(a10, [0.490], label="EPA Partial Base")
hline!(a10, [0],     label="EPA Base")

##############################
# Methane concentration & CH4 emissions
##############################

a111 = plot(times, eng1_ch4,
            label = "CHP Engine 1 Methane Percent",
            ylabel = "(%)", xlabel = "Time",
            legend = :topleft,
            title = "Methane Percent to CHP Engine 1 & 2 Over Time")
plot!(a111, times, eng2_ch4,
      label = "CHP Engine 2 Methane Percent")

# Methane concentration: average across two engines; if both missing assume 60%
ch4_percent = map(eng1_ch4, eng2_ch4) do a, b
    if ismissing(a) && ismissing(b)
        60.0
    elseif ismissing(a)
        b
    elseif ismissing(b)
        a
    else
        (a + b) / 2
    end
end

# Raw data are 0–100%; convert to 0–1 volume fraction
ch4_frac = ch4_percent ./ 100.0

# Leakage fraction equivalent to /0.95*0.05
slip_frac = 0.05 / 0.95

CH4 = out_gas .* 1000 .* slip_frac .* ch4_frac ./ 35.3 .* 0.656 ./ den_safe

a12 = plot(times, CH4, label = "GLSD Data",
           ylabel = "(kg CH4 / m^3 wastewater)",
           xlabel = "Time", legend = :topleft,
           title = "Methane Emission Over Time")
hline!(a12, [5.22e-3], label="EPA Partial Low")
hline!(a12, [8.43e-3], label="EPA Full Low")
hline!(a12, [0.012],   label="EPA Full Base")
hline!(a12, [7.26e-3], label="EPA Partial Base")
hline!(a12, [2.56e-3], label="EPA Base")

##############################
# FeCl3, heat balance, and natural gas LCI
##############################

# FeCl3 LCI (formula matches your original, just kept in this form)
ferric = (0.25/4.13 .* SSO_LCI .+ 1) .* 19230 ./ (23 * 3785.41 * 365)

a13 = plot(times, ferric, label = "GLSD Data",
           ylabel = "(kg/day)", xlabel = "Time",
           legend = :topleft, title = "FeCl3 Over Time")
hline!(a13, [6.67E-4], label="EPA Partial")
hline!(a13, [7.42E-4], label="EPA Full")
hline!(a13, [5.93E-4], label="EPA Base")

P_digestor_heat = [28083839, 36184946, 43385931]  # MJ/year
Current_digestor_heat =
    (P_digestor_heat[3] - P_digestor_heat[1]) / 4.13 .* SSO_LCI .+ P_digestor_heat[1] # MJ/year

const DAYS_PER_YEAR = 365.0

Pellet_Dryer_Heat_Year_Temp = energy_to_pelletization .* DAYS_PER_YEAR
Facility_Heat_Year_Temp     = 13728414.0
Total_Heat_Year_Temp        = Pellet_Dryer_Heat_Year_Temp .+
                              Current_digestor_heat .+
                              Facility_Heat_Year_Temp

natural_gas_energy = Total_Heat_Year_Temp .-
                     Pellet_Dryer_Heat_Year_Temp .-
                     gas_energy_to_CHP .* 0.39 .* DAYS_PER_YEAR

a14 = plot(times, Current_digestor_heat,
           label = "Current Digestor Heat",
           ylabel = "(MJ/year)", xlabel = "Time",
           legend = :topleft, title = "Current Digestor Heat Over Time")
plot!(a14, times, Pellet_Dryer_Heat_Year_Temp,
      label = "Pellet Dryer Heat")
plot!(a14, times, fill(Facility_Heat_Year_Temp, length(times)),
      label = "Facility Heat")
plot!(a14, times, Total_Heat_Year_Temp,
      label = "Total Heat")
plot!(a14, times, natural_gas_energy,
      label = "Natural Gas Energy")

adjusted_CHP(x) = x >= 0 ? x * 0.39 : x
# Convert MJ/year to MJ/day using the same 365-day assumption
ng_energy_daily = Current_digestor_heat ./ DAYS_PER_YEAR .-
                  passmissing(adjusted_CHP).(gas_energy_to_CHP)

ng_LCI = ng_energy_daily ./ den_safe

# Clip negatives to 0
for i in eachindex(ng_LCI)
    if !ismissing(ng_LCI[i]) && ng_LCI[i] <= 0
        ng_LCI[i] = 0.0
    end
end

a15 = plot(times, ng_LCI, label = "Natural Gas LCI",
           ylabel = "(MJ/m^3)", xlabel = "Time",
           legend = :topleft, title = "Natural Gas LCI Over Time")

##############################
# LCI export
##############################

function get_LCI(index::Int)
    if index < 1 || index > length(times)
        error("Index out of bounds")
    end
    println("Index: ", index)
    println("Time: ", times[index])
    println("CHP: ", CHP[index])
    println("Ferric: ", ferric[index])
    println("SSO LCI: ", SSO_LCI[index])
    println("AE: ", AE[index])
    println("AH: ", AH[index])
    println("CH4: ", CH4[index])
    println("Natural Gas Energy: ", ng_LCI[index])
end

function check_data_if_missing()
    index = Int[]
    times_vec = Any[]
    SSO_vec = Float64[]
    AE_vec = Float64[]
    AH_vec = Float64[]
    CH4_vec = Float64[]
    CHP_vec = Float64[]
    ferric_vec = Float64[]
    ng_LCI_vec = Float64[]
    A_SSO_WTE_vec = Float64[]
    A_SSO_LF_vec  = Float64[]
    building_vec  = Float64[]
    concrete_vec  = Float64[]
    earth_vec     = Float64[]
    gravel_vec    = Float64[]
    steel_vec     = Float64[]
    elec_vec      = Float64[]

    for i in eachindex(times)
        if ismissing(SSO_LCI[i]) || ismissing(AE[i]) || ismissing(AH[i]) ||
           ismissing(CH4[i])     || ismissing(CHP[i]) || ismissing(ferric[i]) ||
           ismissing(ng_LCI[i])
            continue
        end

        push!(index, i)
        push!(times_vec, times[i])
        push!(SSO_vec, SSO_LCI[i])
        push!(A_SSO_WTE_vec, A_SSO_WTE[i])
        push!(A_SSO_LF_vec,  A_SSO_LF[i])
        push!(AE_vec, AE[i])
        push!(AH_vec, AH[i])
        push!(CH4_vec, CH4[i])
        push!(CHP_vec, CHP[i])
        push!(ferric_vec, ferric[i])
        push!(ng_LCI_vec, ng_LCI[i])
        push!(building_vec,  building[i])
        push!(concrete_vec,  concrete[i])
        push!(earth_vec,     earth[i])
        push!(gravel_vec,    gravel[i])
        push!(steel_vec,     steel[i])
        push!(elec_vec,      elec[i])
    end

    return index, times_vec, SSO_vec, AE_vec, AH_vec, CH4_vec, CHP_vec,
           ferric_vec, ng_LCI_vec, A_SSO_WTE_vec, A_SSO_LF_vec,
           building_vec, concrete_vec, earth_vec, gravel_vec, steel_vec, elec_vec
end

index, times_vec, SSO_vec, AE_vec, AH_vec, CH4_vec, CHP_vec,
ferric_vec, ng_LCI_vec, A_SSO_WTE_vec, A_SSO_LF_vec,
building_vec, concrete_vec, earth_vec, gravel_vec, steel_vec, elec_vec =
    check_data_if_missing()

function max_consecutive_run_with_index(arr)
    if isempty(arr)
        return 0, 0
    end
    max_len = 1
    current_len = 1
    max_start = 1
    current_start = 1

    for i in 2:length(arr)
        if arr[i] == arr[i-1] + 1
            current_len += 1
        else
            current_len = 1
            current_start = i
        end

        if current_len > max_len
            max_len = current_len
            max_start = current_start
        end
    end

    return max_len, max_start
end

function save_cleaned_data_to_csv(index, times_vec, SSO_vec, AE_vec, AH_vec,
                                  CH4_vec, CHP_vec, ferric_vec, ng_LCI_vec,
                                  A_SSO_WTE_vec, A_SSO_LF_vec,
                                  building_vec, concrete_vec, earth_vec,
                                  gravel_vec, steel_vec, elec_vec;
                                  filename="cleaned_data_per_cubic_AD_feed.csv")
    df = DataFrame(
        index      = index,
        time       = times_vec,
        SSO_LCI    = SSO_vec,
        AE         = AE_vec,
        AH         = AH_vec,
        CH4        = CH4_vec,
        CHP        = CHP_vec,
        ferric     = ferric_vec,
        ng_LCI     = ng_LCI_vec,
        A_SSO_WTE  = A_SSO_WTE_vec,
        A_SSO_LF   = A_SSO_LF_vec,
        building   = building_vec,
        concrete   = concrete_vec,
        earth      = earth_vec,
        gravel     = gravel_vec,
        steel      = steel_vec,
        elec       = elec_vec
    )
    CSV.write(filename, df)
    return filename
end

save_cleaned_data_to_csv(index, times_vec, SSO_vec, AE_vec, AH_vec,
                         CH4_vec, CHP_vec, ferric_vec, ng_LCI_vec,
                         A_SSO_WTE_vec, A_SSO_LF_vec,
                         building_vec, concrete_vec, earth_vec,
           gravel_vec, steel_vec, elec_vec)

# Export aligned raw key inputs (including total feed dig_tfeed) for raw → impact mapping
function save_raw_inputs_with_feed(index_vec)
    raw_df = DataFrame(
        index       = index_vec,
        time        = times[index_vec],
        dig_tfeed   = dig_tfeed[index_vec],                 # kgal/day
        dig_PSSO    = dig_PSSO[index_vec],                  # kgal/day (unit corrected)
        SSO_mass_kg = dig_PSSO[index_vec] .* 1000 .* 3.99,  # kg/day (before volume normalization)
        AE_raw_MJ   = gas_energy_to_CHP[index_vec] .* 0.40, # MJ/day avoided electricity
        AH_raw_MJ   = (energy_to_pelletization[index_vec] .+ gas_energy_to_CHP[index_vec] .* 0.39), # MJ/day avoided heat
        ng_energy_MJ = ng_energy_daily[index_vec],          # MJ/day natural gas
        steel_raw_const = fill(4.91e-5 * 23, length(index_vec)) # CAPEX steel numerator constant (flow-independent)
    )
    CSV.write("raw_inputs_with_feed.csv", raw_df)
end

# Export raw inputs using the indices aligned with cleaned_data
save_raw_inputs_with_feed(index)

# Notes: five inputs used by the RSM surrogate and their sources (excluding constants/defaults)
# - SSO_LCI  : from dig_PSSO (digester SSO feed, gal/kgal→Mgal), normalized by dig_tfeed (total AD feed in) to per m³
# - AE       : avoided electricity, computed from biogas out_gas/out_gas2 as CHP power generation
# - AH       : avoided heat, computed from CHP waste heat (proportional to biogas production)
# - ng_LCI   : natural gas consumption, normalized by dig_tfeed to per m³
# - steel    : steel quantity as CAPEX proxy (steel_vec), read from CAPEX breakdown and normalized by dig_tfeed to per m³

##############################
# Temperature data: Excel → daily averages
##############################

function extract_vectors_celsius_safe(file_path::String)
    xf = XLSX.readxlsx(file_path)
    sheet = xf[1]

    times = DateTime[]
    in1_f  = Union{Missing,Float64}[]
    out1_f = Union{Missing,Float64}[]
    in2_f  = Union{Missing,Float64}[]
    out2_f = Union{Missing,Float64}[]
    in4_f  = Union{Missing,Float64}[]
    out4_f = Union{Missing,Float64}[]
    in3_f  = Union{Missing,Float64}[]
    out3_f = Union{Missing,Float64}[]

    for r in XLSX.eachrow(sheet)
        rn = XLSX.row_number(r)
        rn < 4 && continue

        push!(times, DateTime(r["A"], dateformat"m/d/y H:M"))
        push!(in1_f,  r["B"])
        push!(out1_f, r["C"])
        push!(in2_f,  r["D"])
        push!(out2_f, r["E"])
        push!(in4_f,  r["F"])
        push!(out4_f, r["G"])
        push!(in3_f,  r["H"])
        push!(out3_f, r["I"])
    end

    f_to_c(F) = ismissing(F) ? missing : (F - 32) * 5 / 9

    in1_c  = f_to_c.(in1_f)
    out1_c = f_to_c.(out1_f)
    in2_c  = f_to_c.(in2_f)
    out2_c = f_to_c.(out2_f)
    in4_c  = f_to_c.(in4_f)
    out4_c = f_to_c.(out4_f)
    in3_c  = f_to_c.(in3_f)
    out3_c = f_to_c.(out3_f)

    return times,
           in1_c, out1_c,
           in2_c, out2_c,
           in4_c, out4_c,
           in3_c, out3_c
end

# Pull data from two files
t1, in1_1, out1_1, in2_1, out2_1, in4_1, out4_1, in3_1, out3_1 =
    extract_vectors_celsius_safe("Jan-Dec 2021 Digester Temp Data 30min Readings.xlsx")
t2, in1_2, out1_2, in2_2, out2_2, in4_2, out4_2, in3_2, out3_2 =
    extract_vectors_celsius_safe("Jan-July 2022 Digester Temp Data 30min Readings.xlsx")

# Combine vectors
time_all = vcat(t1, t2)
in1_all  = vcat(in1_1, in1_2)
out1_all = vcat(out1_1, out1_2)
in2_all  = vcat(in2_1, in2_2)
out2_all = vcat(out2_1, out2_2)
in4_all  = vcat(in4_1, in4_2)
out4_all = vcat(out4_1, out4_2)
in3_all  = vcat(in3_1, in3_2)
out3_all = vcat(out3_1, out3_2)

function daily_average(times, temps)
    daily = Dict{Date, Vector{Union{Missing,Float64}}}()

    for (t, v) in zip(times, temps)
        d = Date(t)
        push!(get!(daily, d, Union{Missing,Float64}[]), v)
    end

    dates = sort(collect(keys(daily)))
    means = [mean(skipmissing(daily[d])) for d in dates]

    return dates, means
end

dates_temp, in1_daily   = daily_average(time_all, in1_all)
_,          out1_daily  = daily_average(time_all, out1_all)
_,          in2_daily   = daily_average(time_all, in2_all)
_,          out2_daily  = daily_average(time_all, out2_all)
_,          in4_daily   = daily_average(time_all, in4_all)
_,          out4_daily  = daily_average(time_all, out4_all)
_,          in3_daily   = daily_average(time_all, in3_all)
_,          out3_daily  = daily_average(time_all, out3_all)

df_temp = DataFrame(
    Date      = dates_temp,
    In1_mean  = in1_daily,
    Out1_mean = out1_daily,
    In2_mean  = in2_daily,
    Out2_mean = out2_daily,
    In4_mean  = in4_daily,
    Out4_mean = out4_daily,
    In3_mean  = in3_daily,
    Out3_mean = out3_daily
)

CSV.write("daily_avg_temp.csv", df_temp)

##############################
# Correlation analysis
##############################

function collect_non_missing_pairs(vec1, vec2)
    @assert length(vec1) == length(vec2) "Two vectors must have the same length"

    valid_indices = Int[]
    valid_vec1 = Float64[]
    valid_vec2 = Float64[]

    for i in eachindex(vec1)
        if !ismissing(vec1[i]) && !ismissing(vec2[i])
            push!(valid_indices, i)
            push!(valid_vec1, vec1[i])
            push!(valid_vec2, vec2[i])
        end
    end

    println("Raw data length: ", length(vec1))
    println("Count of valid pairs: ", length(valid_vec1))
    println("Share of valid data: ",
            round(length(valid_vec1) / length(vec1) * 100, digits = 2), "%")

    data_matrix = hcat(valid_vec1, valid_vec2)'
    return data_matrix, valid_indices
end

# out_gas vs dig_PSSO
gas_psso_data, valid_indices = collect_non_missing_pairs(out_gas, dig_PSSO)

println("\nCollected data matrix size: ", size(gas_psso_data))
if size(gas_psso_data, 2) >= 2
    println("out_gas valid range: ",
            minimum(gas_psso_data[1, :]), " - ", maximum(gas_psso_data[1, :]))
    println("dig_PSSO valid range: ",
            minimum(gas_psso_data[2, :]), " - ", maximum(gas_psso_data[2, :]))
    correlation = cor(gas_psso_data[1, :], gas_psso_data[2, :])
    println("Correlation between out_gas and dig_PSSO: ",
            round(correlation, digits = 4))
else
    println("Too few valid data points to compute correlation between out_gas and dig_PSSO")
end

println("\n==============================")
println("Analyzing correlation for centrifuge_cakesolids")

# centrifuge_cakesolids vs out_gas
gas_cake_data, gas_cake_indices = collect_non_missing_pairs(out_gas, centrifuge_cakesolids)
if size(gas_cake_data, 2) >= 2
    correlation_gas_cake = cor(gas_cake_data[1, :], gas_cake_data[2, :])
    println("Correlation between out_gas and centrifuge_cakesolids: ",
            round(correlation_gas_cake, digits = 4))
else
    println("Too few valid data points to compute correlation between out_gas and centrifuge_cakesolids")
end

# centrifuge_cakesolids vs dig_PSSO
psso_cake_data, psso_cake_indices =
    collect_non_missing_pairs(dig_PSSO, centrifuge_cakesolids)
if size(psso_cake_data, 2) >= 2
    correlation_psso_cake = cor(psso_cake_data[1, :], psso_cake_data[2, :])
    println("Correlation between dig_PSSO and centrifuge_cakesolids: ",
            round(correlation_psso_cake, digits = 4))
else
    println("Too few valid data points to compute correlation between dig_PSSO and centrifuge_cakesolids")
end

# Three variables simultaneously valid
function collect_three_variable_pairs(vec1, vec2, vec3)
    @assert length(vec1) == length(vec2) &&
            length(vec2) == length(vec3) "All three vectors must have the same length"

    valid_indices = Int[]
    valid_vec1 = Float64[]
    valid_vec2 = Float64[]
    valid_vec3 = Float64[]

    for i in eachindex(vec1)
        if !ismissing(vec1[i]) && !ismissing(vec2[i]) && !ismissing(vec3[i])
            push!(valid_indices, i)
            push!(valid_vec1, vec1[i])
            push!(valid_vec2, vec2[i])
            push!(valid_vec3, vec3[i])
        end
    end

    println("Count of triplets with all valid data: ", length(valid_vec1))
    println("Share of triplets with valid data: ",
            round(length(valid_vec1) / length(vec1) * 100, digits = 2), "%")

    data_matrix = hcat(valid_vec1, valid_vec2, valid_vec3)'
    return data_matrix, valid_indices
end

three_var_data, three_var_indices =
    collect_three_variable_pairs(out_gas, dig_PSSO, centrifuge_cakesolids)

if size(three_var_data, 2) > 0
    println("\nThree-variable correlation matrix:")
    corr_matrix = cor(three_var_data')
    var_names = ["out_gas", "dig_PSSO", "centrifuge_cakesolids"]

    for i in 1:3
        for j in 1:3
            println("$(var_names[i]) vs $(var_names[j]): ",
                    round(corr_matrix[i, j], digits = 4))
        end
    end
else
    println("No data points found where all three variables are valid simultaneously")
end

# Combined plot: Cake + total gas + PSSO
a16 = plot(times, centrifuge_cakesolids,
           label = "Centrifuge Cake Solids (Short Dry Tons/day)",
           xlabel = "Time", legend = :topleft)
plot!(a16, times, out_gas ./ 100,
      label = "Total Gas Production (100 k ft^3/day)")
plot!(a16, times, dig_PSSO,
      label = "Digester PSSO (k gal/day)",
      ylabel = "y axis",
      title = "Digester PSSO / Total Biogas / Centrifuge Cake Solids Over Time")


include("compute_sankey_numbers.jl")

# Default: average across all days
res = compute_sankey_numbers()

# Quick look at percentages
res.mass_pct
res.gas_mass_split_pct
res.energy_pct

# Print SankeyMATIC-ready text
print_sankeymatic_text(res)
