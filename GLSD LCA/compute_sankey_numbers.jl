# compute_sankey_numbers.jl
#
# Usage:
# 1. Run LCI_Calculation_FU_PerCubicMeter_AD_feedin V2.jl first
#    and ensure these variables exist in scope: times, den_safe, SSO_LCI,
#    out_gas, gas_to_pelletization, energy_to_pelletization,
#    gas_energy_to_CHP, ng_energy_daily, centrifuge_cakesolids, ch4_frac
# 2. Then include this file and call compute_sankey_numbers().

using Statistics

"""
Compute the numbers needed for mass/energy Sankey diagrams.

Arguments:
    idx  — indices of the days to average (defaults to all good days)

Returns:
    A NamedTuple with:
      :mass_abs            absolute mass of each stream (kg/day mean)
      :mass_pct            mass percentages scaled so source/stream each sum to 100
      :gas_mass_split_pct  biogas volume split to dryer/CHP/fugitive (scaled)
      :energy_abs          absolute energy of each stream (MJ/day mean)
      :energy_pct          energy percentages scaled to 100
"""
function compute_sankey_numbers(idx = eachindex(times))

    # Helper: mean over idx while skipping missings
    m(x) = mean(skipmissing(x[idx]))

    # ----------------------
    # 1. Mass basis (kg/day)
    # ----------------------

    # SSO mass: SSO_LCI (kg/m^3) * den_safe (m^3/day)
    m_SSO_day = SSO_LCI .* den_safe     # kg/day

    # Cake dry solids: short ton/day -> kg/day
    cake_kg_day = centrifuge_cakesolids .* 907.185

    # Biogas mass: compute total volume, then convert using CH4/CO2 densities
    V_bg_m3_day = out_gas .* 0.0283168 .* 1000.0  # kft^3 -> ft^3 -> m^3
    ρ_CH4 = 0.656    # kg/m^3
    ρ_CO2 = 1.977    # kg/m^3 (approx.)
    m_CH4_day = V_bg_m3_day .* ch4_frac .* ρ_CH4
    m_CO2_day = V_bg_m3_day .* (1 .- ch4_frac) .* ρ_CO2
    m_bg_day  = m_CH4_day .+ m_CO2_day          # kg/day

    # Approximate total digester outflow as cake + biogas
    # and assume inflow mass ≈ outflow mass
    m_feed_day = cake_kg_day .+ m_bg_day

    # SSO daily mass cannot exceed total mass; clip to avoid negative sludge
    m_SSO_day_clipped = min.(m_SSO_day, m_feed_day)
    m_sludge_day = m_feed_day .- m_SSO_day_clipped

    mass_abs = (
        sludge_in  = m(m_sludge_day),
        SSO_in     = m(m_SSO_day_clipped),
        digestate  = m(cake_kg_day),
        biogas     = m(m_bg_day),
    )

    # Normalize sources and sinks separately to sum to 100 for Sankey plots
    total_in_mass  = mass_abs.sludge_in + mass_abs.SSO_in
    total_out_mass = mass_abs.digestate + mass_abs.biogas

    mass_in_scale  = 100.0 / total_in_mass
    mass_out_scale = 100.0 / total_out_mass

    mass_pct = (
        sludge_in  = mass_abs.sludge_in  * mass_in_scale,
        SSO_in     = mass_abs.SSO_in     * mass_in_scale,
        digestate  = mass_abs.digestate  * mass_out_scale,
        biogas     = mass_abs.biogas     * mass_out_scale,
    )

    # ---- Further split biogas stream by volume into dryer / CHP / fugitive ----

    V_bg_ft3_day = out_gas .* 1000.0
    V_dryer_ft3  = gas_to_pelletization
    V_CHP_ft3    = V_bg_ft3_day .- V_dryer_ft3

    # Assume ~5% of biogas is fugitive/combustion loss (matches CH4 slip_frac)
    fugitive_frac   = 0.05
    V_fugitive_ft3  = V_bg_ft3_day .* fugitive_frac
    V_used_ft3      = V_bg_ft3_day .- V_fugitive_ft3

    # Normalize dryer and CHP volumes to V_used (avoid dryer+CHP > total)
    denom = V_dryer_ft3 .+ V_CHP_ft3
    adj = similar(denom, Float64)
    for i in eachindex(denom)
        if ismissing(denom[i]) || denom[i] <= 0 || ismissing(V_used_ft3[i])
            adj[i] = 1.0
        else
            adj[i] = V_used_ft3[i] / denom[i]
        end
    end

    V_dryer_ft3_adj = V_dryer_ft3 .* adj
    V_CHP_ft3_adj   = V_CHP_ft3   .* adj

    frac_dryer = m(V_dryer_ft3_adj ./ V_bg_ft3_day)
    frac_CHP   = m(V_CHP_ft3_adj   ./ V_bg_ft3_day)
    frac_fug   = m(V_fugitive_ft3  ./ V_bg_ft3_day)

    m_gas_total = mass_abs.biogas

    gas_mass_split_pct = (
        biogas_to_dryer  = m_gas_total * frac_dryer * mass_out_scale,
        biogas_to_CHP    = m_gas_total * frac_CHP   * mass_out_scale,
        biogas_fugitive  = m_gas_total * frac_fug   * mass_out_scale,
    )

    # ----------------------
    # 2. Energy basis (MJ/day)
    # ----------------------

    # Total biogas chemical energy (same 550/948 MJ/ft^3 factor as main script)
    LHV_bg_MJ_per_ft3 = 550.0 / 948.0
    E_bg_total_day = V_bg_ft3_day .* LHV_bg_MJ_per_ft3

    E_dryer_bg_day = energy_to_pelletization           # MJ/day, directly from the script
    E_CHP_bg_day   = gas_energy_to_CHP                 # MJ/day, CHP inlet chemical energy
    E_NG_day       = ng_energy_daily                   # MJ/day, natural gas chemical energy

    # Treat the remainder as fugitives + other combustion losses
    E_fugitive_day = E_bg_total_day .- (E_dryer_bg_day .+ E_CHP_bg_day)
    E_fugitive_day_clamped = max.(E_fugitive_day, 0.0)

    energy_abs = (
        biogas_total = m(E_bg_total_day),
        dryer_bg     = m(E_dryer_bg_day),
        CHP_bg       = m(E_CHP_bg_day),
        fugitive     = m(E_fugitive_day_clamped),
        NG           = m(E_NG_day),
    )

    total_energy_in = energy_abs.biogas_total + energy_abs.NG
    energy_scale = 100.0 / total_energy_in

    energy_pct = (
        biogas_total = energy_abs.biogas_total * energy_scale,
        NG           = energy_abs.NG           * energy_scale,
        dryer_bg     = energy_abs.dryer_bg    * energy_scale,
        CHP_bg       = energy_abs.CHP_bg      * energy_scale,
        fugitive     = energy_abs.fugitive    * energy_scale,
    )

    return (
        mass_abs          = mass_abs,
        mass_pct          = mass_pct,
        gas_mass_split_pct = gas_mass_split_pct,
        energy_abs        = energy_abs,
        energy_pct        = energy_pct,
    )
end

# Small helper: print results as SankeyMATIC-ready text
function print_sankeymatic_text(res)
    mp = res.mass_pct
    gp = res.gas_mass_split_pct
    ep = res.energy_pct

    println("//// Mass Sankey (numbers already normalized to ~100) ////")
    println("Wastewater sludge (primary + secondary) [$(round(mp.sludge_in, digits=1))] Anaerobic digester")
    println("Source-separated organics (SSO) [$(round(mp.SSO_in, digits=1))] Anaerobic digester")
    println("Anaerobic digester [$(round(mp.digestate, digits=1))] Digestate solids")
    println("Anaerobic digester [$(round(mp.biogas, digits=1))] Biogas (CH4 + CO2)")
    println("Biogas (CH4 + CO2) [$(round(gp.biogas_to_dryer, digits=1))] Pelletization dryer (biogas for heat)")
    println("Biogas (CH4 + CO2) [$(round(gp.biogas_to_CHP, digits=1))] CHP (biogas for electricity + heat)")
    println("Biogas (CH4 + CO2) [$(round(gp.biogas_fugitive, digits=1))] Fugitive + combustion emissions")

    println("\n//// Energy Sankey (numbers normalized to ~100) ////")
    println("Biogas chemical energy [$(round(ep.biogas_total, digits=1))] Anaerobic digester energy pool")
    println("Natural gas chemical energy [$(round(ep.NG, digits=1))] Anaerobic digester energy pool")
    println("Anaerobic digester energy pool [$(round(ep.dryer_bg, digits=1))] Pelletization dryer (useful heat)")
    println("Anaerobic digester energy pool [$(round(ep.CHP_bg, digits=1))] CHP (useful electricity + heat)")
    println("Anaerobic digester energy pool [$(round(ep.fugitive, digits=1))] Fugitive + other losses")
end
