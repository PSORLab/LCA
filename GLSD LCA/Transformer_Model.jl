

# This code converts the transformer model into Julia.

# Define some structures for labeling purposes
mutable struct TransformerInput
    date::Vector{String}   # Date
    CODp::Vector{Union{Missing, Real}} # Particulate chemical oxygen demand (COD), gCOD/m^3 [Total COD - Soluble COD]
    CODs_VFA::Vector{Union{Missing, Real}} # Soluble COD minus VFAs, gCOD/m^3
    VFA::Vector{Union{Missing, Real}} # Volatile fatty acid concentration, gCOD/m^3
    TOC::Vector{Union{Missing, Real}} # Total organic carbon, gC/m^3 [Total carbon (TC) - total inorganic carbon (TIC)]
    N_org::Vector{Union{Missing, Real}} # Organic nitrogen, gN/m^3 [Total Kjeldahl nitrogen (TKN) - total ammonia-nitrogen (TAN)]
    TAN::Vector{Union{Missing, Real}} # Total ammonia-nitrogen, gN/m^3 [TKN - N_org]
    P_org::Vector{Union{Missing, Real}} # Organic phosphorous, gP/m^3 [Total phosphorous - Orthophosphate (orthoP)]
    orthoP::Vector{Union{Missing, Real}} # Orthophosphate, gP/m^3 [Total phosphorous - organic phosphorous]
    TIC::Vector{Union{Missing, Real}} # Total inorganic carbon, mol HCO3-/m^3 [mainly bicarbonate, can be estimated from titration of alkalinity] [(TC - TOC)/12, from Nguyen2014]
    S_cat::Vector{Union{Missing, Real}} # Total alkalinity, equ/m^3 [see Bernard2001]
    FS::Vector{Union{Missing, Real}} # Fixed solids, g/m^3 [Total solids (TS) - total volatile solids (TVS)]
    Q::Vector{Union{Missing, Real}} # Flow rate, m^3/day
end

mutable struct ADM1Input #Units should be kgCOD/m^3, but check that Zaher didn't make it gCOD/m^3
    date::Vector{String}   # Date
    S_su::Vector{Union{Missing, Real}}  # Monosaccharides, kgCOD/m^3
    S_aa::Vector{Union{Missing, Real}}  # Amino acids
    S_fa::Vector{Union{Missing, Real}}  # Total LCFA
    S_va::Vector{Union{Missing, Real}}  # Total valerate
    S_bu::Vector{Union{Missing, Real}}  # Total butyrate
    S_pro::Vector{Union{Missing, Real}} # Total propionate
    S_ac::Vector{Union{Missing, Real}}  # Total acetate, kgCOD/m^3
    S_h2::Vector{Union{Missing, Real}}  # Hydrogen
    S_ch4::Vector{Union{Missing, Real}} # Methane
    S_IC::Vector{Union{Missing, Real}}  # Inorganic carbon, kmoleC/m^3
    S_IN::Vector{Union{Missing, Real}}  # Inorganic nitrogen, kmoleN/m^3
    S_I::Vector{Union{Missing, Real}}   # Soluble inerts
    X_c::Vector{Union{Missing, Real}}   # Composite
    X_ch::Vector{Union{Missing, Real}}  # Carbohydrates, kgCOD/m^3
    X_pr::Vector{Union{Missing, Real}}  # Proteins, kgCOD/m^3
    X_li::Vector{Union{Missing, Real}}  # Lipids, kgCOD/m^3
    X_su::Vector{Union{Missing, Real}}  # Biomass (monosaccharides)
    X_aa::Vector{Union{Missing, Real}}  # Biomass (amino acids)
    X_fa::Vector{Union{Missing, Real}}  # Biomass (LCFA)
    X_c4::Vector{Union{Missing, Real}} # Biomass (methane)
    X_pro::Vector{Union{Missing, Real}} # Biomass (propionate)
    X_ac::Vector{Union{Missing, Real}}  # Biomass (acetate)
    X_h2::Vector{Union{Missing, Real}}  # Biomass (hydrogen)
    X_I::Vector{Union{Missing, Real}}   # Particulate inerts, kgCOD/m^3
    S_cat::Vector{Union{Missing, Real}} # Cations, kmole/m^3
    S_an::Vector{Union{Missing, Real}}  # Anions, kmole/m^3
    S_IP::Vector{Union{Missing, Real}}  # Inorganic phosphorous, kmoleP/m^3 (NOT USED IN SAM'S MODEL??)
    Q::Vector{Union{Missing, Real}}     # Flow rate, m^3/day
    # What about S_OH- and S_H+, from the transformer model table?
end

# Define the transform matrices (from the InputData.xls file)
TransStoic = [0	            0	0	0	            0	            -1	0	            0	0	0	0	0	    0	0	0	0	0	0	    0	0	0.00E+00	7.14E-05	0	0	0	    0	        0	        0	0	0	0	0	0	0	0	        0	    -4.16E-20	0.00E+00
              0	            0	0	0	            0	            0	0	            0	-1	0	0	0	    0	0	0	0	0	0	    0	0	1.00E-03	0.00E+00	0	0	0	    0	        0	        0	0	0	0	0	0	0	0	        0	    0.00E+00	0.00E+00
              0	            0	0	0	            0	            0	0	            -1	0	0	0	0	    0	0	0	0	0	0	    0	0	0.00E+00	0.00E+00	0	0	0	    0	        0	        0	0	0	0	0	0	0	0	        0	    5.43E-07	3.17E-05
              0	            0	0	0	            0	            0	0	            0	0	-1	0	0	    0	0	0	0	0	0	    0	0	0.00E+00	0.00E+00	0	0	0	    0	        0	        0	0	0	0	0	0	0	0	        0.001	0.00E+00	0.00E+00
              0	            0	-1	-0.3750012	    0	            0	0	            0	0	0	0	0	    0	0	0	0	0	0.001	0	0	1.00E-10	0.00E+00	0	0	0	    0	        0	        0	0	0	0	0	0	0	0	        0	    -5.00E-11	0.00E+00
              0	            -1	0	-0.375	        0	            0	0	            0	0	0	0	0.001	0	0	0	0	0	0	    0	0	-1.07E-18	0.00E+00	0	0	0	    0	        0	        0	0	0	0	0	0	0	0	        0	    5.47E-19	0.00E+00
              -6.457862448	0	0	-2.712302228	0	            0	-1	            0	0	0	0	0	    0	0	0	0	0	0	    0	0	0.00E+00	0.00E+00	0	0	0	    0	        0.006457862	0	0	0	0	0	0	0	0	        0	    -3.83E-05	0.00E+00
              -6.857142857	0	0	-2.571429771	-1	            0	0	            0	0	0	0	0	    0	0	0	0	0	0	    0	0	1.00E-10	0.00E+00	0	0	0	    0.006857143	0	        0	0	0	0	0	0	0	0	        0	    3.57E-05	0.00E+00
              -1	        0	0	-0.375	        0	            0	0	            0	0	0	0	0	    0	0	0	0	0	0	    0	0	1.07E-18	0.00E+00	0	0	0.001	0	        0	        0	0	0	0	0	0	0	0	        0	    -5.37E-19	0.00E+00
              -0.990616184	0	0	-0.304630328	-0.055920989	0	-0.006387546	0	0	0	-1	0	    0	0	0	0	0	0	    0	0	3.71E-07	5.33E-09	0	0	0	    0	        0	        0	0	0	0	0	0	0	0.000990616	0	    4.30E-07	3.75E-12]

# Adjusted transform matrix to match the Zaher2009 paper:
# TransStoic = [0	            0	0	0	            0	            -1	0	            0	0	0	0	0	    0	0	0	0	0	0	    0	0	0.00E+00	7.14E-05	0	0	0	    0	        0	        0	0	0	0	0	0	0	0	        0	    -4.16E-20	0.00E+00
#               0	            0	0	0	            0	            0	0	            0	-1	0	0	0	    0	0	0	0	0	0	    0	0	1.00E-03	0.00E+00	0	0	0	    0	        0	        0	0	0	0	0	0	0	0	        0	    0.00E+00	0.00E+00
#               0	            0	0	0	            0	            0	0	            -1	0	0	0	0	    0	0	0	0	0	0	    0	0	0.00E+00	0.00E+00	0	0	0	    0	        0	        0	0	0	0	0	0	0	0	        0	    5.43E-07	3.17E-05
#               0	            0	0	0	            0	            0	0	            0	0	-1	0	0	    0	0	0	0	0	0	    0	0	0.00E+00	0.00E+00	0	0	0	    0	        0	        0	0	0	0	0	0	0	0	        0.001	0.00E+00	0.00E+00
#               0	            0	-1	-0.3750012	    0	            0	0	            0	0	0	0	0	    0	0	0	0	0	0.001	0	0	1.00E-10	0.00E+00	0	0	0	    0	        0	        0	0	0	0	0	0	0	0	        0	    -5.00E-11	0.00E+00
#               0	            -1	0	-0.375	        0	            0	0	            0	0	0	0	0.001	0	0	0	0	0	0	    0	0	-1.07E-18	0.00E+00	0	0	0	    0	        0	        0	0	0	0	0	0	0	0	        0	    5.47E-19	0.00E+00
#               -6.457862448	0	0	-2.712302228	0	            0	-1	            0	0	0	0	0	    0	0	0	0	0	0	    0	0	0.00E+00	0.00E+00	0	0	0	    0	        0.006457862	0	0	0	0	0	0	0	0	        0	    -3.83E-05	0.00E+00
#               -6.857142857	0	0	-2.571429771	-1	            0	0	            0	0	0	0	0	    0	0	0	0	0	0	    0	0	1.00E-10	0.00E+00	0	0	0	    0.006857143	0	        0	0	0	0	0	0	0	0	        0	    3.57E-05	0.00E+00
#               -1	        0	0	-0.375	        0	            0	0	            0	0	0	0	0	    0	0	0	0	0	0	    0	0	1.07E-18	0.00E+00	0	0	0.001	0	        0	        0	0	0	0	0	0	0	0	        0	    -5.37E-19	0.00E+00
#               -0.990616184	0	0	-0.404630328	-0.055920989	0	-0.006387546	0	0	0	-1	0	    0	0	0	0	0	0	    0	0	2.85E-07    -1.09E-07	0	0	0	    0	        0	        0	0	0	0	0	0	0	0.000990616	0	    8.18E-07	3.75E-12]

MaxOrder = [9 4 5 6 7 8 3 2 0 1] .+ 1 # Convertion from 0-index to 1-index


# Inputs will be a 12-length vector, OR a [n][12] vector of vectors, OR an [n,12] matrix
function transformer(inputs::Vector{T}) where T <: Real
    # This function does the transformations in the Zaher2007 paper

    # Preallocate outputs and temporary variables
    output = zeros(T, 28) # 27 state variables + 1 flow variable
    flux_in = zeros(T, 11)
    flux_out = zeros(T, 27)
    rho = zeros(T, 10)
    temp = zeros(T, 11)

    # Check the size of inputs
    if size(inputs)!=(12,) # 11 measurements + 1 flow variable
        error("Incorrect input dimensions.")
    end

    # Return a vector of missings if any value is missing
    if count(ismissing.(inputs)) > 0
        out = zeros(Union{Missing,Real}, 28)
        out .= missing
        return out
    end

    # Fill in the flux_in vector
    flux_in = inputs[1:11] .* inputs[12]

    # Step 1) Equation 4: Allocate inlet flux
    for i = 1:10
        Z = MaxOrder[i]
        K = findfirst(x -> x == -1.0, TransStoic[Z,1:11])
        rho[Z] = flux_in[K] / TransStoic[Z,K] # Dividing by -1
        for j = 1:(i-1)
            z = MaxOrder[j]
            rho[Z] -= rho[z] * TransStoic[z,K]/TransStoic[Z,K] # Remove rates that went to other species
        end

        # Step 2 and 3) Equation 5: Make sure we haven't over-allocated inlet flux
        J = 0
        for k = 1:11
            temp[k] = 0.0
            for j = 1:i
                z = MaxOrder[j]
                if (temp[k] + rho[z]*TransStoic[z,k]) <= flux_in[k]
                    temp[k] += rho[z]*TransStoic[z,k]
                elseif (J==0) && (TransStoic[Z,k] != 0.0)
                    J = k
                    rho[Z] = (flux_in[J] - temp[J]) / TransStoic[Z,J]
                elseif (J!=0) && (TransStoic[Z,k] != 0.0) && ((flux_in[J] - temp[J])/TransStoic[Z,J] < ((flux_in[k] - temp[k])/TransStoic[Z,k]))
                    rho[Z] = (flux_in[k] - temp[k])/TransStoic[Z,k]
                end
            end
        end
    end

    # Step 4) Any remaining influxes are assigned to inorganic components
    rho[3] -= (flux_in[7] - temp[7]) # P
    rho[1] -= (flux_in[5] - temp[5]) # N
    rho[2] -= (flux_in[4] - temp[4]) # C

    # Additional COD is assigned to "S_I" (Soluble Inerts) since "X_I" (Particulate Inerts) was estimated from Fixed Solids (?)
    flux_out[12] = (flux_in[1] - temp[1])/1000

    # Calculation of outfluxes (i.e., ADM1 influxes)
    # Solubles
    for i = 1:11
        k = i+11
        for j = 1:10
            flux_out[i] -= rho[j]*TransStoic[j,k]
        end
        if flux_out[i] < 1.0E-12
            flux_out[i] = 0.0
        end
    end
    # Particulates, S_cat, S_an, and SIP
    for i = 13:27
        k = i+11
        for j = 1:10
            flux_out[i] -= rho[j]*TransStoic[j,k]
        end
        if flux_out[i] < 1.0E-12
            flux_out[i] = 0.0
        end
    end
    # Outputs?
    output[28] = inputs[12] # Flow rate carries through
    output[1:27] .= flux_out ./ output[28]

    return output
end

transformer(inputs::Vector{Vector{T}}) where T <: Real = transformer.(inputs)
function transformer(inputs...)
    corrected_input = hcat(inputs...)
    return transformer(corrected_input)
end
function transformer(inputs::TransformerInput)
    output = transformer(inputs.CODp, inputs.CODs_VFA, inputs.VFA, inputs.TOC, inputs.N_org,
                         inputs.TAN, inputs.P_org, inputs.orthoP, inputs.TIC, inputs.S_cat,
                         inputs.FS, inputs.Q)
    return ADM1Input(inputs.date, [output[:,i] for i in 1:28]...)
end



# Inputs will be an [n,12] matrix
function transformer(inputs::Matrix)
    # This function does the transformations in the Zaher2007 paper

    in_size = size(inputs)

    # Preallocate outputs and temporary variables
    output = zeros(Union{Real,Missing}, in_size[1], 28) # 27 state variables + 1 flow variable
    flux_in = zeros(Union{Real,Missing}, in_size[1], 11)
    flux_out = zeros(Union{Real,Missing}, in_size[1], 27)
    rho = zeros(Union{Real,Missing}, in_size[1], 10)
    temp = zeros(Union{Real,Missing}, in_size[1], 11)
    missingflag = zeros(Bool, in_size[1])

    # Check the size of inputs
    if in_size[2]!=12 # 11 measurements + 1 flow variable
        error("Incorrect input dimensions.")
    end

    # Return a vector of missings if any value is missing
    for i in 1:in_size[1]
        if count(ismissing.(inputs[i,:]))>0
            missingflag[i] = true
            inputs[i,:] .= 0.0
        end
    end

    # Fill in the flux_in vector
    flux_in .= inputs[:,1:11] .* inputs[:,12]

    # Step 1) Equation 4: Allocate inlet flux
    for i = 1:10
        Z = MaxOrder[i]
        K = findfirst(x -> x == -1.0, TransStoic[Z,1:11])
        rho[:,Z] .= flux_in[:,K] / TransStoic[Z,K] # Dividing by -1
        for j = 1:(i-1)
            z = MaxOrder[j]
            rho[:,Z] .-= rho[:,z] .* TransStoic[z,K]/TransStoic[Z,K] # Remove rates that went to other species
        end

        # Step 2 and 3) Equation 5: Make sure we haven't over-allocated inlet flux
        J = zeros(Int, in_size[1])
        for k = 1:11
            temp[:,k] .= 0.0
            for j = 1:i
                z = MaxOrder[j]
                first_adjust = ((temp[:,k] .+ rho[:,z].*TransStoic[z,k]) .<= flux_in[:,k])
                temp[first_adjust,k] .+= rho[first_adjust,z].*TransStoic[z,k]
                
                second_adjust = ((first_adjust.==false) .&& (J[:].==0) .&& (TransStoic[Z,k] != 0.0))
                J[second_adjust] .= k
                for m in 1:in_size[1]
                    if second_adjust[m]
                        rho[m,Z] = (flux_in[m,J[m]] - temp[m,J[m]]) / TransStoic[Z,J[m]]
                    end
                end
                
                for m in 1:in_size[1]
                    if first_adjust[m]==false && second_adjust[m]==false
                        if (J[m]!=0) && (TransStoic[Z,k] != 0.0) && ((flux_in[m,J[m]] - temp[m,J[m]])/TransStoic[Z,J[m]] < ((flux_in[m,k] - temp[m,k])/TransStoic[Z,k]))
                            rho[m,Z] = (flux_in[m,k] - temp[m,k])/TransStoic[Z,k]
                        end
                    end
                end
            end
        end
    end

    # Step 4) Any remaining influxes are assigned to inorganic components
    rho[:,3] .-= (flux_in[:,7] .- temp[:,7]) # P
    rho[:,1] .-= (flux_in[:,5] .- temp[:,5]) # N
    rho[:,2] .-= (flux_in[:,4] .- temp[:,4]) # C

    # Additional COD is assigned to "S_I" (Soluble Inerts) since "X_I" (Particulate Inerts) was estimated from Fixed Solids (?)
    flux_out[:,12] .= (flux_in[:,1] .- temp[:,1])./1000

    # Calculation of outfluxes (i.e., ADM1 influxes)
    # Solubles
    for i = 1:11
        k = i+11
        for j = 1:10
            flux_out[:,i] .-= rho[:,j].*TransStoic[j,k]
        end
        flux_out[flux_out.<1.0E-12] .= 0.0
    end
    # Particulates, S_cat, S_an, and SIP
    for i = 13:27
        k = i+11
        for j = 1:10
            flux_out[:,i] .-= rho[:,j].*TransStoic[j,k]
        end
        flux_out[flux_out.<1.0E-12] .= 0.0
    end
    # Outputs?
    output[:,28] .= inputs[:,12] # Flow rate carries through
    output[:, 1:27] .= flux_out ./ output[:,28]


    # Change NaN's to missing
    output[isnan.(output)] .= missing

    # Re-adjust missings that were inputs to missings in the outputs
    output[missingflag.==true,:] .= missing

    return output
end

#####################################################################
############################ USE EXAMPLES ###########################
#####################################################################

# # # Different allowable inputs
# test_mat = [23550	2520.666667	1146	9339.6	598	284	187	32	60.4	60	5397.00	0.000475;
# 23550	2520.666667	1146	9339.6	598	284	187	32	60.4	60	5397.00	0.000475;
# 168400	5000	4146	45000	4077	716	261.1914894	482.8085106	552	300	31000	1.00E-09;
# 109664.8587	5000	4146	45000	4077	716	261.1914894	482.8085106	552	300	31000	1.00E-09]

test_vecs = [[23550,	2520.666667,	1146,	9339.6,	598,	284,	187,	32,	60.4,	60,	5397.00,	0.000475],
[23550,	2520.666667,	1146,	9339.6,	598,	284,	187,	32,	60.4,	60,	5397.00,	0.000475],
[168400,	5000,	4146,	45000,	4077,	716,	261.1914894,	482.8085106,	552,	300,	31000,	1.00E-09],
[109664.8587,	5000,	4146,	45000,	4077,	716,	261.1914894,	482.8085106,	552,	300,	31000,	1.00E-09],
[307122,	3133.9,	3133.9,	101125.0,	3190.8,	289.2,	382.9,	27.1,	172.0,	25.0,	16000.0,	297.04]]

# test_cols = [["Day one", "Day Two", "Day Three", "Day Four"],
# [23550,	23550,	168400,	109664.8587],
# [2520.666667,	2520.666667,	5000,	5000],
# [1146,	1146,	4146,	4146],
# [9339.6,	9339.6,	45000,	45000],
# [598,	598,	4077,	4077],
# [284,	284,	716,	716],
# [187,	187,	261.1914894,	261.1914894],
# [32,	32,	482.8085106,	482.8085106],
# [60.4,	60.4,	552,	552],
# [60,	60,	300,	300],
# [5397.00,	5397.00,	31000,	31000],
# [0.000475,	0.000475,	1.00E-09,	1.00E-09]]

# test_struct = TransformerInput(test_cols...)

# # Different calls to transformer
# transformer(test_mat)
transformer(test_vecs)
# # transformer(test_cols...)
# new = transformer(test_struct)

# Mayowa used inputs of:
new_test = [307122,	3133.9,	3133.9,	101125.0,	3190.8,	289.2,	382.9,	27.1,	172.0,	25.0,	16000.0,	297.04]

# And got outputs of: 
[3.1339,	0,	0,	0,	0,	0,	3.1339,	0,	0,	87.41291,	0.020742,	41.08838,	0,	232.6266,	15.74444,	1.812717,	0,	0,	0,	0,	0,	0,	0,	15.84986,	0.025,	0.078135,	0.00086,	297.04]
# My version has S_IC being 0.1779365, NOT 87.41291.
# Nguyen has S_IC being 0.0919.
# Nguyen has 38.3668 for S_I, Mayowa and I both have 41.08838