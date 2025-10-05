# using LinearAlgebra
# using Plots
# using SpecialFunctions
# plotly()

# function CA(S, α, β, σ, T, Tw, lx, ly, lz)
#     """
#     Calculate the aerodynamic coefficient CA for a blunt body.
#     # Arguments
#     - 'S': Surface area, m²
#     - 'α': Angle of attack, rad
#     - 'β': Angle of sideslip, rad
#     - 'σ': Reflection coefficient
#     - 'T': Temperature, K
#     - 'Tw': Wall temperature, K
#     - 'lx': Length in x-direction, m
#     - 'ly': Length in y-direction, m
#     - 'lz': Length in z-direction, m
#     # Returns
#     - 'CA': Aerodynamic coefficient CA
#     """
#     σN = σ
#     σT = σ

#     cosα = cos(α)
#     cosβ = cos(β)
#     sinβ = sin(β)
#     sinα = sin(α)
#     CA = ((2-σN)/(S√(π))*cosα*cosβ+sign(cosα*cosβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*cosα^2*cosβ^2) +
#             (2-σN)*(cosα^2*cosβ^2+1/(2*S^2)) * (sign(cosα*cosβ)+erf(S*cosα*cosβ)) + 
#             (σN/(2*S)*cosα*cosβ*√(π*Tw/T)) * (1+sign(cosα*cosβ)*erf(S*cosα*cosβ)) +
#             σT*cosα*cosβ*(lx/ly*(1/(S*√(π))*exp(-S^2*sinβ^2)+sinβ*(sign(sinβ)+erf(S*sinβ))) +
#             lx/lz*(1/(S*√(π))*exp(-S^2*sinα^2*cosβ^2)+sinα*cosβ*(sign(sinα*cosβ)+erf(S*sinα*cosβ))))
#     return CA
# end

# function CS(S, α, β, σ, T, Tw, lx, ly, lz)
#     """
#     Calculate the aerodynamic coefficient CN for a blunt body.
#     # Arguments
#     - 'S': Surface area, m²
#     - 'α': Angle of attack, rad
#     - 'β': Angle of sideslip, rad
#     - 'σ': Reflection coefficient
#     - 'T': Temperature, K
#     - 'Tw': Wall temperature, K
#     - 'lx': Length in x-direction, m
#     - 'ly': Length in y-direction, m
#     - 'lz': Length in z-direction, m
#     # Returns
#     - 'CN': Aerodynamic coefficient CN
#     """
#     σN = σ
#     σT = σ

#     cosα = cos(α)
#     cosβ = cos(β)
#     sinβ = sin(β)
#     sinα = sin(α)
    
#     CS = lx/ly*(((2-σN)/(S*√(π))*sinβ+sign(sinβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*sinβ^2) +
#             (2-σN)*(sinβ^2+1/(2*S^2)) * (sign(sinβ)+erf(S*sinβ)) + 
#             (σN/(2*S)*sinβ*√(π*Tw/T)) * (1+sign(sinβ)*erf(S*sinβ))) +
#             σT*sinβ*(1/(S*√(π))*exp(-S^2*cosα^2*cosβ^2) + cosα*cosβ*(sign(cosα*cosβ)+erf(S*cosα*cosβ)) + 
#             lx/lz*(1/(S*√(π))*exp(-S^2*sinα^2*cosβ^2) + sinα*cosβ*(sign(sinα*cosβ)+erf(S*sinα*cosβ))))
#     return CS
# end

# function CN(S, α, β, σ, T, Tw, lx, ly, lz)
#     """
#     Calculate the aerodynamic coefficient CN for a blunt body.
#     # Arguments
#     - 'S': Surface area, m²
#     - 'α': Angle of attack, rad
#     - 'β': Angle of sideslip, rad
#     - 'σ': Reflection coefficient
#     - 'T': Temperature, K
#     - 'Tw': Wall temperature, K
#     - 'lx': Length in x-direction, m
#     - 'ly': Length in y-direction, m
#     - 'lz': Length in z-direction, m
#     # Returns
#     - 'CN': Aerodynamic coefficient CN
#     """
#     σN = σ
#     σT = σ

#     cosα = cos(α)
#     cosβ = cos(β)
#     sinβ = sin(β)
#     sinα = sin(α)

#     CN = lx/lz*((((2-σN)/(S*√(π))*sinα*cosβ+sign(sinα*cosβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*sinα^2*cosβ^2) +
#             (2-σN)*(sinα^2*cosβ^2+1/(2*S^2)) * (sign(sinα*cosβ)+erf(S*sinα*cosβ)) + 
#             (σN/(2*S)*sinα*cosβ*√(π*Tw/T)) * (1+sign(sinα*cosβ)*erf(S*sinα*cosβ)))) +
#             σT*sinα*cosβ*(lx/ly*(1/(S*√(π))*exp(-S^2*sinβ^2) + sinβ*(sign(sinβ)+erf(S*sinβ))) + 
#             (1/(S*√(π))*exp(-S^2*cosα^2*cosβ^2) + cosα*cosβ*(sign(cosα*cosβ)+erf(S*cosα*cosβ))))
#     return CN
# end

# function normalcoefficient(S, aoa, sigma)
#         CN = 1 / (S^2) * ((((2 - sigma) / sqrt(pi)) * S * sin(aoa) + sigma / 2) * exp(-(S * sin(aoa))^2) + ((2 - sigma) * ((S * sin(aoa))^2 + 0.5) + sigma / 2 * sqrt(pi) * (S * sin(aoa))) * (1 + erf(S * sin(aoa))))
#         return CN
# end

# function axialcoefficient(S, aoa, sigma)
#     CA = ((sigma * cos(aoa)) / (sqrt(pi) * S)) * (exp(-(S * sin(aoa))^2) + sqrt(pi) * (S * sin(aoa)) * (1 + erf(S * sin(aoa))))
#     return CA
# end

# T = 973
# V = 7500
# Tw = 300
# α = deg2rad(0.0)
# β = deg2rad(0.0)
# # σ_list = [0, 0.2, 0.4, 0.6, 0.8, 1]
# σ_list = [0.9]
# lx = 0.0001
# ly = 3.76/2
# lz = 1.93
# gas_constant = 287.05  # Specific gas constant for air, J/(kg·K)
# S = V/√(2*gas_constant*T)

# min_angle = -180
# max_angle = 360
# # min_angle = deg2rad(min_angle)
# # max_angle = deg2rad(max_angle)
# p = plot()
# for σ in σ_list
#     ca_list = Float64[]
#     cn_list = Float64[]
#     cs_list = Float64[]
#     for β = min_angle:max_angle
#         β_rad = deg2rad(β)
#         ca = CA(S, α, β_rad, σ, T, Tw, lx, ly, lz)
#         cn = CN(S, α, β_rad, σ, T, Tw, lx, ly, lz)
#         cs = CS(S, α, β_rad, σ, T, Tw, lx, ly, lz)
#         push!(ca_list, ca)
#         push!(cn_list, cn)
#         push!(cs_list, cs)
#     end
#     plot!(p, [min_angle:max_angle], ca_list, label=" CA", xlabel="Angle of sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Sideslip")
#     plot!(p, [min_angle:max_angle], cn_list, label=" CN", xlabel="Angle of sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Sideslip")
#     plot!(p, [min_angle:max_angle], cs_list, label=" CS", xlabel="Angle of sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Sideslip")
# end
# display(p)

# p = plot()
# for σ in σ_list
#     ca_list = Float64[]
#     cn_list = Float64[]
#     cs_list = Float64[]
#     ca_old_list = Float64[]
#     cn_old_list = Float64[]
#     for α = min_angle:max_angle
#         α_rad = deg2rad(α)
#         ca = CA(S, α_rad, β, σ, T, Tw, lx, ly, lz)
#         cn = CN(S, α_rad, β, σ, T, Tw, lx, ly, lz)
#         cs = CS(S, α_rad, β, σ, T, Tw, lx, ly, lz)
#         ca_old = axialcoefficient(S, α_rad, σ)
#         cn_old = normalcoefficient(S, α_rad, σ)
#         push!(ca_list, ca)
#         push!(cn_list, cn)
#         push!(cs_list, cs)
#         push!(ca_old_list, ca_old)
#         push!(cn_old_list, cn_old)
#     end
#     plot!(p, [min_angle:max_angle], ca_old_list, linestyle=:dash, label=" CA (old)", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
#     plot!(p, [min_angle:max_angle], cn_old_list, linestyle=:dash, label=" CN (old)", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
#     plot!(p, [min_angle:max_angle], ca_list, label=" CA", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
#     plot!(p, [min_angle:max_angle], cn_list, label=" CN", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
#     # plot!(p, [min_angle:max_angle], cs_list, label=" CS", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
# end
# display(p)

# p = plot()
# for σ in σ_list
#     cd_list = Float64[]
#     cl_list = Float64[]
#     cs_list = Float64[]
#     cd_old_list = Float64[]
#     cl_old_list = Float64[]
#     for α in min_angle:max_angle
#         α_rad = deg2rad(α)
#         ca = CA(S, α_rad, β, σ, T, Tw, lx, ly, lz)
#         cn = CN(S, α_rad, β, σ, T, Tw, lx, ly, lz)
#         cs = CS(S, α_rad, β, σ, T, Tw, lx, ly, lz)
#         ca_old = axialcoefficient(S, α_rad, σ)
#         cn_old = normalcoefficient(S, α_rad, σ)
#         cd = cos(α_rad)*cos(β)*ca + sin(β)*cs + sin(α_rad)*cos(β)*cn
#         cl = -sin(α_rad)*ca + cos(α_rad)*cn
#         cs = cos(α_rad)*sin(β)*ca - cos(β)*cs - sin(α_rad)*sin(β)*cn
#         cd_old = cos(α_rad)*ca_old + sin(α_rad)*cn_old
#         cl_old = -sin(α_rad)*ca_old + cos(α_rad)*cn_old
#         push!(cd_list, cd)
#         push!(cl_list, cl)
#         push!(cs_list, cs)
#         push!(cd_old_list, cd_old)
#         push!(cl_old_list, cl_old)
#     end
#     plot!(p, [min_angle:max_angle], cd_list, label=" CD", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
#     plot!(p, [min_angle:max_angle], cl_list, label=" CL", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
#     # plot!(p, [min_angle:max_angle], cs_list, label=" CS", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
#     plot!(p, [min_angle:max_angle], cd_old_list, linestyle=:dash, label=" CD (old)", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
#     plot!(p, [min_angle:max_angle], cl_old_list, linestyle=:dash, label=" CL (old)", xlabel="Angle of attack (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Attack")
# end
# display(p)

# p = plot()
# for σ in σ_list
#     cd_list = Float64[]
#     cl_list = Float64[]
#     cs_list = Float64[]
#     for β in min_angle:max_angle
#         β_rad = deg2rad(β)
#         ca = CA(S, α, β_rad, σ, T, Tw, lx, ly, lz)
#         cn = CN(S, α, β_rad, σ, T, Tw, lx, ly, lz)
#         cs = CS(S, α, β_rad, σ, T, Tw, lx, ly, lz)
#         cd = cos(α)*cos(β_rad)*ca + sin(β_rad)*cs - sin(α)*cos(β_rad)*cn
#         cl = sin(α)*ca + cos(α)*cn
#         cs = cos(α)*sin(β_rad)*ca + cos(β_rad)*cs - sin(α)*sin(β_rad)*cn
#         push!(cd_list, cd)
#         push!(cl_list, cl)
#         push!(cs_list, cs)
#     end
#     plot!(p, [min_angle:max_angle], cd_list, label=" CD", xlabel="Angle of sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Sideslip")
#     plot!(p, [min_angle:max_angle], cl_list, label=" CL", xlabel="Angle of sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Sideslip")
#     plot!(p, [min_angle:max_angle], cs_list, label=" CS", xlabel="Angle of sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="Aerodynamic Coefficient vs Angle of Sideslip")
# end
# display(p)

# using Arrow
# using DataFrames
# using UUIDs # For generating unique filenames

# # Create a temporary directory for the dataset.
# # The `try...finally` block ensures the directory is always cleaned up.
# temp_dir_path = mktempdir()
# println("✅ Created temporary dataset directory at: $(temp_dir_path)")

# try
#     ## 1. WRITE THE FIRST FILE
#     println("   Writing part 1...")
#     df1 = DataFrame(id = 1:3, category = ["A", "B", "A"], value = rand(3))
#     Arrow.write(joinpath(temp_dir_path, "part_1"), df1)

#     ## 2. WRITE A SECOND FILE
#     println("   Writing part 2...")
#     df2 = DataFrame(id = 4:5, category = ["C", "B"], value = rand(2))
#     # Using a UUID ensures the filename is always unique.
#     unique_name = "data-1"
#     Arrow.write(joinpath(temp_dir_path, unique_name), df2)

#     println("Files in temporary directory:")
#     # List all files in the temporary directory to confirm they were written.
#     for file in readdir(temp_dir_path)
#         println(" - $(file)")
#     end
#     ## 3. READ THE ENTIRE DIRECTORY AS A SINGLE TABLE
#     println("\n✅ Reading the entire dataset back...")
#     # Arrow.Table automatically finds and combines all .arrow files in the directory.
#     full_table = Arrow.Table(temp_dir_path)
#     full_df = DataFrame(full_table)

#     println("Combined DataFrame:")
#     println(full_df)

# finally
#     # This `finally` block guarantees that the temporary directory and all its
#     # contents are removed, even if an error happened in the `try` block.
#     println("\n✅ Cleaning up temporary directory...")
#     rm(temp_dir_path; recursive=true)
# end

using SpecialFunctions
using LinearAlgebra
const sqrt_π = sqrt(π)
const sqrt_2 = sqrt(2)
"""
    calculate_CA(t1, t2, t3, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf, θ)

Calculates the expression C_A from the provided image.
This version incorporates the clarification that θ is a scalar variable (an angle).

Arguments:
- t1, t2, t3: t₁, t₂, t₃
- cos_α, sin_α, cos_β, sin_β: trig functions of angle of attack and sideslip angles (angles)
- θ: theta (flow angle)
- σ_T, σ_N: sigma_T, sigma_N
- s: s
- Tw: T_w (Wall Temperature)
- T_inf: T_∞ (Temperature at infinity)
"""
function calculate_CA(t1, t2, t3, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf, θ)
    # --- Pre-calculate common values ---
    s_sq = s^2
    sqrt_Tw_Tinf = sqrt(Tw / T_inf)

    # --- Group 1: Terms 1 and 2 ---
    s_sin_β = s * sin_β
    erf_s_sin_β = erf(s_sin_β)
    exp_s_sin_β_sq = exp(-s_sin_β^2)

    paren_term_minus_1 = sqrt_π * s_sin_β * (erf_s_sin_β - 1) + exp_s_sin_β_sq
    paren_term_plus_1 = sqrt_π * s_sin_β * (erf_s_sin_β + 1) + exp_s_sin_β_sq
    
    # Denominator simplifies because sqrt(sec(β)²) * abs(cos(β)) = 1
    den1_2 = sqrt_π * s * t2
    common_factor_1_2 = t1 * cos_α * cos_β * σ_T * θ
    
    term1 = (common_factor_1_2 * (-sin_β) * paren_term_minus_1) / den1_2
    term2_frac = (common_factor_1_2 * sin_β * paren_term_plus_1) / den1_2

    # --- Group 2: Terms 3 and 4 (no denominators) ---
    s_ca_cb = s * cos_α * cos_β
    erf_s_cacb = erf(s_ca_cb)
    exp_s_cacb_sq = exp(-s_ca_cb^2)
    s_sq_ca_sq_cb_sq = s_sq * cos_α^2 * cos_β^2
    
    # Term 3
    exp_bracket3 = (s_ca_cb * (σ_N - 2)) / sqrt_π + 0.5 * σ_N * sqrt_Tw_Tinf
    main_bracket3 = (1 - erf_s_cacb) * ((2 - σ_N) * (s_sq_ca_sq_cb_sq + 0.5) - 
                    0.5 * sqrt_π * s_ca_cb * σ_N * sqrt_Tw_Tinf) + 
                    exp_s_cacb_sq * exp_bracket3
    term3 = 1/s_sq * θ * (-cos_α) * cos_β * main_bracket3
    
    # Term 4
    exp_bracket4 = 0.5 * σ_N * sqrt_Tw_Tinf - (s_ca_cb * (σ_N - 2)) / sqrt_π
    main_bracket4 = (erf_s_cacb + 1) * ((2 - σ_N) * (s_sq_ca_sq_cb_sq + 0.5) +
                    0.5 * sqrt_π * s_ca_cb * σ_N * sqrt_Tw_Tinf) +
                    exp_s_cacb_sq * exp_bracket4
    term4 = 1/s_sq * θ * cos_α * cos_β * main_bracket4

    # --- Group 3: Terms 5 and 6 ---
    s_sa_cb = s * sin_α * cos_β
    erf_s_sacb = erf(s_sa_cb)
    exp_s_sacb_sq = exp(-s_sa_cb^2)

    paren_term5 = sqrt_π * s_sa_cb * (erf_s_sacb - 1) + exp_s_sacb_sq
    paren_term6 = sqrt_π * s_sa_cb * (erf_s_sacb + 1) + exp_s_sacb_sq

    den5_6 = sqrt_π * s * t3
    common_factor_5_6 = t1 * cos_α * cos_β * σ_T * θ
    
    term5 = (common_factor_5_6 * (-cos_β * sin_α) * paren_term5) / den5_6
    term6 = (common_factor_5_6 * (cos_β * sin_α) * paren_term6) / den5_6

    # --- Final Combination ---
    # Based on the signs shown in the image:
    # T1 + T2_frac - 1/s² - T3 + T4 + T5 + T6
    Ca_result = term1 + term2_frac - term3 + term4 + term5 + term6
    
    return Ca_result
end

"""
    calculate_CS(t1, t2, t3, θ, α, β, σ_T, σ_N, s, Tw, T_inf)

Calculates the expression C_S from the provided image.

Arguments:
- t1, t2, t3: t₁, t₂, t₃
- θ, α, β: theta, alpha, beta (angles)
- σ_T, σ_N: sigma_T, sigma_N
- s: s
- Tw: T_w (Wall Temperature)
- T_inf: T_∞ (Temperature at infinity)
"""
function calculate_CS(t1, t2, t3, θ, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf)
    # --- Pre-calculate common values ---
    s_sq = s^2
    sqrt_Tw_Tinf = sqrt(Tw / T_inf)

    # --- Group 1: Terms 1 and 2 ---
    s_sin_β = s * sin_β
    erf_s_sin_β = erf(s_sin_β)
    exp_s_sin_β_sq = exp(-s_sin_β^2)
    s_sq_sin_β_sq = s_sq * sin_β^2

    # Term 1
    exp_bracket1 = (s_sin_β * (σ_N - 2)) / sqrt_π + 0.5 * σ_N * sqrt_Tw_Tinf
    main_bracket1 = (1 - erf_s_sin_β) * ((2 - σ_N) * s_sq_sin_β_sq + 0.5) -
                    0.5 * sqrt_π * s_sin_β * σ_N * sqrt_Tw_Tinf +
                    exp_s_sin_β_sq * exp_bracket1
    num1 = t1 * θ * (-sin_β) * main_bracket1
    den1_2 = s_sq * t2
    term1 = num1 / den1_2

    # Term 2
    exp_bracket2 = 0.5 * σ_N * sqrt_Tw_Tinf - (s_sin_β * (σ_N - 2)) / sqrt_π
    main_bracket2 = (erf_s_sin_β + 1) * ((2 - σ_N) * s_sq_sin_β_sq + 0.5) +
                    0.5 * sqrt_π * s_sin_β * σ_N * sqrt_Tw_Tinf +
                    exp_s_sin_β_sq * exp_bracket2
    num2 = t1 * θ * sin_β * main_bracket2
    term2 = num2 / den1_2

    # --- Group 2: Terms 3 and 4 ---
    s_ca_cb = s * cos_α * cos_β
    erf_s_cacb = erf(s_ca_cb)
    exp_s_cacb_sq = exp(-s_ca_cb^2)
    
    # Term 3
    paren_term3 = sqrt_π * s_ca_cb * (erf_s_cacb - 1) + exp_s_cacb_sq
    num3 = sin_β * σ_T * θ * (-cos_α * cos_β) * paren_term3
    den3_4 = sqrt_π * s
    term3 = num3 / den3_4

    # Term 4
    paren_term4 = sqrt_π * s_ca_cb * (erf_s_cacb + 1) + exp_s_cacb_sq
    num4 = sin_β * σ_T * θ * (cos_α * cos_β) * paren_term4
    term4 = num4 / den3_4

    # --- Group 3: Terms 5 and 6 ---
    s_sa_cb = s * sin_α * cos_β
    erf_s_sacb = erf(s_sa_cb)
    exp_s_sacb_sq = exp(-s_sa_cb^2)

    paren_term5 = sqrt_π * s_sa_cb * (erf_s_sacb - 1) + exp_s_sacb_sq
    paren_term6 = sqrt_π * s_sa_cb * (erf_s_sacb + 1) + exp_s_sacb_sq

    den5_6 = sqrt_π * s * t3
    common_factor_5_6_num = t1 * sin_β * σ_T * θ
    
    num5 = common_factor_5_6_num * (-cos_β * sin_α) * paren_term5
    term5 = num5 / den5_6

    num6 = common_factor_5_6_num * (cos_β * sin_α) * paren_term6
    term6 = num6 / den5_6
    
    # --- Final Combination ---
    # Based on the signs shown at the start and end of each line in the image
    Cs_result = -term1 + term2 + term3 + term4 + term5 + term6
    
    return Cs_result
end

"""
    calculate_CN(t1, t2, t3, θ, α, β, σ_T, σ_N, s, Tw, T_inf)

Calculates the expression C_N from the provided image.

Arguments:
- t1, t2, t3: t₁, t₂, t₃
- θ, α, β: theta, alpha, beta (angles)
- σ_T, σ_N: sigma_T, sigma_N
- s: s
- Tw: T_w (Wall Temperature)
- T_inf: T_∞ (Temperature at infinity)
"""
function calculate_CN(t1, t2, t3, θ, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf)
    # --- Pre-calculate common values ---
    s_sq = s^2
    sqrt_Tw_Tinf = sqrt(Tw / T_inf)

    # --- Group 1: Terms 1 and 2 ---
    s_sin_β = s * sin_β
    erf_s_sin_β = erf(s_sin_β)
    exp_s_sin_β_sq = exp(-s_sin_β^2)

    paren_term_minus_1 = sqrt_π * s_sin_β * (erf_s_sin_β - 1) + exp_s_sin_β_sq
    paren_term_plus_1 = sqrt_π * s_sin_β * (erf_s_sin_β + 1) + exp_s_sin_β_sq
    
    # Denominator simplifies because sqrt(sec(β)²) * abs(cos(β)) = 1
    den1_2 = sqrt_π * s * t2
    common_factor_1_2_num = t1 * sin_α * cos_β * σ_T * θ
    
    num1 = common_factor_1_2_num * (-sin_β) * paren_term_minus_1
    term1 = num1 / den1_2

    num2_frac = common_factor_1_2_num * sin_β * paren_term_plus_1
    term2_frac = num2_frac / den1_2
    term2_sub = 1 / (s_sq * t3)

    # --- Group 2: Terms 3 and 4 ---
    s_sa_cb = s * sin_α * cos_β
    erf_s_sacb = erf(s_sa_cb)
    exp_s_sacb_sq = exp(-s_sa_cb^2)
    s_sq_sa_sq_cb_sq = s_sq * sin_α^2 * cos_β^2

    # Term 3
    exp_bracket3 = (s_sa_cb * (σ_N - 2)) / sqrt_π + 0.5 * σ_N * sqrt_Tw_Tinf + 1 / (s_sq * t3)
    main_bracket3 = (1 - erf_s_sacb) * ((2 - σ_N) * s_sq_sa_sq_cb_sq - 0.5) -
                    0.5 * sqrt_π * s_sa_cb * σ_N * sqrt_Tw_Tinf +
                    exp_s_sacb_sq * exp_bracket3
    term3 = t1 * θ * (-cos_β * sin_α) * main_bracket3
    
    # Term 4
    exp_bracket4 = 0.5 * σ_N * sqrt_Tw_Tinf - (s_sa_cb * (σ_N - 2)) / sqrt_π
    main_bracket4 = (erf_s_sacb + 1) * ((2 - σ_N) * s_sq_sa_sq_cb_sq + 0.5) +
                    0.5 * sqrt_π * s_sa_cb * σ_N * sqrt_Tw_Tinf +
                    exp_s_sacb_sq * exp_bracket4
    term4 = t1 * θ * (cos_β * sin_α) * main_bracket4

    # --- Group 3: Terms 5 and 6 ---
    s_ca_cb = s * cos_α * cos_β
    erf_s_cacb = erf(s_ca_cb)
    exp_s_cacb_sq = exp(-s_ca_cb^2)

    paren_term5 = sqrt_π * s_ca_cb * (erf_s_cacb - 1) + exp_s_cacb_sq
    paren_term6 = sqrt_π * s_ca_cb * (erf_s_cacb + 1) + exp_s_cacb_sq
    
    den5_6 = sqrt_π * s
    common_factor_5_6_num = sin_α * cos_β * σ_T * θ

    num5 = common_factor_5_6_num * (-cos_α * cos_β) * paren_term5
    term5 = num5 / den5_6
    
    num6 = common_factor_5_6_num * (cos_α * cos_β) * paren_term6
    term6 = num6 / den5_6
    
    # --- Final Combination ---
    # Based on the signs shown in the image: 
    # T1 + T2_frac - 1/(s²t₃) - T3 + T4 + T5 + T6
    Cn_result = term1 + term2_frac - term2_sub - term3 + term4 + term5 + term6
    
    return Cn_result
end

"""
    calculate_Cl(σ_T, s, α, β, θ)

Calculates the expression C_L from the provided image.

Arguments:
- σ_T: sigma_T
- s: s
- α, β, θ: alpha, beta, theta (angles)
"""
function calculate_Cl(σ_T, s, sin_α, cos_α, sin_β, cos_β, θ)
    # --- Pre-calculate common values ---

    # --- Calculation for the "Top Part" (first three lines) ---
    # This part is structured as: ( sin(β) * [ ... ] ) - 1/(sqrt(π)*|cos(β)|)
    
    s_sa_cb = s * sin_α * cos_β
    erf_s_sacb = erf(s_sa_cb)
    exp_s_sacb_sq = exp(-s_sa_cb^2)

    # The two θ terms can be combined and simplified algebraically:
    # θ(-x)A + θ(x)B = -θx A + θx B = θx(B - A)
    paren_A = s_sa_cb * erf_s_sacb + exp_s_sacb_sq / sqrt_π - s_sa_cb
    paren_B = -s_sa_cb * erf_s_sacb - exp_s_sacb_sq / sqrt_π - s_sa_cb
    
    term_in_paren_top = paren_B - paren_A # This simplifies to -2*s_sa_cb*erf - 2*exp/sqrt_π
    
    # Note: θ(-cos(β)sin(α)) becomes -θ*cos(β)sin(α) etc.
    sub_part_1 = sin_β * (θ * cos_β * sin_α * term_in_paren_top)
    
    sub_part_2 = 1 / (sqrt_π * abs(cos_β))
    
    top_part = sub_part_1 #- sub_part_2

    # --- Calculation for the "Bottom Part" (last three lines) ---
    # This is structured as: [prefactors] * ( θ(-sinβ)[...] - θ(sinβ)[...] )
    s_sin_β = s * sin_β
    exp_s_sin_β_sq_neg = exp(-s_sin_β^2)
    
    # Note: sqrt(sec(β)²) simplifies to 1/abs(cos(β))
    prefactor_bottom = sin_α * cos_β^3 * (1 / abs(cos_β)) * exp_s_sin_β_sq_neg

    # The bracketed terms can also be simplified:
    exp_s_sin_β_sq_pos = exp(s_sin_β^2)
    erf_s_sin_β = erf(s_sin_β)
    common_in_brackets = sqrt_π * s_sin_β * exp_s_sin_β_sq_pos * erf_s_sin_β
    
    bracket_A = common_in_brackets + sqrt_π * s_sin_β * (-exp_s_sin_β_sq_neg + 1)
    bracket_B = common_in_brackets + sqrt_π * s_sin_β * (exp_s_sin_β_sq_neg + 1)
    
    # θ(-sinβ)A - θ(sinβ)B = -θ*sinβ*A - θ*sinβ*B = -θ*sinβ*(A + B)
    sum_of_brackets = bracket_A + bracket_B # Simplifies to 2*common + 2*sqrt(π)*s_sin_β

    main_paren_bottom = -θ * sin_β * sum_of_brackets

    bottom_part = prefactor_bottom * main_paren_bottom
    
    # --- Final Combination ---
    total_inside_braces = top_part - sub_part_2 * bottom_part
    
    Cl_result = (1 / (2 * s)) * σ_T * total_inside_braces
    
    return Cl_result
end

"""
    calculate_Cm(σ_T, s, α, β, θ)

Calculates the expression C_m from the provided image.

Arguments:
- σ_T: sigma_T
- s: s
- α, β, θ: alpha, beta, theta (angles)
"""
function calculate_Cm(σ_T, s, sin_α, cos_α, sin_β, cos_β, θ)
    # --- Pre-calculate common values ---

    # --- Define arguments for the two main patterns ---
    arg1 = s * sin_α * cos_β
    arg2 = s * cos_α * cos_β

    # --- Term 1 ---
    erf_arg1 = erf(arg1)
    exp_arg1_sq = exp(-arg1^2)
    
    # NOTE: The structure of this parenthesis is unusual compared to the others.
    # It is translated literally as it appears in the image.
    paren1 = sin_α * cos_β - sin_α * cos_β * erf_arg1 - exp_arg1_sq / (sqrt_π * s)
    # term1 = cos_α * (θ * (-cos_β * sin_α)) * paren1

    # --- Term 2 ---
    paren2 = arg1 * erf_arg1 + exp_arg1_sq / sqrt_π + arg1
    # term2 = (1/s) * cos_α * (θ * (cos_β * sin_α)) * paren2

    # --- Term 3 ---
    erf_arg2 = erf(arg2)
    exp_arg2_sq = exp(-arg2^2)
    
    paren3 = arg2 * erf_arg2 + exp_arg2_sq / sqrt_π - arg2
    # term3 = sin_α * (θ * (-cos_α * cos_β)) * paren3

    # --- Term 4 ---
    # NOTE: This term is not prefixed by sin(α), unlike its counterpart, Term 3.
    # This asymmetry is preserved from the image.
    paren4 = -arg2 * erf_arg2 - exp_arg2_sq / sqrt_π - arg2
    term4 = (θ * (cos_α * cos_β)) * paren4

    # --- Final Combination ---
    total_sum = cos_α*θ*(-cos_β*sin_α)*paren1 + 1/s*(cos_α*θ*cos_β*cos_α*paren2 + sin_α*(θ * (-cos_α * cos_β) * paren3 + term4))
    
    Cm_result = 0.5 * cos_β * σ_T * total_sum
    
    return Cm_result
end

"""
    calculate_Cn(σ_T, s, α, β, θ)

Calculates the expression C_n from the provided image.

Arguments:
- σ_T: sigma_T
- s: s
- α, β, θ: alpha, beta, theta (angles)
"""
function calculate_Cn(σ_T, s, sin_α, cos_α, sin_β, cos_β, θ)
    # --- Pre-calculate common values ---
    s_sq = s^2
    abs_cos_β = abs(cos_β)

    # --- Define common exponential terms ---
    exp_arg_cos_term = s_sq * cos_α^2 * cos_β^2
    exp_arg_sin_term = s_sq * sin_β^2
    
    exp_cos_pos = exp(exp_arg_cos_term)
    exp_cos_neg = exp(-exp_arg_cos_term)
    exp_sin_pos = exp(exp_arg_sin_term)
    
    exp_main_pos = exp(exp_arg_cos_term + exp_arg_sin_term)
    exp_main_neg = exp(-(exp_arg_cos_term + exp_arg_sin_term))

    # --- Outer Factor ---
    outer_factor = (1 / (2 * sqrt_π * s * abs_cos_β)) * σ_T * exp_main_neg

    # --- Part 1: Terms factored by (-sin(β)|cos(β)|) ---
    s_ca_cb = s * cos_α * cos_β
    
    # The two θ terms inside the first main parenthesis can be combined and simplified.
    # The structure is prop. to θ*cosαcosβ * (Paren_B - Paren_A)
    # (Paren_B - Paren_A) simplifies to: 2 * sqrt(π) * s_ca_cb * exp(-s²cos²αcos²β)
    diff_paren_AB = 2 * sqrt_π * s_ca_cb * exp_cos_neg
    combo_AB = θ * cos_α * cos_β * diff_paren_AB
    
    part1 = (-sin_β * abs_cos_β) * combo_AB
    
    # --- Part 2: Terms factored by (cos³(β)sqrt(sec²(β))) ---
    s_sin_β = s * sin_β
    erf_s_sinb = erf(s_sin_β)

    # The two θ terms inside the second main parenthesis can also be simplified.
    # The structure is prop. to -θ*sinβ * (Paren_C + Paren_D)
    # (Paren_C + Paren_D) simplifies to: 
    # 2 * (sqrt(π)*s*sinβ*erf(s*sinβ)*exp_main_pos + exp(s²cos²αcos²β))
    sum_paren_CD = 2 * (sqrt_π * s_sin_β * erf_s_sinb * exp_main_pos + exp_cos_pos)
    combo_CD = θ * (-sin_β) * sum_paren_CD
    
    # Prefactor for Part 2 simplifies: cos³(β)sqrt(sec²(β)) = cos³(β)/|cos(β)|
    prefactor2 = cos_β^3 / abs_cos_β
    part2 = prefactor2 * combo_CD

    # --- Final Combination ---
    total_sum_in_brackets = part1 + part2
    
    Cn_result = outer_factor * total_sum_in_brackets
    
    return Cn_result
end

# Test cases
using Plots 
plotly()
# Parameters for testing
Tw = 300.0  # Wall temperature in K
T_inf = 1100.0  # Temperature at infinity in K
lx, ly, lz = 0.1, 0.1, 0.1  # Characteristic lengths in m
S = lx * ly  # Reference area in m²
V_inf = 7500.0  # Free-stream velocity in m/s
α = deg2rad(10.0)  # Angle of attack in radians
β = deg2rad.(0.0:1.0:90.0)   # Angle of sideslip in radians
θ = deg2rad(90.0)  # Flow angle in radians
σ_T = 0.0  # Tangential accommodation coefficient
σ_N = 0.0  # Normal accommodation coefficient

CA_list = zeros(length(β))
CN_list = zeros(length(β))
CS_list = zeros(length(β))
Cl_list = zeros(length(β))
Cm_list = zeros(length(β))
Cn_list = zeros(length(β))
R = 287.0  # Specific gas constant for air in J/(kg·K)
s = V_inf / sqrt(2*R*T_inf)  # Speed ratio
println("CA test: ", calculate_CA(lx, ly, lz, sin(α), cos(α), sin(β[1]), cos(β[1]), σ_T, σ_N, s, Tw, T_inf, θ))
for (i, β_i) in enumerate(β)
    sin_α = sin(α)
    cos_α = cos(α)
    sin_β = sin(β_i)
    cos_β = cos(β_i)

    CA_list[i] = calculate_CA(lx, ly, lz, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf, θ)
    CN_list[i] = calculate_CN(lx, ly, lz, θ, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf)
    CS_list[i] = calculate_CS(lx, ly, lz, θ, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf)
    Cl_list[i] = calculate_Cl(σ_T, s, sin_α, cos_α, sin_β, cos_β, θ)
    Cm_list[i] = calculate_Cm(σ_T, s, sin_α, cos_α, sin_β, cos_β, θ)
    Cn_list[i] = calculate_Cn(σ_T, s, sin_α, cos_α, sin_β, cos_β, θ)
end


display(plot(rad2deg.(β), CA_list, label="C_A", xlabel="Angle of Sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="C_A vs Angle of Sideslip"))
display(plot(rad2deg.(β), CN_list, label="C_N", xlabel="Angle of Sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="C_N vs Angle of Sideslip"))
display(plot(rad2deg.(β), CS_list, label="C_S", xlabel="Angle of Sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="C_S vs Angle of Sideslip"))
display(plot(rad2deg.(β), Cl_list, label="C_l", xlabel="Angle of Sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="C_L vs Angle of Sideslip"))
display(plot(rad2deg.(β), Cm_list, label="C_m", xlabel="Angle of Sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="C_m vs Angle of Sideslip"))
display(plot(rad2deg.(β), Cn_list, label="C_n", xlabel="Angle of Sideslip (degrees)", ylabel="Aerodynamic Coefficient", title="C_n vs Angle of Sideslip"))