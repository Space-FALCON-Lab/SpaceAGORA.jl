include("../physical_models/MonteCarlo_pertrubations.jl")

using SpecialFunctions

const sqrt_π = sqrt(π)
const inv_sqrt_π = 1 / sqrt(π)

function aerodynamic_coefficient_constant(α, body, T, S, args, montecarlo=false)
    """

    """

    CL_body = 0.0
    CD_body = 2 * (2.2 - 0.8)/pi * args.α + 0.8

    if montecarlo == true
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)
    end

    return CL_body, CD_body
end

function aerodynamic_coefficient_fM(α, body, T, S, args, montecarlo=false)
    """

    """

    σ = args.reflection_coefficient
    Tw = T

    function pressure(S, α, ρ_inf, vel, σ)
        """

        """

        p = (ρ_inf*vel^2) / (2*S^2) * ((((2 - σ) / sqrt(pi))*S*sin(α) + sqrt(T/Tw)*σ/2) * exp(-(S*sin(α))^2) + ((2-σ)*((S*sin(α))^2 + 0.5) + (σ/2)*sqrt(pi)*(S*sin(α))) * (1 + erf(S*sin(α))))

        return p
    end

    function τ(S, α, ρ_inf, vel, σ)
        """

        """

        t = ((σ*cos(α)*ρ_inf*vel^2) / (sqrt(pi)*2*S)) * (exp(-(S*sin(α))^2) + sqrt(pi)*(S*sin(α)) * (1 + erf(S*sin(α))))

        return t
    end

    function normalcoefficient(S, aoa, sigma)
        CN = 1 / (S^2) * ((((2 - sigma) / sqrt(pi)) * S * sin(aoa) + sigma / 2) * exp(-(S * sin(aoa))^2) + ((2 - sigma) * ((S * sin(aoa))^2 + 0.5) + sigma / 2 * sqrt(pi) * (S * sin(aoa))) * (1 + erf(S * sin(aoa))))
        return CN
    end

    function axialcoefficient(S, aoa, sigma)
        CA = ((sigma * cos(aoa)) / (sqrt(pi) * S)) * (exp(-(S * sin(aoa))^2) + sqrt(pi) * (S * sin(aoa)) * (1 + erf(S * sin(aoa))))
        return CA
    end

    # println("α: ", α)

    # Solar Panels
    CN_sa = normalcoefficient(S, α, σ)
    CA_sa = axialcoefficient(S, α, σ)
    CL_sa = CN_sa*cos(α) - CA_sa*sin(α)
    CD_sa = CA_sa*cos(α) + CN_sa*sin(α)
    # println("CL_sa: ", CL_sa)
    # println("CD_sa: ", CD_sa)
    # Spacecraft
    # CN_sc = normalcoefficient(S, pi*0.5, σ)
    # CA_sc = axialcoefficient(S, pi*0.5, σ)
    # CL_sc = CN_sc*cos(pi*0.5) - CA_sc*sin(pi*0.5)
    # CD_sc = CA_sc*cos(pi*0.5) + CN_sc*sin(pi*0.5)

    # area_SA = config.get_SA_area(body, body.roots[1])
    # area_SC = config.get_SC_area(body, body.roots[1])
    
    # CD_body = (CD_sa*area_SA + CD_sc*area_SC) / (area_SA + area_SC)
    # CL_body = (CL_sa*area_SA + CL_sc*area_SC) / (area_SA + area_SC)
    if montecarlo == true
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)
    end

    return CL_sa, CD_sa
end

function aerodynamic_coefficient_fM(body, T::Float64, S::Float64, args, montecarlo::Bool=false)
    """
    Calculate the aerodynamic coefficients for a blunt body in ballistic flight using the F.M. model.
        Angle of attack and side slip are calculated as angles between CA (normal to flat plate), and wind-relative velocity vector.
    # Arguments
    - `body`: Body object containing dimensions and properties
    - `T`: Temperature, K
    - `S`: Molecular speed ratio
    - `args`: Dictionary containing additional parameters like reflection coefficient
    - `montecarlo`: Boolean flag to apply Monte Carlo perturbations
    # Returns
    - `CL`: Lift coefficient
    - `CD`: Drag coefficient
    - `CS`: Sideslip coefficient
    # Notes
    - Equations are from Hart et al. (2017, https://doi.org/10.2514/1.A33606), for a rectangular prism. 
    """

    α = body.α
    β = body.β
    θ = body.θ
    # Adjust angles to match model assumptions
    α -= π/2
    σ = args.reflection_coefficient
    Tw = T
    lx = body.dims[1]
    ly = body.dims[2]
    lz = body.dims[3]
    σN = σ
    σT = σ
    cosα = cos(α)
    cosβ = cos(β)
    sinβ = sin(β)
    sinα = sin(α)
    # println("cosα: ", cosα)
    # println("cosβ: ", cosβ)
    # println("sinα: ", sinα)
    # println("sinβ: ", sinβ)
    # Calculate the aerodynamic coefficients in the body frame (flat plate)
    # Axial
    CA = ((2-σN)/(S*sqrt_π)*cosα*cosβ+sign(cosα*cosβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*cosα^2*cosβ^2) +
            (2-σN)*(cosα^2*cosβ^2+1/(2*S^2)) * (sign(cosα*cosβ)+erf(S*cosα*cosβ)) + 
            (σN/(2*S)*cosα*cosβ*√(π*Tw/T)) * (1+sign(cosα*cosβ)*erf(S*cosα*cosβ)) +
            σT*cosα*cosβ*(lx/ly*(1/(S*sqrt_π)*exp(-S^2*sinβ^2)+sinβ*(sign(sinβ)+erf(S*sinβ))) +
            lx/lz*(1/(S*sqrt_π)*exp(-S^2*sinα^2*cosβ^2)+sinα*cosβ*(sign(sinα*cosβ)+erf(S*sinα*cosβ))))
    # Crossflow
    CS = lx/ly*(((2-σN)/(S*sqrt_π)*sinβ+sign(sinβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*sinβ^2) +
                (2-σN)*(sinβ^2+1/(2*S^2)) * (sign(sinβ)+erf(S*sinβ)) + 
                (σN/(2*S)*sinβ*√(π*Tw/T)) * (1+sign(sinβ)*erf(S*sinβ))) +
                σT*sinβ*(1/(S*sqrt_π)*exp(-S^2*cosα^2*cosβ^2) + cosα*cosβ*(sign(cosα*cosβ)+erf(S*cosα*cosβ)) + 
                lx/lz*(1/(S*sqrt_π)*exp(-S^2*sinα^2*cosβ^2) + sinα*cosβ*(sign(sinα*cosβ)+erf(S*sinα*cosβ))))
    # Normal
    CN = lx/lz*((((2-σN)/(S*sqrt_π)*sinα*cosβ+sign(sinα*cosβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*sinα^2*cosβ^2) +
                (2-σN)*(sinα^2*cosβ^2+1/(2*S^2)) * (sign(sinα*cosβ)+erf(S*sinα*cosβ)) + 
                (σN/(2*S)*sinα*cosβ*√(π*Tw/T)) * (1+sign(sinα*cosβ)*erf(S*sinα*cosβ)))) +
                σT*sinα*cosβ*(lx/ly*(1/(S*sqrt_π)*exp(-S^2*sinβ^2) + sinβ*(sign(sinβ)+erf(S*sinβ))) + 
                (1/(S*sqrt_π)*exp(-S^2*cosα^2*cosβ^2) + cosα*cosβ*(sign(cosα*cosβ)+erf(S*cosα*cosβ))))

    # println("CA: ", CA)
    # println("CS: ", CS)
    # println("CN: ", CN)
    # Calculate the aerodynamic coefficients in the body frame (box)
    # Axial
    # CA = calculate_CA(lx, ly, lz, sinα, cosα, sinβ, cosβ, σT, σN, S, Tw, T, θ)
    # # Crossflow
    # CS = calculate_CS(lx, ly, lz, θ, sinα, cosα, sinβ, cosβ, σT, σN, S, Tw, T)
    # # Normal
    # CN = calculate_CN(lx, ly, lz, θ, sinα, cosα, sinβ, cosβ, σT, σN, S, Tw, T)

    # # Calculate moment coefficients
    # Cl = calculate_Cl(σT, S, sinα, cosα, sinβ, cosβ, θ)
    # Cm = calculate_Cm(σT, S, sinα, cosα, sinβ, cosβ, θ)
    # Cn = calculate_Cn(σT, S, sinα, cosα, sinβ, cosβ, θ)

    # Calculate the aerodynamic coefficients in the wind frame
    CL = -sinα*CA + cosα*CN
    CD = cosα*cosβ*CA + sinβ*CS + sinα*cosβ*CN
    CS = -cosα*sinβ*CA + cosβ*CS - sinα*sinβ*CN

    # If doing Monte Carlo simulations, apply perturbations to the coefficients
    if montecarlo == true
        CL, CD = monte_carlo_aerodynamics(CL, CD, args)
    end

    # return CL, CD, CS, 0.0, 0.0, 0.0#, Cl, Cm, Cn
    return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
end

function aerodynamic_coefficient_no_ballistic_flight(α, body, args, T=0, S=0, a=0, montecarlo=false)
    """

    """

    # Newtonian flow
    NoseRadius = body.nose_radius
    BaseRadius = body.base_radius

    k = NoseRadius / BaseRadius

    Cp_max = 2
    δ = body.δ

    CA_body = (1 - sin(δ)^4)*k^2 + (2*sin(δ)^2*cos(α)^2 + cos(δ)^2*sin(α)^2) * (1 - (k*cos(δ))^2)
    CN_body = (1 - (k*cos(δ))^2) * cos(δ^2) *sin(2*α)
    
    CD_body = CA_body*cos(α) + CN_body*sin(α) - 0.15
    CL_body = CN_body*cos(α) - CA_body*sin(α)

    if montecarlo == true
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)
    end

    return CL_body, CD_body
end

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