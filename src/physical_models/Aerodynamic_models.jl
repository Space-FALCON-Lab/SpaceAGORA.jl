include("../physical_models/MonteCarlo_pertrubations.jl")

using SpecialFunctions

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

function aerodynamic_coefficient_fM(α::Float64, β::Float64, body, T::Float64, S::Float64, args, montecarlo::Bool=false)
    """
    Calculate the aerodynamic coefficients for a blunt body in ballistic flight using the F.M. model.
        Angle of attack and side slip are calculated as angles between CA (normal to flat plate), and wind-relative velocity vector.
    # Arguments
    - `α`: Angle of attack, rad
    - `β`: Angle of sideslip, rad
    - `body`: Body object containing dimensions and properties
    - `T`: Temperature, K
    - `S`: Surface area, m²
    - `args`: Dictionary containing additional parameters like reflection coefficient
    - `montecarlo`: Boolean flag to apply Monte Carlo perturbations
    # Returns
    - `CL`: Lift coefficient
    - `CD`: Drag coefficient
    - `CS`: Sideslip coefficient
    # Notes
    - Equations are from Hart et al. (2017, https://doi.org/10.2514/1.A33606), for a rectangular prism. 
    """

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

    # Calculate the aerodynamic coefficients in the body frame
    # Axial
    CA = ((2-σN)/(S*√(π))*cosα*cosβ+sign(cosα*cosβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*cosα^2*cosβ^2) +
            (2-σN)*(cosα^2*cosβ^2+1/(2*S^2)) * (sign(cosα*cosβ)+erf(S*cosα*cosβ)) + 
            (σN/(2*S)*cosα*cosβ*√(π*Tw/T)) * (1+sign(cosα*cosβ)*erf(S*cosα*cosβ)) +
            σT*cosα*cosβ*(lx/ly*(1/(S*√(π))*exp(-S^2*sinβ^2)+sinβ*(sign(sinβ)+erf(S*sinβ))) +
            lx/lz*(1/(S*√(π))*exp(-S^2*sinα^2*cosβ^2)+sinα*cosβ*(sign(sinα*cosβ)+erf(S*sinα*cosβ))))
    # Crossflow
    CS = lx/ly*(((2-σN)/(S*√(π))*sinβ+sign(sinβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*sinβ^2) +
                (2-σN)*(sinβ^2+1/(2*S^2)) * (sign(sinβ)+erf(S*sinβ)) + 
                (σN/(2*S)*sinβ*√(π*Tw/T)) * (1+sign(sinβ)*erf(S*sinβ))) +
                σT*sinβ*(1/(S*√(π))*exp(-S^2*cosα^2*cosβ^2) + cosα*cosβ*(sign(cosα*cosβ)+erf(S*cosα*cosβ)) + 
                lx/lz*(1/(S*√(π))*exp(-S^2*sinα^2*cosβ^2) + sinα*cosβ*(sign(sinα*cosβ)+erf(S*sinα*cosβ))))
    # Normal
    CN = lx/lz*((((2-σN)/(S*√(π))*sinα*cosβ+sign(sinα*cosβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*sinα^2*cosβ^2) +
                (2-σN)*(sinα^2*cosβ^2+1/(2*S^2)) * (sign(sinα*cosβ)+erf(S*sinα*cosβ)) + 
                (σN/(2*S)*sinα*cosβ*√(π*Tw/T)) * (1+sign(sinα*cosβ)*erf(S*sinα*cosβ)))) +
                σT*sinα*cosβ*(lx/ly*(1/(S*√(π))*exp(-S^2*sinβ^2) + sinβ*(sign(sinβ)+erf(S*sinβ))) + 
                (1/(S*√(π))*exp(-S^2*cosα^2*cosβ^2) + cosα*cosβ*(sign(cosα*cosβ)+erf(S*cosα*cosβ))))

    # Calculate the aerodynamic coefficients in the wind frame
    CL = -sinα*CA + cosα*CN
    CD = cosα*cosβ*CA + sinβ*CS + sinα*cosβ*CN
    CS = cosα*sinβ*CA + cosβ*CS - sinα*sinβ*CN

    # If doing Monte Carlo simulations, apply perturbations to the coefficients
    if montecarlo == true
        CL, CD = monte_carlo_aerodynamics(CL, CD, args)
    end

    return CL, CD, CS
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