include("../physical_models/MonteCarlo_pertrubations.jl")

using SpecialFunctions

import .config

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

    # Solar Panels
    CN_sa = normalcoefficient(S, α, σ)
    CA_sa = axialcoefficient(S, α, σ)
    CL_sa = CN_sa*cos(α) - CA_sa*sin(α)
    CD_sa = CA_sa*cos(α) + CN_sa*sin(α)

    # Spacecraft
    CN_sc = normalcoefficient(S, pi*0.5, σ)
    CA_sc = axialcoefficient(S, pi*0.5, σ)
    CL_sc = CN_sc*cos(pi*0.5) - CA_sc*sin(pi*0.5)
    CD_sc = CA_sc*cos(pi*0.5) + CN_sc*sin(pi*0.5)

    CD_body = (CD_sa*body.area_SA + CD_sc*body.area_SC) / (body.area_SA + body.area_SC)
    CL_body = (CL_sa*body.area_SA + CL_sc*body.area_SC) / (body.area_SA + body.area_SC)

    if montecarlo == true
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)
    end

    return CL_body, CD_body
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