using Interpolations
using SpecialFunctions
using ..AbstractTypes: AbstractThermalModel, AbstractPlanet

@kwdef struct MaxwellianHeat{P <: AbstractPlanet} <: AbstractThermalModel
    thermal_accomodation_factor::Float64
    planet::P
    thermal_contact::Bool = false
end

# Blunt body
function heatrate_convective(S, T, m, ρ, v, α)
    """
    ...
    # Arguments
    - 'S': Surface area, m²
    - 'T': Temperature, K
    - 'm': Model, Model
    - 'ρ': Density, kg/m³
    - 'v': Velocity, m/s
    - 'α': Angle of attack, rad
    ...

    Calculates the convective heat rate for a blunt body using ...

    # Returns
    - 'q_conv': Convective heat rate, W/m²
    """

    rn     = m.body.nose_radius + 0.25
    k      = m.planet.k
    q_conv = (k * sqrt(ρ/rn) * v^3) * 1e-4

    return q_conv
end

function heatrate_radiative(S, T, m, ρ, v, α)
    """
    ...
    # Arguments
    - 'S': Surface area, m²
    - 'T': Temperature, K
    - 'm': Model, Model
    - 'ρ': Density, kg/m³
    - 'v': Velocity, m/s
    - 'α': Angle of attack, rad
    ...

    Calculates the radiative heat rate for a blunt body using ...

    # Returns
    - 'q_rad': Convective heat rate, W/m²
    """

    rn = m.body.nose_radius
    C  = 4.736 * 1e4
    b  = 1.22
    fv = [2040, 1780, 1550, 1313, 1065, 850, 660, 495, 359, 238, 151, 115, 81, 55, 35, 19.5, 9.7, 4.3, 1.5, 0]
    vf = [16000,15500, 15000, 14500,14000, 13500, 13000, 12500, 12000, 11500, 11000, 10750, 10500, 10250, 10000, 9750, 9500, 9250, 9000, 0] # m/s

    fn = linear_interpolation(sort(vf), sort(fv)) # check interpolation

    a  = 1.072 * 1e6 * v^(-1.88) * ρ^(-0.325)
    f  = fn(v)

    q_rad = (C * rn^a * ρ^b * f)

    return q_rad
end

function heatrate_convective_radiative(S, T, m, ρ, v, α)
    """
    ...
    # Arguments
    - 'S': Surface area, m²
    - 'T': Temperature, K
    - 'm': Model, Model
    - 'ρ': Density, kg/m³
    - 'v': Velocity, m/s
    - 'α': Angle of attack, rad
    ...

    Calculates the total heatrate as the sum of convective and radiative heat rate for a blunt body.

    # Returns
    - 'q_conv + q_rad': Convective and radiative heat rate, W/m²
    """

    q_conv = heatrate_convective(S, T, m, ρ, v, α)
    q_rad  = heatrate_radiative(S, T, m, ρ, v, α)

    return q_conv # + q_rad
end

function getHeatRate(model::MaxwellianHeat, S::Float64, T::Float64, ρ::Float64, v::Float64, α::Float64)::Float64
    if !model.thermal_contact
        r_prime  = (1 / S^2) * (2*S^2 + 1 - 1 / (1 + sqrt(pi) * S * sin(α) * erf(S * sin(α) * exp((S * sin(α))^2))))
        St_prime = (1 / (4 * sqrt(pi) * S)) * (exp(-(S * sin(α))^2) + (sqrt(pi)) * (S * sin(α)) * erf(S * sin(α)))
    else
        r_prime  = (1 / S^2) * (2*S^2 + 1 - 1 / (1 + sqrt(pi) * S * sin(α) * (1 + erf(S * sin(α))) * exp((S * sin(α))^2)))
        St_prime = (1 / (4 * sqrt(pi) * S)) * (exp(-(S * sin(α))^2) + (sqrt(pi)) * (S * sin(α)) * (1 + erf(S * sin(α))))
    end

    γ = model.planet.γ
    T_0 = T * (1 + ((γ - 1) / γ) * S^2)
    T_r = T + (γ/(γ + 1)) * r_prime * (T_0 - T)
    T_p = T
    T_w = T_p

    # heat_rate = (m.aerodynamics.thermal_accomodation_factor * ρ * m.planet.R * T_p) * 
    #             (sqrt(m.planet.R * T_p / (2 * pi))) * 
    #             ((S^2 + (γ) / (γ - 1) - (γ + 1) / (2 * (γ - 1)) * (T_w / T_p)) * 
    #             (exp(-(S * sin(α))^2) + sqrt(pi) * (S * sin(α)) * 
    #             (1 + erf(S * sin(α)))) - 0.5 * exp(-(S * sin(α))^2)) * 1e-4  # W/cm^2
    
    heat_rate = (model.thermal_accomodation_factor * ρ * model.planet.R * T_p) * 
                (sqrt(model.planet.R * T_p / (2 * pi))) * (
                (S^2 + (γ) / (γ - 1) - (γ + 1) / (2 * (γ - 1)) * (T_w / T_p)) * 
                (exp(-(S * sin(α))^2) + sqrt(pi) * (S * sin(α)) *
                (1 + erf(S * sin(α)))) - 0.5 * exp(-(S * sin(α))^2)) * 1e-4  # W/cm^2

    return heat_rate
end