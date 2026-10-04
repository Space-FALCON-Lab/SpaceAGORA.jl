using CSV
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
    return q_conv
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
    
    heat_rate = (model.thermal_accomodation_factor * ρ * model.planet.R * T_p) * 
                (sqrt(model.planet.R * T_p / (2 * pi))) * (
                (S^2 + (γ) / (γ - 1) - (γ + 1) / (2 * (γ - 1)) * (T_w / T_p)) * 
                (exp(-(S * sin(α))^2) + sqrt(pi) * (S * sin(α)) *
                (1 + erf(S * sin(α)))) - 0.5 * exp(-(S * sin(α))^2)) * 1e-4  # W/cm^2

    return heat_rate
end

"""
    SuttonGravesHeat(; planet, nose_radius_m, k=planet.k)

Sutton–Graves stagnation-point convective heating (Sutton and Graves, NASA TR
R-376, 1971),

    q = k √(ρ / r_n) v³,

returned in W/cm² and applied to every thermal link (panel). `k` is the
atmosphere's Sutton–Graves coefficient in kg^0.5/m and defaults to the planet's
(`planet.k`); `nose_radius_m` is the effective nose radius `r_n`. The free-stream
speed ratio, temperature and incidence passed to [`getHeatRate`](@ref) do not
enter this correlation.
"""
struct SuttonGravesHeat <: AbstractThermalModel
    nose_radius_m::Float64
    k::Float64
    function SuttonGravesHeat(nose_radius_m::Real, k::Real)
        r_f, k_f = Float64(nose_radius_m), Float64(k)
        (isfinite(r_f) && r_f > 0.0) || throw(ArgumentError("SuttonGravesHeat.nose_radius_m must be > 0 m, got $r_f."))
        (isfinite(k_f) && k_f >= 0.0) || throw(ArgumentError("SuttonGravesHeat.k must be >= 0 kg^0.5/m, got $k_f."))
        return new(r_f, k_f)
    end
end
SuttonGravesHeat(; planet::AbstractPlanet, nose_radius_m::Real, k::Real=planet.k) = SuttonGravesHeat(nose_radius_m, k)

function getHeatRate(model::SuttonGravesHeat, S::Float64, T::Float64, ρ::Float64, v::Float64, α::Float64)::Float64
    (ρ > 0.0 && v > 0.0) || return 0.0
    return model.k * sqrt(ρ / model.nose_radius_m) * v^3 * 1e-4  # W/m^2 -> W/cm^2
end

"""
    TabularHeat(velocities_m_s, densities_kg_m3, heat_rates_W_cm2)
    TabularHeat(path::AbstractString)

Vehicle-level heat flux from an aerothermal database tabulated on a velocity ×
density grid. `heat_rates_W_cm2[i, j]` is the flux at `velocities_m_s[i]` and
`densities_kg_m3[j]`; both axes must be strictly increasing, and densities
positive. The flux is interpolated bilinearly in velocity and log density and
held at the grid's edge in velocity and above its highest density. Below the
lowest tabulated density it scales in proportion to density, as free-molecular
heating does, so the rarefied limit goes to zero rather than to the edge value.
The same flux is applied to every thermal link.

The file form reads a CSV with columns `velocity_m_s`, `density_kg_m3` and
`heat_rate_W_cm2`, one row per grid point, rows in any order, every combination
of the listed velocities and densities present exactly once.
"""
struct TabularHeat <: AbstractThermalModel
    velocities::Vector{Float64}
    densities::Vector{Float64}
    log_densities::Vector{Float64}
    heat_rates::Matrix{Float64}
    function TabularHeat(velocities::AbstractVector{<:Real}, densities::AbstractVector{<:Real},
                         heat_rates::AbstractMatrix{<:Real})
        v = Float64.(collect(velocities)); d = Float64.(collect(densities)); q = Float64.(collect(heat_rates))
        (length(v) >= 2 && length(d) >= 2) || throw(ArgumentError("TabularHeat needs at least two velocities and two densities."))
        size(q) == (length(v), length(d)) || throw(ArgumentError(
            "TabularHeat.heat_rates must be $(length(v)) x $(length(d)) (velocities x densities), got $(size(q))."))
        (all(isfinite, v) && issorted(v; lt=<=)) || throw(ArgumentError("TabularHeat velocities must be finite and strictly increasing."))
        (all(x -> isfinite(x) && x > 0.0, d) && issorted(d; lt=<=)) ||
            throw(ArgumentError("TabularHeat densities must be positive and strictly increasing."))
        all(x -> isfinite(x) && x >= 0.0, q) || throw(ArgumentError("TabularHeat heat rates must be finite and >= 0."))
        return new(v, d, log.(d), q)
    end
end

function TabularHeat(path::AbstractString)
    rows = CSV.File(path)
    cols = propertynames(rows)
    for c in (:velocity_m_s, :density_kg_m3, :heat_rate_W_cm2)
        c in cols || throw(ArgumentError("TabularHeat file $(path) has no column $(c)."))
    end
    vs = sort(unique(Float64.(rows.velocity_m_s)))
    ds = sort(unique(Float64.(rows.density_kg_m3)))
    q = fill(NaN, length(vs), length(ds))
    for r in rows
        i = searchsortedfirst(vs, Float64(r.velocity_m_s)); j = searchsortedfirst(ds, Float64(r.density_kg_m3))
        isnan(q[i, j]) || throw(ArgumentError("TabularHeat file $(path) repeats grid point ($(vs[i]), $(ds[j]))."))
        q[i, j] = Float64(r.heat_rate_W_cm2)
    end
    any(isnan, q) && throw(ArgumentError("TabularHeat file $(path) does not cover the full velocity x density grid."))
    return TabularHeat(vs, ds, q)
end

@inline function _tabular_bracket(axis::Vector{Float64}, x::Float64)::Tuple{Int, Float64}
    x <= axis[1] && return 1, 0.0
    x >= axis[end] && return length(axis) - 1, 1.0
    i = searchsortedlast(axis, x)
    return i, (x - axis[i]) / (axis[i + 1] - axis[i])
end

function getHeatRate(model::TabularHeat, S::Float64, T::Float64, ρ::Float64, v::Float64, α::Float64)::Float64
    ρ > 0.0 || return 0.0
    scale = 1.0
    ρ_eval = ρ
    if ρ < model.densities[1]
        scale = ρ / model.densities[1]   # free-molecular: flux proportional to density
        ρ_eval = model.densities[1]
    end
    i, fv = _tabular_bracket(model.velocities, v)
    j, fd = _tabular_bracket(model.log_densities, log(ρ_eval))
    q = model.heat_rates
    @inbounds q_ij = (1 - fv) * (1 - fd) * q[i, j] + fv * (1 - fd) * q[i + 1, j] +
                     (1 - fv) * fd * q[i, j + 1] + fv * fd * q[i + 1, j + 1]
    return scale * q_ij
end
