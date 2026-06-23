module PhysicsAdapter

using ..KeplerianPropagator
using LinearAlgebra

"""
    get_planet_from_config(config_dict::AbstractDict)

Create a planet object from configuration dictionary for use with Keplerian propagator.

# Arguments
- `config_dict::AbstractDict`: Configuration dictionary

# Returns
- Planet-like object with μ, J2, and Rp_e fields
"""
function get_planet_from_config(config_dict::AbstractDict)
    phys = config_dict["physical_constants"]
    
    # Create a simple planet structure
    planet = (
        μ = Float64(phys["mu"]),
        J2 = Float64(phys["J2"]),
        Rp_e = Float64(phys["radius"]),
    )
    
    return planet
end

"""
    kepler_satellite_state_adapter(k::Real, dt::Real, a::Real, e::Real, inc::Real, Ω::Real, ω::Real, ν0::Real, planet; u=0.0, tstep=0.0, use_j2::Bool=false)

Adapter wrapper for SpaceAGORA's Keplerian propagator to match CAPOConstellation interface.

# Arguments
- `k::Real`: Time step index
- `dt::Real`: Time step size [s]
- `a::Real`: Semi-major axis [m]
- `e::Real`: Eccentricity
- `inc::Real`: Inclination [rad]
- `Ω::Real`: Right ascension of ascending node [rad]
- `ω::Real`: Argument of periapsis [rad]
- `ν0::Real`: Initial true anomaly [rad]
- `planet`: Planet object with μ, J2, Rp_e fields
- `u::Real=0.0`: Low-thrust acceleration magnitude [m/s²]
- `tstep::Real=0.0`: Thrust application time step [s]
- `use_j2::Bool=false`: Whether to include J2 perturbation

# Returns
- `SVector{6, Float64}`: Cartesian state [rx, ry, rz, vx, vy, vz]
"""
function kepler_satellite_state_adapter(
    k::Real,
    dt::Real,
    a::Real,
    e::Real,
    inc::Real,
    Ω::Real,
    ω::Real,
    ν0::Real,
    planet;
    u = 0.0,
    tstep = 0.0,
    use_j2::Bool = false,
)
    return kepler_satellite_state(k, dt, a, e, inc, Ω, ω, ν0, planet; u=u, tstep=tstep, use_j2=use_j2)
end

"""
    kepler_client_state_adapter(k::Real, dt::Real, a::Real, e::Real, inc::Real, Ω::Real, ω::Real, ν0::Real, planet; use_j2::Bool=false)

Adapter wrapper for client satellite propagation (no thrust support).

# Arguments
- `k::Real`: Time step index
- `dt::Real`: Time step size [s]
- `a::Real`: Semi-major axis [m]
- `e::Real`: Eccentricity
- `inc::Real`: Inclination [rad]
- `Ω::Real`: Right ascension of ascending node [rad]
- `ω::Real`: Argument of periapsis [rad]
- `ν0::Real`: Initial true anomaly [rad]
- `planet`: Planet object with μ, J2, Rp_e fields
- `use_j2::Bool=false`: Whether to include J2 perturbation

# Returns
- `SVector{6, Float64}`: Cartesian state [rx, ry, rz, vx, vy, vz]
"""
function kepler_client_state_adapter(
    k::Real,
    dt::Real,
    a::Real,
    e::Real,
    inc::Real,
    Ω::Real,
    ω::Real,
    ν0::Real,
    planet;
    use_j2::Bool = false,
)
    return kepler_client_state(k, dt, a, e, inc, Ω, ω, ν0, planet; use_j2=use_j2)
end

"""
    state_to_cartesian(oe::AbstractVector{<:Real}, planet) -> Vector{Float64}

Convert orbital elements to Cartesian state.

# Arguments
- `oe::AbstractVector{<:Real}`: Orbital elements [a, e, inc, Ω, ω, ν]
- `planet`: Planet object with μ field

# Returns
- `Vector{Float64}`: Cartesian state [rx, ry, rz, vx, vy, vz]
"""
function state_to_cartesian(oe::AbstractVector{<:Real}, planet)
    a, e, inc, Ω, ω, ν = oe
    return oe2cart(a, e, inc, Ω, ω, ν, planet)
end

"""
    cartesian_to_state(state::AbstractVector{<:Real}, μ::Real) -> Vector{Float64}

Convert Cartesian state to orbital elements.

# Arguments
- `state::AbstractVector{<:Real}`: Cartesian state [rx, ry, rz, vx, vy, vz]
- `μ::Real`: Gravitational parameter

# Returns
- `Vector{Float64}`: Orbital elements [a, e, inc, Ω, ω, ν]
"""
function cartesian_to_state(state::AbstractVector{<:Real}, μ::Real)
    r = state[1:3]
    v = state[4:6]
    
    h = cross(r, v)
    h_norm = norm(h)
    
    e_vec = (cross(v, h) / μ) - (r / norm(r))
    e = norm(e_vec)
    
    r_norm = norm(r)
    v_norm = norm(v)
    
    # Specific orbital energy
    E = (v_norm^2 / 2) - (μ / r_norm)
    a = -μ / (2E)
    
    # Inclination
    inc = acos(h[3] / h_norm)
    
    # Right ascension of ascending node
    n = cross([0, 0, 1], h)
    n_norm = norm(n)
    if n_norm > 1e-12
        Ω = acos(n[1] / n_norm)
        if n[2] < 0
            Ω = 2π - Ω
        end
    else
        Ω = 0.0
    end
    
    # Argument of periapsis
    if e > 1e-12
        if n_norm > 1e-12
            ω = acos(dot(n, e_vec) / (n_norm * e))
            if e_vec[3] < 0
                ω = 2π - ω
            end
        else
            ω = 0.0
        end
    else
        ω = 0.0
    end
    
    # True anomaly
    if e > 1e-12
        ν = acos(dot(e_vec, r) / (e * r_norm))
        if dot(r, v) < 0
            ν = 2π - ν
        end
    else
        ν = 0.0
    end
    
    return [a, e, inc, Ω, ω, ν]
end

export get_planet_from_config, kepler_satellite_state_adapter, kepler_client_state_adapter,
       state_to_cartesian, cartesian_to_state

end # module PhysicsAdapter
