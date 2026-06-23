module KeplerianPropagator

# AbstractPlanet will be available when loaded through SpaceAGORA.jl
# We accept any planet type with μ, J2, and Rp_e fields

"""
    oe2cart(a, e, inc, Ω, ω, ν, planet)

Convert classical orbital elements to Cartesian state (position and velocity) in the
inertial frame.

# Arguments
- `a::Real`: Semi-major axis [m]
- `e::Real`: Eccentricity
- `inc::Real`: Inclination [rad]
- `Ω::Real`: Right ascension of ascending node [rad]
- `ω::Real`: Argument of periapsis [rad]
- `ν::Real`: True anomaly [rad]
- `planet`: Planet object providing gravitational parameter μ (must have .μ field)

# Returns
- `SVector{6, T}`: Cartesian state [rx, ry, rz, vx, vy, vz] in inertial frame
"""
function oe2cart(a::Real, e::Real, inc::Real, Ω::Real, ω::Real, ν::Real, planet)
    μ = planet.μ
    return _oe2cart_promoted(promote(a, e, inc, Ω, ω, ν, μ)...)
end

function _oe2cart_promoted(a::T, e::T, inc::T, Ω::T, ω::T, ν::T, μ::T) where T
    # Compute eccentric anomaly E from ν
    cosE = (e + cos(ν)) / (1 + e*cos(ν))
    sinE = sqrt(1 - e^2) * sin(ν) / (1 + e*cos(ν))
    
    # Rotation matrix Q columns
    cosΩ, sinΩ = cos(Ω), sin(Ω)
    cosω, sinω = cos(ω), sin(ω)
    cosi, sini = cos(inc), sin(inc)
    
    q1 = T[
        cosω*cosΩ - sinω*cosi*sinΩ,
        cosω*sinΩ + sinω*cosi*cosΩ,
        sinω*sini
    ]
    q2 = T[
        -sinω*cosΩ - cosω*cosi*sinΩ,
        -sinω*sinΩ + cosω*cosi*cosΩ,
        cosω*sini
    ]
    
    # Premultiplied vectors
    u = a * q1
    v = a * sqrt(1 - e^2) * q2
    c = a * e * q1
    
    # Position
    r_eci = c + u * cosE + v * sinE
    
    # Mean motion and anomaly rate
    n = sqrt(μ / a^3)
    E_dot = n / (1 - e*cosE)
    
    # Velocity
    v_eci = E_dot * (v * cosE - u * sinE)
    
    return T[r_eci[1], r_eci[2], r_eci[3], v_eci[1], v_eci[2], v_eci[3]]
end

function _kepler_j2_raan(Ω::Real, a::Real, e::Real, inc::Real, planet, t::Real; use_j2::Bool=false)
    use_j2 || return Ω
    aF = Float64(a)
    eF = Float64(e)
    p = aF * (1.0 - eF^2)
    p > 0.0 || return Ω
    n = sqrt(Float64(planet.μ) / aF^3)
    Ωdot = -1.5 * n * Float64(planet.J2) * (Float64(planet.Rp_e) / p)^2 * cos(Float64(inc))
    return Ω + Ωdot * t
end

"""
    kepler_satellite_state(k, dt, a, e, inc, Ω, ω, ν0, planet; u=0.0, tstep=0.0, use_j2=false)

Propagate a satellite's Keplerian orbit forward in time, optionally including J2 perturbation
and low-thrust acceleration.

# Arguments
- `k::Real`: Time step index (integer multiplier of dt)
- `dt::Real`: Time step size [s]
- `a::Real`: Semi-major axis [m]
- `e::Real`: Eccentricity
- `inc::Real`: Inclination [rad]
- `Ω::Real`: Right ascension of ascending node [rad]
- `ω::Real`: Argument of periapsis [rad]
- `ν0::Real`: Initial true anomaly [rad]
- `planet`: Planet object providing μ, J2, and radius (must have .μ, .J2, .Rp_e fields)
- `u::Real=0.0`: Low-thrust acceleration magnitude [m/s²] (along-track direction)
- `tstep::Real=0.0`: Thrust application time step [s]
- `use_j2::Bool=false`: Whether to include J2 perturbation

# Returns
- `SVector{6, T}`: Cartesian state [rx, ry, rz, vx, vy, vz] at time k*dt
"""
function kepler_satellite_state(k::Real, dt::Real,
    a::Real, e::Real, inc::Real,
    Ω::Real, ω::Real, ν0::Real, planet; u=0.0, tstep=0.0, use_j2::Bool=false)
    μ = planet.μ
    T = promote_type(typeof(k), typeof(dt), typeof(a), typeof(e), typeof(inc), typeof(Ω), typeof(ω), typeof(ν0), typeof(μ), typeof(u), typeof(tstep))
    kT, dtT, aT, eT, incT, ΩT, ωT, ν0T, μT, uT, tstepT = T(k), T(dt), T(a), T(e), T(inc), T(Ω), T(ω), T(ν0), T(μ), T(u), T(tstep)
    n = sqrt(μT / aT^3)
    t = kT * dtT
    Mk = _kepler_true_to_mean_anomaly(ν0T, eT) + n * t
    νk = _kepler_mean_to_true_anomaly(Mk, eT)
    # u is assumed in units of m/s^2, applied only in along track or retrograde directions
    aT += (2 * uT) / sqrt(μT) * aT^(3 / 2) * tstepT # low thrust approximation
    Ωk = _kepler_j2_raan(ΩT, aT, eT, incT, planet, t; use_j2=use_j2)
    return oe2cart(aT, eT, incT, Ωk, ωT, νk, planet)
end

function _kepler_true_to_mean_anomaly(ν::T, e::T) where T
    β = sqrt(max(zero(T), one(T) - e^2))
    E = atan(β * sin(ν), e + cos(ν))
    return E - e * sin(E)
end

function _kepler_mean_to_true_anomaly(M::T, e::T) where T
    E = M
    for _ in 1:30
        E -= (E - e * sin(E) - M) / (one(T) - e * cos(E))
    end
    return 2 * atan(
        sqrt(max(zero(T), one(T) + e)) * sin(E / 2),
        sqrt(max(zero(T), one(T) - e)) * cos(E / 2),
    )
end

"""
    kepler_prop_state(k, dt, oe, planet; u=0.0, tstep=0.0, use_j2=false)

Propagate orbital elements forward in time, optionally including J2 perturbation
and low-thrust acceleration.

# Arguments
- `k::Real`: Time step index (integer multiplier of dt)
- `dt::Real`: Time step size [s]
- `oe::Vector{<:Real}`: Orbital elements [a, e, inc, Ω, ω, ν0]
- `planet`: Planet object providing μ, J2, and radius (must have .μ, .J2, .Rp_e fields)
- `u::Real=0.0`: Low-thrust acceleration magnitude [m/s²]
- `tstep::Real=0.0`: Thrust application time step [s]
- `use_j2::Bool=false`: Whether to include J2 perturbation

# Returns
- `Vector{T}`: Propagated orbital elements [a, e, inc, Ω, ω, ν]
"""
function kepler_prop_state(k::Real, dt::Real, oe::Vector{<:Real}, planet; u=0.0, tstep=0.0, use_j2::Bool=false)
    a, e, inc, Ω, ω, ν0 = oe
    μ = planet.μ
    T = promote_type(typeof(k), typeof(dt), eltype(oe), typeof(μ), typeof(u), typeof(tstep))
    kT, dtT, aT, eT, incT, ΩT, ωT, ν0T, μT, uT, tstepT = T(k), T(dt), T(a), T(e), T(inc), T(Ω), T(ω), T(ν0), T(μ), T(u), T(tstep)
    n = sqrt(μT / aT^3)
    t = kT * dtT
    Mk = _kepler_true_to_mean_anomaly(ν0T, eT) + n * t
    νk = _kepler_mean_to_true_anomaly(Mk, eT)
    # u is assumed in units of m/s^2, applied only in along track or retrograde directions
    aT += (2 * uT) / sqrt(μT) * aT^(3 / 2) * tstepT # low thrust approximation
    Ωk = _kepler_j2_raan(ΩT, aT, eT, incT, planet, t; use_j2=use_j2)
    return [aT, eT, incT, Ωk, ωT, νk]
end

"""
    kepler_client_state(k, dt, a, e, inc, Ω, ω, ν0, planet; use_j2=false)

Propagate a client satellite's Keplerian orbit forward in time (no thrust support),
optionally including J2 perturbation.

# Arguments
- `k::Real`: Time step index (integer multiplier of dt)
- `dt::Real`: Time step size [s]
- `a::Real`: Semi-major axis [m]
- `e::Real`: Eccentricity
- `inc::Real`: Inclination [rad]
- `Ω::Real`: Right ascension of ascending node [rad]
- `ω::Real`: Argument of periapsis [rad]
- `ν0::Real`: Initial true anomaly [rad]
- `planet`: Planet object providing μ, J2, and radius (must have .μ, .J2, .Rp_e fields)
- `use_j2::Bool=false`: Whether to include J2 perturbation

# Returns
- `SVector{6, T}`: Cartesian state [rx, ry, rz, vx, vy, vz] at time k*dt
"""
function kepler_client_state(k::Real, dt::Real,
    a::Real, e::Real, inc::Real,
    Ω::Real, ω::Real, ν0::Real, planet; use_j2::Bool=false)
    μ = planet.μ
    T = promote_type(typeof(k), typeof(dt), typeof(a), typeof(e), typeof(inc), typeof(Ω), typeof(ω), typeof(ν0), typeof(μ))
    kT, dtT, aT, eT, incT, ΩT, ωT, ν0T, μT = T(k), T(dt), T(a), T(e), T(inc), T(Ω), T(ω), T(ν0), T(μ)
    n = sqrt(μT / aT^3)
    t = kT * dtT
    Mk = _kepler_true_to_mean_anomaly(ν0T, eT) + n * t
    νk = _kepler_mean_to_true_anomaly(Mk, eT)
    Ωk = _kepler_j2_raan(ΩT, aT, eT, incT, planet, t; use_j2=use_j2)
    return oe2cart(aT, eT, incT, Ωk, ωT, νk, planet)
end

export oe2cart, kepler_satellite_state, kepler_prop_state, kepler_client_state

end # module KeplerianPropagator
