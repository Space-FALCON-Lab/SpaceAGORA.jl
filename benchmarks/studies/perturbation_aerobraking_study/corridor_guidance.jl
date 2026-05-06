using ComponentArrays
using StaticArrays
using LinearAlgebra

import SpaceAGORA.SimulationModel: calcGuidanceEffect!

# Orbit-by-orbit corridor guidance. At each new orbit the guidance estimates the
# current periapsis altitude from the osculating orbital elements, computes the J2
# ellipsoid surface radius at the sub-periapsis latitude, and schedules a correction
# burn if the predicted altitude falls outside [h_min_m, h_max_m]. The correction
# targets the corridor midpoint. The control layer handles the apoapsis burn timing.

mutable struct CorridorManeuverGuidanceModel <: AbstractGuidanceModel
    h_min_m::Float64
    h_max_m::Float64
    h_mid_m::Float64   # (h_min + h_max) / 2, precomputed
    # Mutable runtime state
    _last_orbit::Int64
    _last_dv::Float64   # m/s, 0 if no correction this orbit
    _last_dir::Float64  # 0.0 = prograde (raise), π = retrograde (lower)
end

function CorridorManeuverGuidanceModel(h_min_m::Float64, h_max_m::Float64)
    return CorridorManeuverGuidanceModel(
        h_min_m, h_max_m, (h_min_m + h_max_m) / 2.0,
        -1, 0.0, 0.0
    )
end

# ---------------------------------------------------------------------------
# Orbital mechanics helpers (self-contained, no dependency on SpaceAGORA internals)
# ---------------------------------------------------------------------------

function _ellipsoid_radius_at_geocentric_lat(Rp_e::Float64, Rp_p::Float64, lat_rad::Float64)::Float64
    ca, sa = cos(lat_rad), sin(lat_rad)
    return Rp_e * Rp_p / sqrt((Rp_p * ca)^2 + (Rp_e * sa)^2)
end

# Returns (a, e, inc_rad, omega_rad) from ECI position/velocity + gravitational parameter.
# Returns (NaN, NaN, NaN, NaN) for degenerate cases.
function _keplerian_elements(pos::SVector{3,Float64}, vel::SVector{3,Float64}, μ::Float64)
    r = norm(pos)
    v = norm(vel)
    if !isfinite(r) || r <= 0.0 || !isfinite(v)
        return NaN, NaN, NaN, NaN
    end

    h_vec = cross(pos, vel)
    h = norm(h_vec)
    if !isfinite(h) || h <= 0.0
        return NaN, NaN, NaN, NaN
    end

    e_vec = cross(vel, h_vec) / μ - pos / r
    e = norm(e_vec)

    ε = v^2 / 2.0 - μ / r
    if ε >= 0.0     # parabolic or hyperbolic
        return NaN, NaN, NaN, NaN
    end
    a = -μ / (2.0 * ε)
    if !isfinite(a) || a <= 0.0
        return NaN, NaN, NaN, NaN
    end

    inc = acos(clamp(h_vec[3] / h, -1.0, 1.0))

    # Ascending node vector: N = ẑ × h
    N_vec = SVector{3,Float64}(-h_vec[2], h_vec[1], 0.0)
    N = norm(N_vec)

    if N > 1e-10 && e > 1e-8
        omega = acos(clamp(dot(N_vec, e_vec) / (N * e), -1.0, 1.0))
        if e_vec[3] < 0.0
            omega = 2π - omega
        end
    else
        omega = 0.0
    end

    return a, e, inc, omega
end

# Periapsis altitude above the J2 ellipsoid.
# lat_peri = arcsin(sin(i) * sin(ω)) — sub-periapsis geocentric latitude.
function _periapsis_altitude_m(a::Float64, e::Float64, inc_rad::Float64, omega_rad::Float64, planet)::Float64
    r_peri = a * (1.0 - e)
    lat_peri = asin(clamp(sin(inc_rad) * sin(omega_rad), -1.0, 1.0))
    R_peri = _ellipsoid_radius_at_geocentric_lat(planet.Rp_e, planet.Rp_p, lat_peri)
    return r_peri - R_peri
end

# ΔV at apoapsis to shift periapsis from current elements to h_target_m altitude.
# Returns (|ΔV| m/s, direction_rad) where direction_rad is 0.0 (prograde) or π (retrograde).
function _corridor_correction_dv(a::Float64, e::Float64, inc_rad::Float64, omega_rad::Float64, planet, h_target_m::Float64)
    μ = planet.μ
    lat_peri = asin(clamp(sin(inc_rad) * sin(omega_rad), -1.0, 1.0))
    R_peri = _ellipsoid_radius_at_geocentric_lat(planet.Rp_e, planet.Rp_p, lat_peri)
    r_apo = a * (1.0 + e)
    r_peri_target = R_peri + h_target_m
    a_target = (r_apo + r_peri_target) / 2.0

    if a_target <= 0.0 || !isfinite(a_target)
        return 0.0, 0.0
    end

    v_current = sqrt(μ * (2.0 / r_apo - 1.0 / a))
    v_target   = sqrt(μ * (2.0 / r_apo - 1.0 / a_target))

    if !isfinite(v_current) || !isfinite(v_target)
        return 0.0, 0.0
    end

    dv = v_target - v_current
    return abs(dv), dv >= 0.0 ? 0.0 : π
end

# ---------------------------------------------------------------------------
# Guidance effector callback — called at guidance_rates Hz by the simulation
# ---------------------------------------------------------------------------

function calcGuidanceEffect!(
    model::CorridorManeuverGuidanceModel,
    u::ComponentVector,
    p,
    ::Float64,
    i::Int64
)::Nothing
    if i < 1 || i > length(p.shared_buffers.maneuver_commands)
        return nothing
    end

    orbit_counter = p.orbit_counter[i]

    if orbit_counter != model._last_orbit
        model._last_orbit = orbit_counter
        model._last_dv    = 0.0
        model._last_dir   = 0.0

        planet = p.args.environment_model.planet
        pos = SVector{3,Float64}(u.sc[i].pos)
        vel = SVector{3,Float64}(u.sc[i].vel)

        a, e, inc, omega = _keplerian_elements(pos, vel, planet.μ)

        if isfinite(a) && isfinite(e) && 0.0 <= e < 1.0
            h_peri = _periapsis_altitude_m(a, e, inc, omega, planet)
            if isfinite(h_peri) && (h_peri < model.h_min_m || h_peri > model.h_max_m)
                dv, dir = _corridor_correction_dv(a, e, inc, omega, planet, model.h_mid_m)
                if isfinite(dv) && dv > 0.0
                    model._last_dv  = dv
                    model._last_dir = dir
                end
            end
        end
    end

    p.shared_buffers.maneuver_commands[i] = PropulsiveManeuverCommand(
        valid=true,
        delta_v_mps=model._last_dv,
        direction_rad=model._last_dir,
        source_orbit=orbit_counter,
    )
    return nothing
end
