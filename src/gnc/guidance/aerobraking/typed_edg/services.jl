# Private SpaceAGORA projections and atmosphere/frame services for typed EDG.
# The distinct query conventions here are intentional compatibility contracts.
module EDGServices
using ..ConfigTypes: ODEParams
using ..EnvironmentModels: getDensity, _environment_wind
using ..EphemeridesModels: planet_frame_lpi
using ..FrameTransforms: r_intor_p!, rtolatlong, latlongtoNED, rvtoorbitalelement
using LinearAlgebra, StaticArrays

@inline function _edg_control_sat_state(u, i::Int)
    return hasproperty(u, :sc) ? u.sc[i] : u
end

@inline function _edg_control_pos_vel_mass(sc)
    pos = hasproperty(sc, :pos) ? SVector{3, Float64}(sc.pos) : SVector{3, Float64}(sc[1], sc[2], sc[3])
    vel = hasproperty(sc, :vel) ? SVector{3, Float64}(sc.vel) : SVector{3, Float64}(sc[4], sc[5], sc[6])
    mass = hasproperty(sc, :mass) ? Float64(sc.mass) : (length(sc) >= 7 ? Float64(sc[7]) : NaN)
    return pos, vel, mass
end

function _edg_environment_state(u, p::ODEParams, t::Float64, i::Int)
    sc = _edg_control_sat_state(u, i)
    pos, vel, _ = _edg_control_pos_vel_mass(sc)
    args = p.args
    planet = args.environment_model.planet
    et = _edg_ephemeris_time(p, t)
    pos_pp, vel_pp = r_intor_p!(pos, vel, planet, et, args.environment_model.ephemerides_model)
    lla = rtolatlong(pos_pp, planet)
    rho, temperature, wind = getDensity(args.environment_model.density_model, lla[1], lla[2], lla[3], t, args.environment_model.wind, p)
    wind = _environment_wind(args.environment_model.wind, wind)
    uD, uN, uE = latlongtoNED(lla)
    wE, wN, wU = wind
    wind_pp = wN * uN + wE * uE - wU * uD
    # `wind` is the air's own velocity relative to the rotating planet in local
    # east/north/up components (GRAM ewWind/nsWind/verticalWind, "Eastward
    # Wind" positive toward east), and `vel_pp` is the spacecraft's velocity
    # relative to the rotating planet, so the airspeed is their difference. The
    # aerodynamic wrench and the T-EDG predictor use the same convention.
    vel_pp_rw = vel_pp - wind_pp
    speed = norm(vel_pp_rw)
    sound_speed = sqrt(max(0.0, planet.γ * planet.R * temperature))
    molecular_speed_ratio = sound_speed > 0.0 ? sqrt(0.5 * planet.γ) * speed / sound_speed : 0.0
    dynamic_pressure = 0.5 * max(0.0, rho) * speed^2
    return (
        altitude_m=Float64(lla[1]),
        rho=max(0.0, rho),
        temperature=max(temperature, eps(Float64)),
        speed=speed,
        molecular_speed_ratio=molecular_speed_ratio,
        dynamic_pressure=dynamic_pressure,
    )
end

@inline function _edg_in_drag_passage(p::ODEParams, env)::Bool
    ei_m = 1e3 * Float64(p.args.environment_model.EI)
    return isfinite(env.altitude_m) && isfinite(ei_m) && env.altitude_m <= ei_m
end

@inline function _edg_ephemeris_time(p::ODEParams, t_abs::Float64)::Float64
    if hasproperty(p, :shared_buffers) && hasproperty(p.shared_buffers, :et_start)
        return p.shared_buffers.et_start[] + t_abs
    end
    return t_abs
end

function _edg_planet_frame_lpi(p::ODEParams, t_abs::Float64)
    planet = p.args.environment_model.planet
    ephemerides_model = p.args.environment_model.ephemerides_model
    return planet_frame_lpi(planet, _edg_ephemeris_time(p, t_abs), ephemerides_model)
end

function _edg_targeting_prediction_environment(p::ODEParams, r::SVector{3, Float64}, v::SVector{3, Float64}, t_abs::Float64)
    planet = p.args.environment_model.planet
    ephemerides_model = p.args.environment_model.ephemerides_model
    et = _edg_ephemeris_time(p, t_abs)
    l_pi = _edg_planet_frame_lpi(p, t_abs)
    pos_pp, vel_pp = r_intor_p!(r, v, planet, et, ephemerides_model)
    lla = rtolatlong(pos_pp, planet)
    rho, temperature, wind = getDensity(
        p.args.environment_model.density_model,
        lla[1],
        lla[2],
        lla[3],
        t_abs,
        p.args.environment_model.wind,
        p,
    )
    wind = _environment_wind(p.args.environment_model.wind, wind)
    uD, uN, uE = latlongtoNED(lla)
    wE, wN, wU = wind
    wind_pp = wN * uN + wE * uE - wU * uD
    # Airspeed: spacecraft velocity minus the air's velocity (see
    # `_edg_environment_state`).
    vel_pp_rw = vel_pp - wind_pp
    speed = norm(vel_pp_rw)
    sound_speed = sqrt(max(0.0, planet.γ * planet.R * temperature))
    speed_ratio = sound_speed > 0.0 ? sqrt(0.5 * planet.γ) * speed / sound_speed : 0.0
    return (
        l_pi=l_pi,
        pos_pp=pos_pp,
        vel_pp=vel_pp,
        vel_pp_rw=vel_pp_rw,
        altitude_m=Float64(lla[1]),
        rho=max(0.0, rho),
        temperature=max(temperature, eps(Float64)),
        speed=speed,
        molecular_speed_ratio=speed_ratio,
        dynamic_pressure=0.5 * max(0.0, rho) * speed^2,
    )
end

function _edg_sample_prediction_atmosphere(p::ODEParams, altitude::Float64, t_abs::Float64)
    rho, temperature, _ = getDensity(
        p.args.environment_model.density_model,
        max(0.0, altitude),
        0.0,
        0.0,
        t_abs,
        p.args.environment_model.wind,
        p,
    )
    return max(0.0, rho), max(temperature, eps(Float64))
end

end # module EDGServices
