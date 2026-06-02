using ComponentArrays

const _APO_TARGET_IDLE = Int64(0)
const _APO_TARGET_COMMAND_ISSUED = Int64(1)
const _APO_TARGET_DISCARDED = Int64(2)

@inline function _ensure_apo_target_state!(guidanceAlg::ApoapsisTargetPeriapsisRaiseGuidanceModel, i::Int64)
    while length(guidanceAlg.command_state) < i
        push!(guidanceAlg.command_state, _APO_TARGET_IDLE)
    end
    return nothing
end

@inline function _wrap_2pi_guidance(θ::Float64)::Float64
    θw = mod(θ, 2pi)
    return θw < 0.0 ? θw + 2pi : θw
end

function _osculating_elements_and_periapsis_direction(pos::SVector{3, Float64}, vel::SVector{3, Float64}, planet)
    rmag = norm(pos)
    vmag2 = dot(vel, vel)
    h = cross(pos, vel)
    hmag = norm(h)
    if !(isfinite(rmag) && rmag > 0.0 && isfinite(vmag2) && isfinite(hmag) && hmag > 0.0)
        return nothing
    end

    specific_energy = 0.5 * vmag2 - planet.μ / rmag
    specific_energy < 0.0 || return nothing
    a = -planet.μ / (2.0 * specific_energy)
    e_vec = cross(vel, h) / planet.μ - pos / rmag
    e = norm(e_vec)
    if !(isfinite(a) && a > 0.0 && isfinite(e) && 0.0 <= e < 1.0)
        return nothing
    end

    true_anomaly = if e > 1e-10
        ν = acos(clamp(dot(e_vec, pos) / (e * rmag), -1.0, 1.0))
        dot(pos, vel) < 0.0 ? 2pi - ν : ν
    else
        0.0
    end
    periapsis_direction = e > 1e-10 ? SVector{3, Float64}(e_vec / e) : SVector{3, Float64}(-pos / rmag)
    return (; a, e, true_anomaly, periapsis_direction)
end

@inline function _oblate_altitude_from_radius(radius::Float64, u_pp::SVector{3, Float64}, planet)::Float64
    x = radius * u_pp[1]
    y = radius * u_pp[2]
    z = radius * u_pp[3]

    f = (planet.Rp_e - planet.Rp_p) / planet.Rp_e
    e2 = 1.0 - (1.0 - f)^2
    ep2 = e2 / (1.0 - e2)
    p_xy = sqrt(x^2 + y^2)
    θ = atan(z * planet.Rp_e, p_xy * planet.Rp_p)
    lat = atan(z + ep2 * planet.Rp_p * sin(θ)^3, p_xy - e2 * planet.Rp_e * cos(θ)^3)
    N = planet.Rp_e / sqrt(1.0 - e2 * sin(lat)^2)
    return p_xy * cos(lat) + (z + e2 * N * sin(lat)^2) * sin(lat) - N
end

@inline function _oblate_surface_radius(u_pp::SVector{3, Float64}, planet)::Float64
    return inv(sqrt((u_pp[1]^2 + u_pp[2]^2) / planet.Rp_e^2 + u_pp[3]^2 / planet.Rp_p^2))
end

function _radius_for_oblate_altitude(target_altitude_m::Float64, u_pp::SVector{3, Float64}, planet)::Float64
    target_altitude_m >= 0.0 || return NaN
    lo = _oblate_surface_radius(u_pp, planet)
    hi = lo + target_altitude_m + abs(planet.Rp_e - planet.Rp_p) + 1.0
    while _oblate_altitude_from_radius(hi, u_pp, planet) < target_altitude_m
        hi += max(target_altitude_m, abs(planet.Rp_e - planet.Rp_p), 1.0)
    end
    for _ in 1:80
        mid = 0.5 * (lo + hi)
        if _oblate_altitude_from_radius(mid, u_pp, planet) < target_altitude_m
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

function calcGuidanceEffect!(guidanceAlg::AerobrakingCampaignPropulsiveManeuverGuidanceModel, u::ComponentVector, p::ODEParams, t::Float64, i::Int64)
    if i < 1 || i > length(p.shared_buffers.maneuver_commands)
        return nothing
    end

    orbit_counter = p.orbit_counter[i]
    maneuver_idx = findfirst(==(orbit_counter), guidanceAlg.maneuver_orbit_number)

    # Campaign guidance currently emits at most one burn command per orbit by
    # looking up the precomputed passage-to-delta-v table.
    command = if maneuver_idx === nothing
        PropulsiveManeuverCommand(valid=true, source_orbit=orbit_counter)
    else
        delta_v_cmd = Float64(guidanceAlg.maneuver_Δv[maneuver_idx])
        PropulsiveManeuverCommand(
            valid=true,
            delta_v_mps=abs(delta_v_cmd),
            direction_rad=(delta_v_cmd > 0.0 ? 0.0 : π),
            source_orbit=orbit_counter
        )
    end

    p.shared_buffers.maneuver_commands[i] = command
    return nothing
end

function calcGuidanceEffect!(guidanceAlg::ApoapsisTargetPeriapsisRaiseGuidanceModel, u::ComponentVector, p::ODEParams, t::Float64, i::Int64)
    if i < 1 || i > length(p.shared_buffers.maneuver_commands)
        return nothing
    end
    _ensure_apo_target_state!(guidanceAlg, i)

    if guidanceAlg.command_state[i] == _APO_TARGET_DISCARDED
        p.shared_buffers.maneuver_commands[i] = PropulsiveManeuverCommand()
        return nothing
    end

    if guidanceAlg.command_state[i] == _APO_TARGET_COMMAND_ISSUED
        plans = p.shared_buffers.maneuver_burn_plans
        if i <= length(plans) && plans[i].valid
            p.shared_buffers.maneuver_commands[i] = PropulsiveManeuverCommand()
            guidanceAlg.command_state[i] = _APO_TARGET_DISCARDED
        end
        return nothing
    end

    pos = SVector{3, Float64}(u.sc[i].pos)
    vel = SVector{3, Float64}(u.sc[i].vel)
    planet = p.args.environment_model.planet
    elements = _osculating_elements_and_periapsis_direction(pos, vel, planet)
    elements === nothing && return nothing

    a = elements.a
    e = elements.e
    current_apoapsis_radius_m = a * (1.0 + e)
    target_apoapsis_radius_m = Float64(guidanceAlg.target_apoapsis_radius_m)
    if !(isfinite(current_apoapsis_radius_m) && isfinite(target_apoapsis_radius_m))
        return nothing
    end
    current_apoapsis_radius_m <= target_apoapsis_radius_m + guidanceAlg.apoapsis_tolerance_m || return nothing

    ν = _wrap_2pi_guidance(Float64(elements.true_anomaly))
    pre_apoapsis = ν <= π
    distance_to_apoapsis = π - ν
    pre_apoapsis && distance_to_apoapsis <= guidanceAlg.apoapsis_window_rad || return nothing

    ephemerides_model = p.args.environment_model.ephemerides_model
    et = p.shared_buffers.et_start[] + Float64(t)
    l_pi = planet_frame_lpi(planet, et, ephemerides_model)
    periapsis_direction_pp = SVector{3, Float64}(l_pi * elements.periapsis_direction)
    periapsis_direction_pp = periapsis_direction_pp / norm(periapsis_direction_pp)
    target_periapsis_radius_m = _radius_for_oblate_altitude(
        Float64(guidanceAlg.target_periapsis_altitude_m),
        periapsis_direction_pp,
        planet
    )
    if !(isfinite(target_periapsis_radius_m) && target_periapsis_radius_m < current_apoapsis_radius_m)
        return nothing
    end

    target_a = 0.5 * (current_apoapsis_radius_m + target_periapsis_radius_m)
    v_current_apoapsis = sqrt(planet.μ * (2.0 / current_apoapsis_radius_m - 1.0 / a))
    v_target_apoapsis = sqrt(planet.μ * (2.0 / current_apoapsis_radius_m - 1.0 / target_a))
    delta_v_mps = v_target_apoapsis - v_current_apoapsis
    if !(isfinite(delta_v_mps) && delta_v_mps > 0.0)
        p.shared_buffers.maneuver_commands[i] = PropulsiveManeuverCommand()
        guidanceAlg.command_state[i] = _APO_TARGET_DISCARDED
        return nothing
    end

    p.shared_buffers.maneuver_commands[i] = PropulsiveManeuverCommand(
        valid=true,
        delta_v_mps=delta_v_mps,
        direction_rad=0.0,
        source_orbit=Int64(p.orbit_counter[i])
    )
    guidanceAlg.command_state[i] = _APO_TARGET_COMMAND_ISSUED
    return nothing
end
