using Roots

"""
    SolarPanelAngleOfAttackControlModel(; controlled_panel_links=(2, 3))

Control effector that articulates the listed solar-panel links to realize a
commanded angle of attack during aerobraking passes. Link indices must be
positive and at least one link is required.
"""
struct SolarPanelAngleOfAttackControlModel <: AbstractControlEffectorModel
    controlled_panel_links::Tuple{Vararg{Int}}
end

@inline function _edg_panel_link_tuple(controlled_panel_links)
    links = tuple((Int(idx) for idx in controlled_panel_links)...)
    isempty(links) && throw(ArgumentError("SolarPanelAngleOfAttackControlModel requires at least one controlled link."))
    any(<=(0), links) && throw(ArgumentError("Controlled panel link indices must be positive."))
    return links
end

function SolarPanelAngleOfAttackControlModel(; controlled_panel_links=(2, 3))
    links = _edg_panel_link_tuple(controlled_panel_links)
    return SolarPanelAngleOfAttackControlModel(links)
end

"""
    AerobrakingEnergyDepletionControlModel

Control-side companion of the energy-depletion guidance strategy: tracks the
guidance-selected mode and drives the panel angle-of-attack effector under
the configured heat and structural limits. Holds the shared
[`AerobrakingEnergyDepletionConfig`](@ref) / [`AerobrakingEnergyDepletionState`](@ref)
and a [`SolarPanelAngleOfAttackControlModel`](@ref).
"""
struct AerobrakingEnergyDepletionControlModel <: AbstractControlEffectorModel
    config::AerobrakingEnergyDepletionConfig
    state::AerobrakingEnergyDepletionState
    aoa_effector::SolarPanelAngleOfAttackControlModel
end

function AerobrakingEnergyDepletionControlModel(
    config::AerobrakingEnergyDepletionConfig,
    state::AerobrakingEnergyDepletionState;
    aoa_effector::SolarPanelAngleOfAttackControlModel=SolarPanelAngleOfAttackControlModel(config.controlled_panel_links),
)
    return AerobrakingEnergyDepletionControlModel(config, state, aoa_effector)
end

@inline function _edg_control_state_index_ok(state::AerobrakingEnergyDepletionState, i::Int)::Bool
    return 1 <= i <= length(state.selected_mode)
end

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
    uD, uN, uE = latlongtoNED(lla)
    wE, wN, wU = wind
    wind_pp = wN * uN + wE * uE - wU * uD
    vel_pp_rw = vel_pp + wind_pp
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

function _edg_recompute_switches!(
    model::AerobrakingEnergyDepletionControlModel,
    p::ODEParams,
    env,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass::Float64,
    heat_load_j_cm2::Float64,
    t::Float64,
    i::Int,
)
    config = model.config
    state = model.state
    if state.selected_mode[i] == :targeting && state.targeting_active[i]
        isfinite(state.targeting_switch_s[i]) && return nothing
        _edg_in_drag_passage(p, env) || return nothing
        state.targeting_switch_s[i] = _edg_solve_targeting_switch(
            config,
            state,
            p,
            spacecraft,
            pos,
            vel,
            mass,
            t,
            i;
            heat_load_j_cm2=heat_load_j_cm2,
            heat_rate_control=(:heat_rate in config.max_energy_submodes),
            structural_control=(:structural_load in config.max_energy_submodes),
        )
        state.last_switch_solve_t[i] = t
    elseif state.selected_mode[i] == :max_energy_depletion && (:heat_load in config.max_energy_submodes)
        if !_edg_in_drag_passage(p, env)
            if state.heat_load_drag_passage_active[i]
                state.heat_load_switch_solved[i] = false
                state.heat_load_drag_passage_active[i] = false
            end
            return nothing
        end
        state.heat_load_drag_passage_active[i] = true
        state.heat_load_switch_solved[i] && return nothing
        state.heat_load_switches_s[i] = _edg_solve_heat_load_switches(
            config,
            p,
            spacecraft,
            pos,
            vel,
            mass,
            env,
            heat_load_j_cm2,
            t;
            heat_rate_control=(:heat_rate in config.max_energy_submodes),
            structural_control=(:structural_load in config.max_energy_submodes),
        )
        state.heat_load_switch_solved[i] = true
        state.last_switch_solve_t[i] = t
    end
    return nothing
end

@inline function _edg_base_alpha(model::AerobrakingEnergyDepletionControlModel, t::Float64, i::Int)::Float64
    config = model.config
    state = model.state
    mode = state.selected_mode[i]
    if mode == :safe_low_drag || state.safe_low_drag[i]
        return config.min_alpha_rad
    elseif mode == :targeting && state.targeting_active[i]
        return t >= state.targeting_switch_s[i] ? config.min_alpha_rad : config.max_alpha_rad
    elseif mode == :max_energy_depletion
        if _edg_heat_load_low_drag_active(model, t, i)
            return config.min_alpha_rad
        end
    end
    return config.max_alpha_rad
end

function _edg_command_alpha!(
    model::AerobrakingEnergyDepletionControlModel,
    p::ODEParams,
    u,
    env,
    spacecraft,
    base_alpha::Float64,
    heat_load_j_cm2::Float64,
    heat_load_low_drag_active::Bool,
    i::Int,
)
    config = model.config
    state = model.state
    alpha = clamp(base_alpha, config.min_alpha_rad, config.max_alpha_rad)
    alpha_hr = alpha
    alpha_struct = alpha
    heat_rate_limit = config.heat_rate_limit_w_cm2
    heat_rate_active = (:heat_rate in config.max_energy_submodes) &&
        !heat_load_low_drag_active
    if heat_rate_active
        alpha_past = isfinite(state.last_alpha_rad[i]) ? state.last_alpha_rad[i] : alpha
        alpha_hr = _edg_heat_rate_alpha(config, p, env, alpha; limit_override=heat_rate_limit, alpha_past=alpha_past)
    end
    structural_active = (:structural_load in config.max_energy_submodes) &&
        !heat_load_low_drag_active
    if structural_active
        alpha_struct = _edg_structural_alpha(
            config,
            p,
            env,
            spacecraft,
            model.aoa_effector.controlled_panel_links,
            alpha,
        )
    end
    if heat_rate_active &&
            structural_active
        alpha = min(alpha_hr, alpha_struct)
    elseif heat_rate_active
        alpha = alpha_hr
    elseif structural_active
        alpha = alpha_struct
    end
    if (:heat_load in config.max_energy_submodes) && heat_load_j_cm2 >= config.heat_load_limit_j_cm2
        alpha = min(alpha, config.min_alpha_rad)
    end
    alpha = clamp(alpha, config.min_alpha_rad, config.max_alpha_rad)
    state.last_alpha_rad[i] = alpha
    state.last_alpha_heat_rate_rad[i] = alpha_hr
    state.last_alpha_structural_rad[i] = alpha_struct
    state.last_heat_rate_w_cm2[i] = _edg_maxwellian_heat_rate(p, env, alpha)
    state.last_heat_load_j_cm2[i] = heat_load_j_cm2
    state.last_dynamic_pressure_pa[i] = env.dynamic_pressure
    return alpha
end

function _apply_solar_panel_aoa!(
    effector::SolarPanelAngleOfAttackControlModel,
    spacecraft,
    alpha::Float64,
)
    links = spacecraft.links
    for idx in effector.controlled_panel_links
        1 <= idx <= length(links) || throw(ArgumentError("Controlled panel link index $(idx) is out of bounds for spacecraft with $(length(links)) links."))
        link = links[idx]
        link.root && throw(ArgumentError("Controlled panel link $(idx) is a root link and cannot be rotated."))
        axis = SVector{3, Float64}(abs.(link.r))
        rotate_link(link, axis, pi / 2 - alpha)
        link.α = alpha
    end
    return nothing
end

function calcControlEffect!(
    model::AerobrakingEnergyDepletionControlModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int64,
)
    state = model.state
    _edg_control_state_index_ok(state, i) || return nothing
    config = model.config
    if state.selected_mode[i] == :inactive
        state.selected_mode[i] = (:max_energy_depletion in config.guidance_modes) ? :max_energy_depletion : :safe_low_drag
    end
    sc = _edg_control_sat_state(u, i)
    env = _edg_environment_state(u, p, Float64(t), i)
    spacecraft = p.args.dynamics_model.spacecraft[i]
    pos, vel, mass = _edg_control_pos_vel_mass(sc)
    heat_load = _edg_max_heat_load_for_links(sc, model.aoa_effector.controlled_panel_links)
    _edg_recompute_switches!(model, p, env, spacecraft, pos, vel, mass, heat_load, Float64(t), i)
    heat_load_low_drag_active = _edg_heat_load_low_drag_active(model, Float64(t), i)
    base_alpha = _edg_base_alpha(model, Float64(t), i)
    alpha = _edg_command_alpha!(model, p, u, env, spacecraft, base_alpha, heat_load, heat_load_low_drag_active, i)
    _apply_solar_panel_aoa!(model.aoa_effector, spacecraft, alpha)
    return nothing
end

function calcControlForceTorque(
    model::AerobrakingEnergyDepletionControlModel,
    u::AbstractVector,
    p::ODEParams,
    i::Int64,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function calcControlForceTorque(
    model::SolarPanelAngleOfAttackControlModel,
    u::AbstractVector,
    p::ODEParams,
    i::Int64,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function _edg_orbit_metrics_from_rv(pos::SVector{3, Float64}, vel::SVector{3, Float64}, mass::Float64, planet)
    radius = norm(pos)
    energy = 0.5 * dot(vel, vel) - planet.μ / radius
    if !(isfinite(energy) && energy < 0.0)
        return (energy=energy, periapsis=NaN, apoapsis=Inf)
    end

    oe = rvtoorbitalelement(pos, vel, mass, planet)
    a, e = oe[1], oe[2]
    if !(isfinite(a) && isfinite(e))
        return (energy=energy, periapsis=NaN, apoapsis=NaN)
    end
    return (energy=energy, periapsis=a * (1.0 - e), apoapsis=a * (1.0 + e))
end

function _edg_target_energy_from_apoapsis(planet, target_apoapsis_radius_m::Float64, periapsis_radius_m::Float64)
    if !(isfinite(target_apoapsis_radius_m) && target_apoapsis_radius_m > 0.0 &&
            isfinite(periapsis_radius_m) && periapsis_radius_m > 0.0)
        return NaN
    end
    return -planet.μ / (target_apoapsis_radius_m + periapsis_radius_m)
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
    uD, uN, uE = latlongtoNED(lla)
    wE, wN, wU = wind
    wind_pp = wN * uN + wE * uE - wU * uD
    vel_pp_rw = vel_pp + wind_pp
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

function _edg_targeting_constrained_alpha(
    config::AerobrakingEnergyDepletionConfig,
    p::ODEParams,
    spacecraft,
    env,
    base_alpha::Float64,
    alpha_past::Float64;
    heat_rate_control::Bool,
    structural_control::Bool,
)::Float64
    base = clamp(base_alpha, config.min_alpha_rad, config.max_alpha_rad)
    base <= config.min_alpha_rad + 1e-8 && return config.min_alpha_rad

    alpha_hr = base
    alpha_struct = base
    if heat_rate_control
        thermal_model = p.args.environment_model.thermal_model
        taf = hasproperty(thermal_model, :thermal_accomodation_factor) ? Float64(thermal_model.thermal_accomodation_factor) : 1.0
        planet = p.args.environment_model.planet
        alpha_hr = _energy_depletion_heatrate_root_alpha(
            taf=taf,
            rho=env.rho,
            T_p=env.temperature,
            R=planet.R,
            gamma=planet.γ,
            S=env.molecular_speed_ratio,
            max_alpha=base,
            min_alpha=config.min_alpha_rad,
            heat_rate_limit=config.heat_rate_limit_w_cm2,
            alpha_past=alpha_past,
        )
    end
    if structural_control
        alpha_struct = _energy_depletion_struct_load_root_alpha(
            config,
            env,
            spacecraft,
            config.controlled_panel_links,
            base,
        )
    end
    return clamp(min(alpha_hr, alpha_struct), config.min_alpha_rad, config.max_alpha_rad)
end

function _edg_targeting_prediction_time_grid(duration::Float64)
    step = 0.1
    n = clamp(ceil(Int, duration / step) + 1, 64, 20_000)
    return collect(range(0.0, duration; length=n))
end

function _edg_targeting_aero_acceleration(
    config::AerobrakingEnergyDepletionConfig,
    p::ODEParams,
    spacecraft,
    r::SVector{3, Float64},
    v::SVector{3, Float64},
    mass::Float64,
    t_abs::Float64,
    alpha::Float64,
)::SVector{3, Float64}
    planet = p.args.environment_model.planet
    env = _edg_targeting_prediction_environment(p, r, v, t_abs)
    if !(env.rho > 0.0 && env.speed > eps(Float64))
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    h_pp = cross(env.pos_pp, env.vel_pp)
    h_norm = norm(h_pp)
    h_norm <= eps(Float64) && return SVector{3, Float64}(0.0, 0.0, 0.0)
    vel_hat = env.vel_pp_rw / env.speed
    lift_vec = cross(h_pp / h_norm, vel_hat)
    lift_norm = norm(lift_vec)
    lift_hat = lift_norm > eps(Float64) ? lift_vec / lift_norm : SVector{3, Float64}(0.0, 0.0, 0.0)
    drag_hat = -vel_hat
    q = env.dynamic_pressure
    controlled = Set{Int}(config.controlled_panel_links)
    aero_module = getfield(getfield(parentmodule(@__MODULE__), :DynamicEffectors), :AerodynamicEffectors)
    aero_coeff = getfield(aero_module, :aerodynamic_coefficient_fM)

    force_pp = MVector{3, Float64}(0.0, 0.0, 0.0)
    for (idx, link) in pairs(spacecraft.links)
        area = max(0.0, Float64(link.ref_area))
        area == 0.0 && continue
        link_alpha = link.root ? (pi / 2) :
            (idx in controlled ? alpha : clamp(Float64(link.α), config.min_alpha_rad, config.max_alpha_rad))
        coeffs = aero_coeff(link, env.temperature, max(env.molecular_speed_ratio, eps(Float64)), link_alpha, Float64(link.β), Float64(link.θ))
        cl = Float64(coeffs[1])
        cd = max(0.0, Float64(coeffs[2]))
        force_pp .+= q * area * (cd * drag_hat + cl * lift_hat)
    end

    return SVector{3, Float64}(env.l_pi' * SVector{3, Float64}(force_pp)) / mass
end

function _edg_integrated_targeting_trajectory(
    config::AerobrakingEnergyDepletionConfig,
    p::ODEParams,
    spacecraft,
    pos0::SVector{3, Float64},
    vel0::SVector{3, Float64},
    mass::Float64,
    t::Float64,
    times::Vector{Float64},
    switch_time_s::Float64;
    heat_rate_control::Bool,
    structural_control::Bool,
)
    planet = p.args.environment_model.planet
    n = length(times)
    positions = Vector{SVector{3, Float64}}(undef, n)
    velocities = Vector{SVector{3, Float64}}(undef, n)
    alpha_profile = Vector{Float64}(undef, n)
    positions[1] = pos0
    velocities[1] = vel0
    alpha_past = config.max_alpha_rad

    function acceleration(r, v, tau, alpha)
        radius = max(norm(r), eps(Float64))
        gravity = -planet.μ * r / radius^3
        aero = _edg_targeting_aero_acceleration(config, p, spacecraft, r, v, mass, t + tau, alpha)
        return gravity + aero
    end

    for j in 1:(n - 1)
        tau = times[j]
        r = positions[j]
        v = velocities[j]
        base_alpha = (t + tau) >= switch_time_s ? config.min_alpha_rad : config.max_alpha_rad
        env = _edg_targeting_prediction_environment(p, r, v, t + tau)
        alpha = _edg_targeting_constrained_alpha(
            config,
            p,
            spacecraft,
            env,
            base_alpha,
            alpha_past;
            heat_rate_control=heat_rate_control,
            structural_control=structural_control,
        )
        alpha_profile[j] = alpha
        alpha_past = alpha

        dt = times[j + 1] - times[j]
        a1 = acceleration(r, v, tau, alpha)
        k1r, k1v = v, a1
        a2 = acceleration(r + 0.5 * dt * k1r, v + 0.5 * dt * k1v, tau + 0.5 * dt, alpha)
        k2r, k2v = v + 0.5 * dt * k1v, a2
        a3 = acceleration(r + 0.5 * dt * k2r, v + 0.5 * dt * k2v, tau + 0.5 * dt, alpha)
        k3r, k3v = v + 0.5 * dt * k2v, a3
        a4 = acceleration(r + dt * k3r, v + dt * k3v, times[j + 1], alpha)
        k4r, k4v = v + dt * k3v, a4
        positions[j + 1] = r + dt * (k1r + 2.0 * k2r + 2.0 * k3r + k4r) / 6.0
        velocities[j + 1] = v + dt * (k1v + 2.0 * k2v + 2.0 * k3v + k4v) / 6.0
    end

    final_tau = times[end]
    final_base = (t + final_tau) >= switch_time_s ? config.min_alpha_rad : config.max_alpha_rad
    final_env = _edg_targeting_prediction_environment(p, positions[end], velocities[end], t + final_tau)
    alpha_profile[end] = _edg_targeting_constrained_alpha(
        config,
        p,
        spacecraft,
        final_env,
        final_base,
        alpha_past;
        heat_rate_control=heat_rate_control,
        structural_control=structural_control,
    )

    h = Vector{Float64}(undef, n)
    gamma = Vector{Float64}(undef, n)
    speed = Vector{Float64}(undef, n)
    rho = Vector{Float64}(undef, n)
    temperature = Vector{Float64}(undef, n)
    speed_ratio = Vector{Float64}(undef, n)
    for j in eachindex(times)
        env = _edg_targeting_prediction_environment(p, positions[j], velocities[j], t + times[j])
        radius = max(norm(positions[j]), eps(Float64))
        vel_norm = max(norm(velocities[j]), eps(Float64))
        h[j] = env.altitude_m
        gamma[j] = asin(clamp(dot(positions[j], velocities[j]) / max(radius * vel_norm, eps(Float64)), -1.0, 1.0))
        speed[j] = env.speed
        rho[j] = env.rho
        temperature[j] = env.temperature
        speed_ratio[j] = env.molecular_speed_ratio
    end
    return (
        time=times,
        h=h,
        gamma=gamma,
        speed=speed,
        rho=rho,
        temperature=temperature,
        speed_ratio=speed_ratio,
        positions=positions,
        velocities=velocities,
        alpha_profile=alpha_profile,
    )
end

function _edg_integrated_max_energy_depletion_trajectory(
    config::AerobrakingEnergyDepletionConfig,
    p::ODEParams,
    spacecraft,
    pos0::SVector{3, Float64},
    vel0::SVector{3, Float64},
    mass::Float64,
    t::Float64,
    times::Vector{Float64},
    heat_load_switches::NTuple{2, Float64};
    heat_rate_control::Bool,
    structural_control::Bool,
)
    planet = p.args.environment_model.planet
    n = length(times)
    positions = Vector{SVector{3, Float64}}(undef, n)
    velocities = Vector{SVector{3, Float64}}(undef, n)
    alpha_profile = Vector{Float64}(undef, n)
    positions[1] = pos0
    velocities[1] = vel0
    alpha_past = config.max_alpha_rad

    function acceleration(r, v, tau, alpha)
        radius = max(norm(r), eps(Float64))
        gravity = -planet.μ * r / radius^3
        aero = _edg_targeting_aero_acceleration(config, p, spacecraft, r, v, mass, t + tau, alpha)
        return gravity + aero
    end

    function max_energy_alpha(r, v, tau)
        t_abs = t + tau
        heat_load_low_drag = (:heat_load in config.max_energy_submodes) &&
            isfinite(heat_load_switches[1]) &&
            isfinite(heat_load_switches[2]) &&
            heat_load_switches[1] <= t_abs <= heat_load_switches[2]
        heat_load_low_drag && return config.min_alpha_rad

        base_alpha = config.max_alpha_rad
        env = _edg_targeting_prediction_environment(p, r, v, t_abs)
        alpha = _edg_targeting_constrained_alpha(
            config,
            p,
            spacecraft,
            env,
            base_alpha,
            alpha_past;
            heat_rate_control=heat_rate_control,
            structural_control=structural_control,
        )
        return alpha
    end

    for j in 1:(n - 1)
        tau = times[j]
        r = positions[j]
        v = velocities[j]
        alpha = max_energy_alpha(r, v, tau)
        alpha_profile[j] = alpha
        alpha_past = alpha

        dt = times[j + 1] - times[j]
        a1 = acceleration(r, v, tau, alpha)
        k1r, k1v = v, a1
        a2 = acceleration(r + 0.5 * dt * k1r, v + 0.5 * dt * k1v, tau + 0.5 * dt, alpha)
        k2r, k2v = v + 0.5 * dt * k1v, a2
        a3 = acceleration(r + 0.5 * dt * k2r, v + 0.5 * dt * k2v, tau + 0.5 * dt, alpha)
        k3r, k3v = v + 0.5 * dt * k2v, a3
        a4 = acceleration(r + dt * k3r, v + dt * k3v, times[j + 1], alpha)
        k4r, k4v = v + dt * k3v, a4
        positions[j + 1] = r + dt * (k1r + 2.0 * k2r + 2.0 * k3r + k4r) / 6.0
        velocities[j + 1] = v + dt * (k1v + 2.0 * k2v + 2.0 * k3v + k4v) / 6.0
    end

    final_tau = times[end]
    alpha_profile[end] = max_energy_alpha(positions[end], velocities[end], final_tau)

    h = Vector{Float64}(undef, n)
    gamma = Vector{Float64}(undef, n)
    speed = Vector{Float64}(undef, n)
    rho = Vector{Float64}(undef, n)
    temperature = Vector{Float64}(undef, n)
    speed_ratio = Vector{Float64}(undef, n)
    for j in eachindex(times)
        env = _edg_targeting_prediction_environment(p, positions[j], velocities[j], t + times[j])
        radius = max(norm(positions[j]), eps(Float64))
        vel_norm = max(norm(velocities[j]), eps(Float64))
        h[j] = env.altitude_m
        gamma[j] = asin(clamp(dot(positions[j], velocities[j]) / max(radius * vel_norm, eps(Float64)), -1.0, 1.0))
        speed[j] = env.speed
        rho[j] = env.rho
        temperature[j] = env.temperature
        speed_ratio[j] = env.molecular_speed_ratio
    end
    return (
        time=times,
        h=h,
        gamma=gamma,
        speed=speed,
        rho=rho,
        temperature=temperature,
        speed_ratio=speed_ratio,
        positions=positions,
        velocities=velocities,
        alpha_profile=alpha_profile,
    )
end

function _edg_predict_targeting_outcome(
    config::AerobrakingEnergyDepletionConfig,
    p::ODEParams,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass_state::Float64,
    t::Float64,
    switch_time_s::Float64;
    heat_rate_control::Bool,
    structural_control::Bool,
)
    mass = _edg_predict_mass(spacecraft, mass_state)
    duration = _edg_drag_passage_duration(config, p, pos, vel, mass)
    times = _edg_targeting_prediction_time_grid(duration)
    track = _edg_integrated_targeting_trajectory(
        config,
        p,
        spacecraft,
        pos,
        vel,
        mass,
        t,
        times,
        switch_time_s;
        heat_rate_control=heat_rate_control,
        structural_control=structural_control,
    )
    planet = p.args.environment_model.planet
    final_metrics = _edg_orbit_metrics_from_rv(track.positions[end], track.velocities[end], mass, planet)
    return (
        switch_time_s=switch_time_s,
        duration_s=last(track.time),
        energy_jkg=final_metrics.energy,
        periapsis_radius_m=final_metrics.periapsis,
        apoapsis_radius_m=final_metrics.apoapsis,
        track=track,
        alpha_profile=track.alpha_profile,
    )
end

function _edg_predict_max_energy_depletion_outcome(
    config::AerobrakingEnergyDepletionConfig,
    p::ODEParams,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass_state::Float64,
    t::Float64,
    env,
    heat_load_j_cm2::Float64;
    heat_rate_control::Bool,
    structural_control::Bool,
)
    mass = _edg_predict_mass(spacecraft, mass_state)
    duration = _edg_drag_passage_duration(config, p, pos, vel, mass)
    times = _edg_targeting_prediction_time_grid(duration)
    heat_load_switches = (:heat_load in config.max_energy_submodes) ?
        _edg_solve_heat_load_switches(
            config,
            p,
            spacecraft,
            pos,
            vel,
            mass,
            env,
            heat_load_j_cm2,
            t;
            heat_rate_control=heat_rate_control,
            structural_control=structural_control,
        ) :
        (Inf, Inf)
    track = _edg_integrated_max_energy_depletion_trajectory(
        config,
        p,
        spacecraft,
        pos,
        vel,
        mass,
        t,
        times,
        heat_load_switches;
        heat_rate_control=heat_rate_control,
        structural_control=structural_control,
    )
    planet = p.args.environment_model.planet
    final_metrics = _edg_orbit_metrics_from_rv(track.positions[end], track.velocities[end], mass, planet)
    return (
        switch_time_s=NaN,
        heat_load_switches_s=heat_load_switches,
        duration_s=last(track.time),
        energy_jkg=final_metrics.energy,
        periapsis_radius_m=final_metrics.periapsis,
        apoapsis_radius_m=final_metrics.apoapsis,
        track=track,
        alpha_profile=track.alpha_profile,
    )
end

function _edg_targeting_switch_outcomes(
    config::AerobrakingEnergyDepletionConfig,
    p::ODEParams,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass::Float64,
    t::Float64;
    heat_rate_control::Bool,
    structural_control::Bool,
)
    low_drag = _edg_predict_targeting_outcome(
        config,
        p,
        spacecraft,
        pos,
        vel,
        mass,
        t,
        t;
        heat_rate_control=heat_rate_control,
        structural_control=structural_control,
    )
    high_drag = _edg_predict_targeting_outcome(
        config,
        p,
        spacecraft,
        pos,
        vel,
        mass,
        t,
        t + low_drag.duration_s + 1.0;
        heat_rate_control=heat_rate_control,
        structural_control=structural_control,
    )
    return low_drag, high_drag
end

function _edg_targeting_bracket_outcomes(
    config::AerobrakingEnergyDepletionConfig,
    p::ODEParams,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass::Float64,
    t::Float64;
    heat_load_j_cm2::Float64=0.0,
    heat_rate_control::Bool,
    structural_control::Bool,
)
    low_drag = _edg_predict_targeting_outcome(
        config,
        p,
        spacecraft,
        pos,
        vel,
        mass,
        t,
        t;
        heat_rate_control=heat_rate_control,
        structural_control=structural_control,
    )
    env = _edg_targeting_prediction_environment(p, pos, vel, t)
    max_energy_depletion = _edg_predict_max_energy_depletion_outcome(
        config,
        p,
        spacecraft,
        pos,
        vel,
        mass,
        t,
        env,
        heat_load_j_cm2;
        heat_rate_control=heat_rate_control,
        structural_control=structural_control,
    )
    return low_drag, max_energy_depletion
end

function _edg_targeting_outcome_with_heat_load(
    config::AerobrakingEnergyDepletionConfig,
    p::ODEParams,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass::Float64,
    t::Float64,
    switch_time_s::Float64,
    accumulated_heat_load_j_cm2::Float64;
    heat_rate_control::Bool,
    structural_control::Bool,
)
    outcome = _edg_predict_targeting_outcome(
        config,
        p,
        spacecraft,
        pos,
        vel,
        mass,
        t,
        switch_time_s;
        heat_rate_control=heat_rate_control,
        structural_control=structural_control,
    )
    future_heat_load = _edg_profile_heat_load(
        config,
        p,
        outcome.track,
        outcome.alpha_profile;
        heat_rate_control=false,
    )
    return merge(outcome, (heat_load_j_cm2=accumulated_heat_load_j_cm2 + future_heat_load,))
end

function _edg_certify_targeting_candidates(
    candidate_times::AbstractVector{<:Real},
    evaluate_candidate;
    heat_load_limit_j_cm2::Real,
    energy_order_tolerance_jkg::Real,
    heat_load_tolerance_j_cm2::Real,
)
    isempty(candidate_times) && throw(ArgumentError("candidate_times must not be empty."))
    energy_tolerance = Float64(energy_order_tolerance_jkg)
    heat_tolerance = Float64(heat_load_tolerance_j_cm2)
    heat_limit = Float64(heat_load_limit_j_cm2)
    energy_tolerance >= 0.0 || throw(ArgumentError("energy_order_tolerance_jkg must be >= 0.0."))
    heat_tolerance >= 0.0 || throw(ArgumentError("heat_load_tolerance_j_cm2 must be >= 0.0."))

    first_time = Float64(first(candidate_times))
    first_outcome = evaluate_candidate(first_time)
    certified_times = Float64[]
    certified_outcomes = typeof(first_outcome)[]

    if !isfinite(first_outcome.energy_jkg)
        return (times=certified_times, outcomes=certified_outcomes, failure=:nonfinite_energy, failure_time_s=first_time)
    end
    if !(isfinite(first_outcome.heat_load_j_cm2) && first_outcome.heat_load_j_cm2 <= heat_limit + heat_tolerance)
        return (times=certified_times, outcomes=certified_outcomes, failure=:heat_load, failure_time_s=first_time)
    end

    push!(certified_times, first_time)
    push!(certified_outcomes, first_outcome)
    previous_energy = first_outcome.energy_jkg

    for candidate_time in Iterators.drop(candidate_times, 1)
        time_s = Float64(candidate_time)
        outcome = evaluate_candidate(time_s)
        if !isfinite(outcome.energy_jkg)
            return (times=certified_times, outcomes=certified_outcomes, failure=:nonfinite_energy, failure_time_s=time_s)
        end
        if !(isfinite(outcome.heat_load_j_cm2) && outcome.heat_load_j_cm2 <= heat_limit + heat_tolerance)
            return (times=certified_times, outcomes=certified_outcomes, failure=:heat_load, failure_time_s=time_s)
        end
        if !(outcome.energy_jkg < previous_energy - energy_tolerance)
            return (times=certified_times, outcomes=certified_outcomes, failure=:energy_order, failure_time_s=time_s)
        end
        push!(certified_times, time_s)
        push!(certified_outcomes, outcome)
        previous_energy = outcome.energy_jkg
    end

    return (times=certified_times, outcomes=certified_outcomes, failure=:none, failure_time_s=Inf)
end

function _edg_disable_uncertified_targeting!(
    config::AerobrakingEnergyDepletionConfig,
    state::AerobrakingEnergyDepletionState,
    i::Int;
    prefer_max_energy_depletion::Bool,
)
    state.targeting_active[i] = false
    if prefer_max_energy_depletion && (:max_energy_depletion in config.guidance_modes)
        state.selected_mode[i] = :max_energy_depletion
        state.safe_low_drag[i] = false
    else
        state.selected_mode[i] = :safe_low_drag
        state.safe_low_drag[i] = true
    end
    return Inf
end

function _edg_solve_targeting_switch(
    config::AerobrakingEnergyDepletionConfig,
    state::AerobrakingEnergyDepletionState,
    p::ODEParams,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass::Float64,
    t::Float64,
    i::Int;
    heat_load_j_cm2::Float64=0.0,
    heat_rate_control::Bool,
    structural_control::Bool,
)
    isfinite(state.target_energy_jkg[i]) || return t + 0.5 * config.planning_horizon_s

    predicted_mass = _edg_predict_mass(spacecraft, mass)
    duration = _edg_drag_passage_duration(config, p, pos, vel, predicted_mass)
    t_low = t
    nominal_t_high = t + duration + 1.0
    candidate_times = collect(range(t_low, nominal_t_high; length=config.targeting_certification_samples))
    heat_load_limit = (:heat_load in config.max_energy_submodes) ? config.heat_load_limit_j_cm2 : Inf

    evaluate_candidate(t_switch) = _edg_targeting_outcome_with_heat_load(
        config,
        p,
        spacecraft,
        pos,
        vel,
        mass,
        t,
        t_switch,
        heat_load_j_cm2;
        heat_rate_control=heat_rate_control,
        structural_control=structural_control,
    )
    certification = _edg_certify_targeting_candidates(
        candidate_times,
        evaluate_candidate;
        heat_load_limit_j_cm2=heat_load_limit,
        energy_order_tolerance_jkg=config.targeting_energy_order_tolerance_jkg,
        heat_load_tolerance_j_cm2=config.targeting_heat_load_tolerance_j_cm2,
    )
    if length(certification.times) < 2
        return _edg_disable_uncertified_targeting!(
            config,
            state,
            i;
            prefer_max_energy_depletion=false,
        )
    end

    low_drag = first(certification.outcomes)
    certified_end = last(certification.outcomes)
    t_high = last(certification.times)
    target_energy = state.target_energy_jkg[i]
    target_apoapsis = config.target_apoapsis_radius_m
    state.bracket_min_energy_jkg[i] = certified_end.energy_jkg
    state.bracket_max_energy_jkg[i] = low_drag.energy_jkg

    energy_tolerance = config.targeting_energy_order_tolerance_jkg
    if target_energy > low_drag.energy_jkg + energy_tolerance
        return _edg_disable_uncertified_targeting!(
            config,
            state,
            i;
            prefer_max_energy_depletion=false,
        )
    elseif target_energy < certified_end.energy_jkg - energy_tolerance
        return _edg_disable_uncertified_targeting!(
            config,
            state,
            i;
            prefer_max_energy_depletion=true,
        )
    end

    function energy_residual(t_switch)
        outcome = evaluate_candidate(t_switch)
        return outcome.energy_jkg - target_energy
    end

    function apoapsis_residual(t_switch)
        outcome = evaluate_candidate(t_switch)
        return outcome.apoapsis_radius_m - target_apoapsis
    end

    function solve_energy_switch()
        f_low = low_drag.energy_jkg - target_energy
        f_high = certified_end.energy_jkg - target_energy
        if !(isfinite(f_low) && isfinite(f_high)) || f_low * f_high > 0.0
            denom = certified_end.energy_jkg - low_drag.energy_jkg
            frac = abs(denom) > eps(Float64) ? (target_energy - low_drag.energy_jkg) / denom : 0.5
            return t_low + clamp(frac, 0.0, 1.0) * (t_high - t_low)
        end
        return Roots.find_zero(energy_residual, (t_low, t_high), Roots.Brent(); rtol=1e-7)
    end

    function solve_apoapsis_switch()
        isfinite(target_apoapsis) && target_apoapsis > 0.0 || return solve_energy_switch()
        f_low = low_drag.apoapsis_radius_m - target_apoapsis
        f_high = certified_end.apoapsis_radius_m - target_apoapsis
        if !(isfinite(f_low) && isfinite(f_high)) || f_low * f_high > 0.0
            denom = certified_end.apoapsis_radius_m - low_drag.apoapsis_radius_m
            if isfinite(denom) && abs(denom) > eps(Float64)
                frac = (target_apoapsis - low_drag.apoapsis_radius_m) / denom
                return t_low + clamp(frac, 0.0, 1.0) * (t_high - t_low)
            end
            return solve_energy_switch()
        end
        return Roots.find_zero(apoapsis_residual, (t_low, t_high), Roots.Brent(); rtol=1e-7)
    end

    t_switch = solve_apoapsis_switch()
    for _ in 1:2
        outcome = evaluate_candidate(t_switch)
        apo_error = outcome.apoapsis_radius_m - config.target_apoapsis_radius_m
        state.target_energy_jkg[i] = outcome.energy_jkg
        if isfinite(apo_error) && abs(apo_error) <= 25.0
            break
        end
        denom = outcome.apoapsis_radius_m + outcome.periapsis_radius_m
        if !(isfinite(apo_error) && isfinite(denom) && denom > 0.0)
            break
        end
        energy_correction = p.args.environment_model.planet.μ / denom^2 * apo_error
        isfinite(energy_correction) || break
        target_energy -= energy_correction
        state.target_energy_jkg[i] = target_energy
        t_switch = solve_apoapsis_switch()
    end

    return clamp(t_switch, t_low, t_high)
end
