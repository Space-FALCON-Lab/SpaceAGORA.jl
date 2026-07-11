using Roots

@inline function _edg_heat_load_scale_height(p::ODEParams)
    density_model = p.args.environment_model.density_model
    planet = p.args.environment_model.planet
    if hasproperty(density_model, :H)
        return max(Float64(density_model.H), 1.0)
    end
    return hasproperty(planet, :H) ? max(Float64(planet.H), 1.0) : 7_000.0
end

@inline function _edg_total_ref_area(spacecraft)::Float64
    area = 0.0
    for link in spacecraft.links
        link_area = Float64(link.ref_area)
        isfinite(link_area) && link_area > 0.0 && (area += link_area)
    end
    return max(area, eps(Float64))
end

@inline function _edg_predict_mass(spacecraft, mass::Float64)::Float64
    if isfinite(mass) && mass > 0.0
        return mass
    end
    total = 0.0
    for link in spacecraft.links
        total += max(0.0, Float64(link.m))
    end
    if hasproperty(spacecraft, :prop_mass)
        total += max(0.0, Float64(spacecraft.prop_mass))
    end
    return max(total, 1.0)
end

@inline function _edg_max_heat_load_for_links(sc, links::Tuple{Vararg{Int}})::Float64
    hasproperty(sc, :heat_loads) || return 0.0
    heat_loads = sc.heat_loads
    value = 0.0
    for idx in links
        if 1 <= idx <= length(heat_loads)
            candidate = Float64(heat_loads[idx])
            if isfinite(candidate)
                value = max(value, max(0.0, candidate))
            end
        end
    end
    return value
end

function _edg_weighted_aero_coefficients(spacecraft, temperature::Float64, speed_ratio::Float64, alpha::Float64)
    area = _edg_total_ref_area(spacecraft)
    lift_area = 0.0
    drag_area = 0.0
    aero_module = getfield(getfield(parentmodule(@__MODULE__), :DynamicEffectors), :AerodynamicEffectors)
    aero_coeff = getfield(aero_module, :aerodynamic_coefficient_fM)
    for link in spacecraft.links
        link_area = max(0.0, Float64(link.ref_area))
        link_area == 0.0 && continue
        coeffs = aero_coeff(link, temperature, speed_ratio, alpha, Float64(link.β), Float64(link.θ))
        lift_area += Float64(coeffs[1]) * link_area
        drag_area += max(0.0, Float64(coeffs[2])) * link_area
    end
    return lift_area / area, max(drag_area / area, eps(Float64))
end

function _edg_heat_load_coefficients(config, p::ODEParams, spacecraft, env)
    temperature = env.temperature
    speed_ratio = max(env.molecular_speed_ratio, eps(Float64))
    cl_low, cd_low = _edg_weighted_aero_coefficients(spacecraft, temperature, speed_ratio, config.min_alpha_rad)
    _, cd_high = _edg_weighted_aero_coefficients(spacecraft, temperature, speed_ratio, config.max_alpha_rad)
    cd_slope = (cd_high - cd_low) / max(config.max_alpha_rad - config.min_alpha_rad, eps(Float64))
    if !(isfinite(cd_slope) && abs(cd_slope) > eps(Float64))
        cd_slope = (2.2 - 0.8) / (pi / 2)
        cd_low = 0.8
    end
    return (cd_slope=cd_slope, cl_low=cl_low, cd_low=cd_low)
end

@inline function _edg_eccentric_anomaly_from_true(nu::Float64, e::Float64)::Float64
    E = atan(sqrt(max(0.0, 1.0 - e^2)) * sin(nu), e + cos(nu))
    return mod(E, 2pi)
end

@inline function _edg_mean_anomaly_from_true(nu::Float64, e::Float64)::Float64
    E = _edg_eccentric_anomaly_from_true(nu, e)
    return mod(E - e * sin(E), 2pi)
end

function _edg_drag_passage_duration(config, p::ODEParams, pos::SVector{3, Float64}, vel::SVector{3, Float64}, mass::Float64)
    planet = p.args.environment_model.planet
    oe = rvtoorbitalelement(pos, vel, mass, planet)
    a, e, nu = oe[1], oe[2], oe[6]
    if !(isfinite(a) && isfinite(e) && 0.0 <= e < 1.0 && a > 0.0)
        return config.planning_horizon_s
    end
    nu_exit = nu > pi ? 2pi - nu : nu
    n = sqrt(planet.μ / a^3)
    M0 = _edg_mean_anomaly_from_true(nu, e)
    M_exit = _edg_mean_anomaly_from_true(nu_exit, e)
    duration = mod(M_exit - M0, 2pi) / n
    if !(isfinite(duration) && duration > 1.0)
        return min(config.planning_horizon_s, 1_000.0)
    end
    return min(duration, config.planning_horizon_s)
end

function _edg_prediction_time_grid(duration::Float64)
    step = 1.0
    n = clamp(ceil(Int, duration / step) + 1, 64, 2_000)
    return collect(range(0.0, duration; length=n))
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

function _edg_closed_form_heat_load_trajectory(
    config,
    p::ODEParams,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass::Float64,
    t::Float64,
    env,
)
    planet = p.args.environment_model.planet
    duration = _edg_drag_passage_duration(config, p, pos, vel, mass)
    times = _edg_prediction_time_grid(duration)
    n = length(times)
    h = Vector{Float64}(undef, n)
    gamma = Vector{Float64}(undef, n)
    speed = Vector{Float64}(undef, n)
    rho = Vector{Float64}(undef, n)
    temperature = Vector{Float64}(undef, n)
    speed_ratio = Vector{Float64}(undef, n)

    r0 = norm(pos)
    v0 = max(norm(vel), eps(Float64))
    gamma0 = asin(clamp(dot(pos, vel) / max(r0 * v0, eps(Float64)), -1.0, 1.0))
    altitude0 = env.altitude_m
    oe = rvtoorbitalelement(pos, vel, mass, planet)
    a = isfinite(oe[1]) && oe[1] > 0.0 ? oe[1] : r0
    t_peri = max(0.5 * duration, eps(Float64))

    for j in eachindex(times)
        tau = times[j]
        altitude = altitude0 + v0 * gamma0 * (tau - tau^2 / (2.0 * t_peri))
        altitude = max(0.0, altitude)
        radius = planet.Rp_e + altitude
        v = sqrt(max(0.0, planet.μ * (2.0 / radius - 1.0 / a)))
        hdot = v0 * gamma0 * (1.0 - tau / t_peri)
        g = asin(clamp(hdot / max(v, eps(Float64)), -1.0, 1.0))
        r, T = _edg_sample_prediction_atmosphere(p, altitude, t + tau)
        sound_speed = sqrt(max(0.0, planet.γ * planet.R * T))
        h[j] = altitude
        gamma[j] = g
        speed[j] = max(v, eps(Float64))
        rho[j] = r
        temperature[j] = T
        speed_ratio[j] = sound_speed > 0.0 ? sqrt(0.5 * planet.γ) * speed[j] / sound_speed : 0.0
    end

    return (time=times, h=h, gamma=gamma, speed=speed, rho=rho, temperature=temperature, speed_ratio=speed_ratio)
end

function _edg_integrated_heat_load_trajectory(
    config,
    p::ODEParams,
    spacecraft,
    pos0::SVector{3, Float64},
    vel0::SVector{3, Float64},
    mass::Float64,
    t::Float64,
    times::Vector{Float64},
    alpha_profile::Vector{Float64},
    altitude_offset_m::Float64=0.0,
)
    planet = p.args.environment_model.planet
    area = _edg_total_ref_area(spacecraft)
    n = length(times)
    positions = Vector{SVector{3, Float64}}(undef, n)
    velocities = Vector{SVector{3, Float64}}(undef, n)
    positions[1] = pos0
    velocities[1] = vel0

    function acceleration(r, v, tau, alpha)
        radius = max(norm(r), eps(Float64))
        gravity = -planet.μ * r / radius^3
        altitude = radius - planet.Rp_e + altitude_offset_m
        rho, T = _edg_sample_prediction_atmosphere(p, altitude, t + tau)
        speed = norm(v)
        if !(rho > 0.0 && speed > 0.0)
            return gravity
        end
        sound_speed = sqrt(max(0.0, planet.γ * planet.R * T))
        S = sound_speed > 0.0 ? sqrt(0.5 * planet.γ) * speed / sound_speed : 0.0
        _, cd = _edg_weighted_aero_coefficients(spacecraft, T, S, alpha)
        drag = -0.5 * rho * speed^2 * cd * area / mass * (v / speed)
        return gravity + drag
    end

    for j in 1:(n - 1)
        dt = times[j + 1] - times[j]
        alpha = alpha_profile[min(j, length(alpha_profile))]
        r = positions[j]
        v = velocities[j]
        a1 = acceleration(r, v, times[j], alpha)
        k1r, k1v = v, a1
        a2 = acceleration(r + 0.5 * dt * k1r, v + 0.5 * dt * k1v, times[j] + 0.5 * dt, alpha)
        k2r, k2v = v + 0.5 * dt * k1v, a2
        a3 = acceleration(r + 0.5 * dt * k2r, v + 0.5 * dt * k2v, times[j] + 0.5 * dt, alpha)
        k3r, k3v = v + 0.5 * dt * k2v, a3
        a4 = acceleration(r + dt * k3r, v + dt * k3v, times[j + 1], alpha)
        k4r, k4v = v + dt * k3v, a4
        positions[j + 1] = r + dt * (k1r + 2.0 * k2r + 2.0 * k3r + k4r) / 6.0
        velocities[j + 1] = v + dt * (k1v + 2.0 * k2v + 2.0 * k3v + k4v) / 6.0
    end

    h = Vector{Float64}(undef, n)
    gamma = Vector{Float64}(undef, n)
    speed = Vector{Float64}(undef, n)
    rho = Vector{Float64}(undef, n)
    temperature = Vector{Float64}(undef, n)
    speed_ratio = Vector{Float64}(undef, n)
    for j in eachindex(times)
        r = positions[j]
        v = velocities[j]
        radius = max(norm(r), eps(Float64))
        vel_norm = max(norm(v), eps(Float64))
        altitude = radius - planet.Rp_e + altitude_offset_m
        r_atm, T = _edg_sample_prediction_atmosphere(p, altitude, t + times[j])
        sound_speed = sqrt(max(0.0, planet.γ * planet.R * T))
        h[j] = altitude
        gamma[j] = asin(clamp(dot(r, v) / max(radius * vel_norm, eps(Float64)), -1.0, 1.0))
        speed[j] = vel_norm
        rho[j] = r_atm
        temperature[j] = T
        speed_ratio[j] = sound_speed > 0.0 ? sqrt(0.5 * planet.γ) * vel_norm / sound_speed : 0.0
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
    )
end

function _edg_heat_load_lambdas(track, alpha_profile::Vector{Float64}, coeffs, mass::Float64, area::Float64, planet, scale_height::Float64, k::Float64)
    n = length(track.time)
    lambda_v = zeros(Float64, n)
    lambda_gamma = zeros(Float64, n)
    lambda_h = zeros(Float64, n)
    lambda_switch = Vector{Float64}(undef, n)

    lambda_v[end] = track.speed[end]
    lambda_h[end] = planet.μ / max((planet.Rp_e + track.h[end])^2, eps(Float64))

    for j in n:-1:2
        jm = j - 1
        dt = track.time[j] - track.time[jm]
        v = max(track.speed[jm], eps(Float64))
        h = track.h[jm]
        gamma = track.gamma[jm]
        rho = max(track.rho[jm], 0.0)
        alpha = alpha_profile[jm]
        radius = max(planet.Rp_e + h, eps(Float64))
        g = planet.g_ref * planet.Rp_e^2 / radius^2
        cd = coeffs.cd_low + alpha * coeffs.cd_slope

        lambda_v_dot =
            -3.0 * k * rho * v^2 * alpha / pi +
            lambda_v[j] * (rho * area * cd * v) / mass -
            lambda_gamma[j] * ((rho * area * coeffs.cl_low) / (2.0 * mass) + g / v^2 + 1.0 / radius) -
            lambda_h[j] * gamma

        lambda_gamma_dot = lambda_v[j] * g - lambda_h[j] * v

        lambda_h_dot =
            k * rho * v^3 * alpha / (pi * scale_height) -
            lambda_v[j] * ((rho * area * cd * v^2) / (2.0 * mass * scale_height) + 2.0 * g * gamma / radius) +
            lambda_gamma[j] * (rho * area * coeffs.cl_low * v / (2.0 * mass * scale_height) - 2.0 * g / (radius * v) + v / radius^2)

        lambda_v[jm] = lambda_v[j] - lambda_v_dot * dt
        lambda_gamma[jm] = lambda_gamma[j] - lambda_gamma_dot * dt
        lambda_h[jm] = lambda_h[j] - lambda_h_dot * dt
    end

    denom = max(abs(area * coeffs.cd_slope * pi), eps(Float64))
    for j in eachindex(lambda_switch)
        lambda_switch[j] = k * 2.0 * mass * track.speed[j] / denom
    end
    return lambda_switch, lambda_v
end

function _edg_heat_load_alpha_profile(config, track, coeffs, mass::Float64, area::Float64, planet, scale_height::Float64, k::Float64, seed_profile::Vector{Float64})
    lambda_switch, lambda_v = _edg_heat_load_lambdas(track, seed_profile, coeffs, mass, area, planet, scale_height, k)
    profile = similar(seed_profile)
    for j in eachindex(profile)
        profile[j] = lambda_v[j] >= lambda_switch[j] ? config.max_alpha_rad : config.min_alpha_rad
    end
    return profile, lambda_switch, lambda_v
end

function _edg_heat_load_track_env(track, j::Int)
    speed = max(track.speed[j], 0.0)
    rho = max(track.rho[j], 0.0)
    return (
        altitude_m=track.h[j],
        rho=rho,
        temperature=track.temperature[j],
        speed=speed,
        molecular_speed_ratio=track.speed_ratio[j],
        dynamic_pressure=0.5 * rho * speed^2,
    )
end

function _edg_constrained_heat_load_alpha_profile(
    config,
    p::ODEParams,
    spacecraft,
    controlled_panel_links::Tuple{Vararg{Int}},
    track,
    raw_profile::Vector{Float64};
    heat_rate_control::Bool,
    structural_control::Bool,
)
    profile = similar(raw_profile)
    alpha_past = config.max_alpha_rad
    thermal_model = p.args.environment_model.thermal_model
    taf = hasproperty(thermal_model, :thermal_accomodation_factor) ? Float64(thermal_model.thermal_accomodation_factor) : 1.0
    planet = p.args.environment_model.planet

    for j in eachindex(raw_profile)
        base_alpha = clamp(raw_profile[j], config.min_alpha_rad, config.max_alpha_rad)
        if base_alpha <= config.min_alpha_rad + 1e-8
            profile[j] = config.min_alpha_rad
            alpha_past = profile[j]
            continue
        end

        env = _edg_heat_load_track_env(track, j)
        alpha_hr = base_alpha
        alpha_struct = base_alpha

        if heat_rate_control
            alpha_hr = _energy_depletion_heatrate_root_alpha(
                taf=taf,
                rho=env.rho,
                T_p=env.temperature,
                R=planet.R,
                gamma=planet.γ,
                S=env.molecular_speed_ratio,
                max_alpha=base_alpha,
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
                controlled_panel_links,
                base_alpha,
            )
        end

        profile[j] = clamp(min(alpha_hr, alpha_struct), config.min_alpha_rad, config.max_alpha_rad)
        alpha_past = profile[j]
    end
    return profile
end

function _edg_profile_heat_rates(config, p::ODEParams, track, alpha_profile::Vector{Float64}; heat_rate_control::Bool)
    thermal_model = p.args.environment_model.thermal_model
    taf = hasproperty(thermal_model, :thermal_accomodation_factor) ? Float64(thermal_model.thermal_accomodation_factor) : 1.0
    planet = p.args.environment_model.planet
    heat_rate_limit = config.heat_rate_limit_w_cm2
    qdot = Vector{Float64}(undef, length(track.time))
    for j in eachindex(track.time)
        rate = _energy_depletion_heat_rate_calc(
            taf,
            max(0.0, track.rho[j]),
            track.temperature[j],
            track.temperature[j],
            planet.R,
            planet.γ,
            max(0.0, track.speed_ratio[j]),
            alpha_profile[j],
        )
        if heat_rate_control && alpha_profile[j] > config.min_alpha_rad + 1e-8 &&
                isfinite(heat_rate_limit) && heat_rate_limit > 0.0
            rate = min(rate, heat_rate_limit)
        end
        qdot[j] = isfinite(rate) && rate > 0.0 ? rate : 0.0
    end
    return qdot
end

function _edg_integrate_series(time::Vector{Float64}, values::Vector{Float64})
    total = 0.0
    for j in 1:(length(time) - 1)
        total += 0.5 * (values[j] + values[j + 1]) * (time[j + 1] - time[j])
    end
    return total
end

function _edg_profile_heat_load(config, p::ODEParams, track, alpha_profile::Vector{Float64}; heat_rate_control::Bool)
    qdot = _edg_profile_heat_rates(config, p, track, alpha_profile; heat_rate_control=heat_rate_control)
    return _edg_integrate_series(track.time, qdot)
end

function _edg_first_low_alpha_interval_indices(config, alpha_profile::Vector{Float64})
    low = map(alpha -> alpha <= config.min_alpha_rad + 1e-8, alpha_profile)
    first_low = findfirst(identity, low)
    first_low === nothing && return nothing

    last_low = first_low
    while last_low < lastindex(low) && low[last_low + 1]
        last_low += 1
    end
    return (first_low, last_low)
end

function _edg_first_two_switch_alpha_profile(config, alpha_profile::Vector{Float64}, high_profile::Vector{Float64}=fill(config.max_alpha_rad, length(alpha_profile)))
    applied_profile = copy(high_profile)
    interval = _edg_first_low_alpha_interval_indices(config, alpha_profile)
    interval === nothing && return applied_profile

    first_low, last_low = interval
    applied_profile[first_low:last_low] .= config.min_alpha_rad
    return applied_profile
end

function _edg_low_alpha_switch_window(config, t::Float64, track, alpha_profile::Vector{Float64})
    interval = _edg_first_low_alpha_interval_indices(config, alpha_profile)
    interval === nothing && return (Inf, Inf)

    first_low, last_low = interval
    if first_low == last_low
        return (Inf, Inf)
    end
    return (t + track.time[first_low], t + track.time[last_low])
end

function _edg_balanced_tpbvp_heat_load_window(
    config,
    p::ODEParams,
    t::Float64,
    track,
    spacecraft,
    controlled_panel_links::Tuple{Vararg{Int}},
    heat_load_j_cm2::Float64,
    target_load::Float64;
    heat_rate_control::Bool,
    structural_control::Bool,
)
    n = length(track.time)
    n < 3 && return (Inf, Inf)

    high_profile = _edg_constrained_heat_load_alpha_profile(
        config,
        p,
        spacecraft,
        controlled_panel_links,
        track,
        fill(config.max_alpha_rad, n);
        heat_rate_control=heat_rate_control,
        structural_control=structural_control,
    )
    low_profile = fill(config.min_alpha_rad, n)
    qdot_high = _edg_profile_heat_rates(config, p, track, high_profile; heat_rate_control=false)
    qdot_low = _edg_profile_heat_rates(config, p, track, low_profile; heat_rate_control=false)
    total_high = heat_load_j_cm2 + _edg_integrate_series(track.time, qdot_high)
    required_saving = total_high - target_load
    required_saving <= 0.0 && return (Inf, Inf)

    savings = max.(qdot_high .- qdot_low, 0.0)
    maximum(savings) <= 0.0 && return (Inf, Inf)

    cumulative = zeros(Float64, n)
    for j in 1:(n - 1)
        cumulative[j + 1] = cumulative[j] + 0.5 * (savings[j] + savings[j + 1]) * (track.time[j + 1] - track.time[j])
    end
    if cumulative[end] < required_saving
        return (t + first(track.time), t + last(track.time))
    end

    best_start = 0
    best_stop = 0
    best_over = Inf
    over_tolerance = 0.02
    for start_idx in 1:(n - 1)
        target_cumulative = cumulative[start_idx] + required_saving
        stop_idx = searchsortedfirst(cumulative, target_cumulative)
        if start_idx < stop_idx <= n
            over = cumulative[stop_idx] - cumulative[start_idx] - required_saving
            if over < best_over - over_tolerance || (over <= best_over + over_tolerance && stop_idx > best_stop)
                best_start = start_idx
                best_stop = stop_idx
                best_over = over
            end
        end
    end

    if best_start == 0 || best_stop == 0 || best_start == best_stop
        return (Inf, Inf)
    end
    return (t + track.time[best_start], t + track.time[best_stop])
end

@inline function _edg_padded_heat_load_window(config, t::Float64, track, window::Tuple{Float64, Float64})
    all(isfinite, window) || return window
    if config.heat_load_switch_solver == :tpbvp_integration
        return (
            max(t + first(track.time), window[1]),
            min(t + last(track.time), window[2]),
        )
    end

    duration = last(track.time) - first(track.time)
    return (
        max(t + first(track.time), window[1] - 10.0),
        min(t + last(track.time), window[2] + max(10.0, 0.15 * duration)),
    )
end

function _edg_heat_load_profile_for_k(
    config,
    p::ODEParams,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass::Float64,
    t::Float64,
    env,
    coeffs,
    k::Float64,
    controlled_panel_links::Tuple{Vararg{Int}},
    heat_rate_control::Bool,
    structural_control::Bool,
)
    planet = p.args.environment_model.planet
    area = _edg_total_ref_area(spacecraft)
    scale_height = _edg_heat_load_scale_height(p)
    base_track = _edg_closed_form_heat_load_trajectory(config, p, pos, vel, mass, t, env)
    alpha_profile = fill(config.max_alpha_rad, length(base_track.time))

    if config.heat_load_switch_solver == :tpbvp_integration
        altitude_offset_m = env.altitude_m - (norm(pos) - planet.Rp_e)
        alpha_profile = _edg_constrained_heat_load_alpha_profile(
            config,
            p,
            spacecraft,
            controlled_panel_links,
            base_track,
            alpha_profile;
            heat_rate_control=heat_rate_control,
            structural_control=structural_control,
        )
        track = _edg_integrated_heat_load_trajectory(
            config, p, spacecraft, pos, vel, mass, t, base_track.time, alpha_profile, altitude_offset_m
        )
        for _ in 1:2
            raw_profile, _, _ = _edg_heat_load_alpha_profile(config, track, coeffs, mass, area, planet, scale_height, k, alpha_profile)
            alpha_profile = _edg_constrained_heat_load_alpha_profile(
                config,
                p,
                spacecraft,
                controlled_panel_links,
                track,
                raw_profile;
                heat_rate_control=heat_rate_control,
                structural_control=structural_control,
            )
            track = _edg_integrated_heat_load_trajectory(
                config, p, spacecraft, pos, vel, mass, t, base_track.time, alpha_profile, altitude_offset_m
            )
        end
        raw_profile, _, _ = _edg_heat_load_alpha_profile(config, track, coeffs, mass, area, planet, scale_height, k, alpha_profile)
        alpha_profile = _edg_constrained_heat_load_alpha_profile(
            config,
            p,
            spacecraft,
            controlled_panel_links,
            track,
            raw_profile;
            heat_rate_control=heat_rate_control,
            structural_control=structural_control,
        )
        return track, alpha_profile
    end

    track = base_track
    for _ in 1:3
        alpha_profile, _, _ = _edg_heat_load_alpha_profile(config, track, coeffs, mass, area, planet, scale_height, k, alpha_profile)
    end
    return track, alpha_profile
end

function _edg_solve_heat_load_switches(
    config,
    p::ODEParams,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass_state::Float64,
    env,
    heat_load_j_cm2::Float64,
    t::Float64;
    heat_rate_control::Bool,
    structural_control::Bool,
)
    if !(isfinite(config.heat_load_limit_j_cm2) && config.heat_load_limit_j_cm2 > 0.0)
        return (Inf, Inf)
    end

    mass = _edg_predict_mass(spacecraft, mass_state)
    coeffs = _edg_heat_load_coefficients(config, p, spacecraft, env)
    target_margin = config.heat_load_switch_solver == :tpbvp_integration ?
        max(2.0, 6.5e-2 * config.heat_load_limit_j_cm2) :
        max(0.3, 1e-2 * config.heat_load_limit_j_cm2)
    target_load = config.heat_load_limit_j_cm2 - target_margin
    target_load = max(0.0, target_load)
    low_k = 0.0
    high_k = lowercase(String(p.args.environment_model.planet.name)) == "mars" ? 0.1 : 10.0

    last_track = Ref{Any}(nothing)
    last_profile = Ref{Vector{Float64}}(Float64[])
    best_under_track = Ref{Any}(nothing)
    best_under_profile = Ref{Vector{Float64}}(Float64[])
    best_under_residual = Ref(Inf)

    function residual(k)
        track, raw_profile = _edg_heat_load_profile_for_k(
            config,
            p,
            spacecraft,
            pos,
            vel,
            mass,
            t,
            env,
            coeffs,
            k,
            config.controlled_panel_links,
            heat_rate_control,
            structural_control,
        )
        high_profile = _edg_constrained_heat_load_alpha_profile(
            config,
            p,
            spacecraft,
            config.controlled_panel_links,
            track,
            fill(config.max_alpha_rad, length(raw_profile));
            heat_rate_control=heat_rate_control,
            structural_control=structural_control,
        )
        profile = _edg_first_two_switch_alpha_profile(config, raw_profile, high_profile)
        last_track[] = track
        last_profile[] = profile
        predicted_load = heat_load_j_cm2 + _edg_profile_heat_load(config, p, track, profile; heat_rate_control=false)
        f = predicted_load - target_load
        if f <= 0.0 && abs(f) < best_under_residual[]
            best_under_track[] = track
            best_under_profile[] = copy(profile)
            best_under_residual[] = abs(f)
        end
        return f
    end

    f_low = residual(low_k)
    low_track = last_track[]
    low_profile = copy(last_profile[])
    if f_low <= 0.0
        return (Inf, Inf)
    end

    f_high = residual(high_k)
    if f_high > 0.0
        return (t + first(last_track[].time), t + last(last_track[].time))
    end

    bracket = (low_k, high_k)
    previous_k = low_k
    previous_f = f_low
    for k in range(low_k, high_k; length=21)[2:end]
        f = residual(k)
        if previous_f * f <= 0.0
            bracket = (previous_k, k)
            break
        end
        previous_k = k
        previous_f = f
    end

    k_sol = try
        Roots.find_zero(residual, bracket, Roots.Brent())
    catch
        Roots.find_zero(residual, bracket, Roots.Bisection())
    end
    root_residual = residual(k_sol)
    track = last_track[]
    profile = copy(last_profile[])
    if root_residual > 0.0 && best_under_track[] !== nothing
        track = best_under_track[]
        profile = copy(best_under_profile[])
    end
    if config.heat_load_switch_solver == :tpbvp_integration
        balanced_window = _edg_balanced_tpbvp_heat_load_window(
            config,
            p,
            t,
            low_track,
            spacecraft,
            config.controlled_panel_links,
            heat_load_j_cm2,
            target_load;
            heat_rate_control=heat_rate_control,
            structural_control=structural_control,
        )
        if all(isfinite, balanced_window) && balanced_window[1] < balanced_window[2]
            return balanced_window
        end
    end
    window = _edg_low_alpha_switch_window(config, t, track, profile)
    if all(isfinite, window) && window[1] < window[2]
        return _edg_padded_heat_load_window(config, t, track, window)
    end

    fallback = _edg_low_alpha_switch_window(config, t, low_track, low_profile)
    if all(isfinite, fallback) && fallback[1] < fallback[2]
        return _edg_padded_heat_load_window(config, t, low_track, fallback)
    end
    duration = _edg_drag_passage_duration(config, p, pos, vel, mass)
    return (t + 0.25 * duration, t + 0.75 * duration)
end

@inline function _edg_heat_load_low_drag_active(model, t::Float64, i::Int)::Bool
    config = model.config
    state = model.state
    switches = state.heat_load_switches_s[i]
    return state.selected_mode[i] == :max_energy_depletion &&
        (:heat_load in config.max_energy_submodes) &&
        isfinite(switches[1]) &&
        isfinite(switches[2]) &&
        t >= switches[1] &&
        t <= switches[2]
end
