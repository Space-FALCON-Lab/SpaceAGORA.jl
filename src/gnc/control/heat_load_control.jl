using Roots
using NLsolve

const _EDG_HEAT_LOAD_SHOOTING_RESIDUAL_TOL = 1e-3

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

@inline function _edg_controlled_ref_area(spacecraft, controlled_panel_links::Tuple{Vararg{Int}})::Float64
    area = 0.0
    for idx in controlled_panel_links
        if 1 <= idx <= length(spacecraft.links)
            link_area = Float64(spacecraft.links[idx].ref_area)
            isfinite(link_area) && link_area > 0.0 && (area += link_area)
        end
    end
    return max(0.0, area)
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
    panel_links = Tuple(idx for (idx, link) in pairs(spacecraft.links) if !link.root)
    cl, cd = _edg_legacy_spacecraft_aero_coefficients(spacecraft, speed_ratio, panel_links, alpha)
    return cl, max(cd, eps(Float64))
end

function _edg_heat_load_coefficients(config, p::ODEParams, spacecraft, env)
    planet = p.args.environment_model.planet
    temperature = Float64(planet.T)
    speed_ratio = max(env.speed / sqrt(2.0 * planet.R * temperature), eps(Float64))
    cl_low, cd_low = _edg_weighted_aero_coefficients(spacecraft, temperature, speed_ratio, 0.0)
    _, cd_high = _edg_weighted_aero_coefficients(spacecraft, temperature, speed_ratio, config.max_alpha_rad)
    cd_slope = (cd_high - cd_low) / (pi / 2)
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
    if !(isfinite(a) && isfinite(e) && 0.0 < e < 1.0 && a > 0.0)
        return config.planning_horizon_s
    end
    nu_exit = if dot(pos, vel) < 0.0
        nu > pi ? 2pi - nu : nu
    else
        interface_radius = Float64(planet.Rp_e) + 1e3 * Float64(p.args.environment_model.EI)
        semi_latus_rectum = a * (1.0 - e^2)
        cos_nu_exit = (semi_latus_rectum / interface_radius - 1.0) / e
        isfinite(cos_nu_exit) || return config.planning_horizon_s
        acos(clamp(cos_nu_exit, -1.0, 1.0))
    end
    n = sqrt(planet.μ / a^3)
    M0 = _edg_mean_anomaly_from_true(nu, e)
    M_exit = _edg_mean_anomaly_from_true(nu_exit, e)
    duration = mod(M_exit - M0, 2pi) / n
    if !(isfinite(duration) && duration > 1.0)
        return min(config.planning_horizon_s, 1_000.0)
    end
    return min(duration, config.planning_horizon_s)
end

function _edg_prediction_time_grid(duration::Float64, step::Float64=0.1)
    n = clamp(ceil(Int, duration / step) + 1, 64, 20_000)
    return collect(range(0.0, duration; length=n))
end

@inline function _edg_chebyshev_t0_t4(x::Float64)
    return (1.0, x, 2.0x^2 - 1.0, 4.0x^3 - 3.0x, 8.0x^4 - 8.0x^2 + 1.0)
end

@inline function _edg_chebyshev4(x::Float64, coefficients::NTuple{5, Float64})
    basis = _edg_chebyshev_t0_t4(x)
    return sum(coefficients[j] * basis[j] for j in 1:5)
end

function _edg_chebyshev44(gamma_hat::Float64, velocity_hat::Float64, coefficients::NTuple{15, Float64})
    tg = _edg_chebyshev_t0_t4(gamma_hat)
    tv = _edg_chebyshev_t0_t4(velocity_hat)
    orders = ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2),
              (3, 0), (2, 1), (1, 2), (0, 3), (4, 0), (3, 1), (2, 2), (1, 3), (0, 4))
    return sum(coefficients[j] * tg[orders[j][1] + 1] * tv[orders[j][2] + 1] for j in eachindex(orders))
end

function _edg_closed_form_correction_fit(planet)
    name = lowercase(String(planet.name))
    if name == "mars"
        return (4350.0, 450.0, -5.25, 2.25,
            (4.560850, -1.715925, -0.375150, 0.376025, -0.254500),
            (-4.74628450, -3.69824500, -6.06348525, 1.07688200, 3.54740275,
             1.27053200, -0.26357900, -1.12508600, -0.91570600, -0.17288575,
             0.045299125, 0.280380250, 0.224278500, -0.004521000, -0.025270125))
    elseif name == "venus"
        return (9000.0, 600.0, -4.375, 1.375,
            (16.3939125, -0.6316500, -0.3870500, 0.3974500, -0.2209625),
            (-11.578333625, -6.131747500, -3.504510250, 1.440936250, 2.801654000,
             0.536784250, 0.049912500, -0.281389000, -0.520214000, -0.099717750,
             -0.094991375, -0.214584750, -0.037281250, 0.078927750, 0.016912750))
    elseif name == "earth"
        return (9075.0, 725.0, -5.125, 2.125,
            (19.807700, -3.164075, 0.592100, -0.113625, 0.043900),
            (-1.989471500, 0.523520500, -1.992902250, 0.462889750, 1.965380250,
             0.591760750, -0.182282500, -0.886610000, -0.871433000, -0.200818750,
             0.024188625, 0.286623750, 0.468237250, 0.332824000, 0.043238125))
    elseif name == "titan"
        return (2050.0, 150.0, -12.5, 6.5,
            (148.480400, 10.480925, 0.030500, -0.003725, 0.000600),
            (-8.397098625, -4.739801250, -1.473748250, 0.671454000, 1.024378750,
             0.184619500, -0.031104750, -0.118313000, -0.145654000, -0.027513750,
             -0.034697125, -0.048001000, 0.005169000, 0.019228250, 0.003789250))
    end
    throw(ArgumentError("No deprecated closed-form correction fit is defined for planet $(planet.name)."))
end

function _edg_closed_form_f1_f2(planet, v0::Float64, gamma0::Float64, times::Vector{Float64}, t_peri::Float64)
    velocity_center, velocity_half_width, gamma_center, gamma_half_width, f1_coefficients, f2_coefficients =
        _edg_closed_form_correction_fit(planet)
    velocity_hat = (v0 - velocity_center) / velocity_half_width
    gamma_hat = (rad2deg(gamma0) - gamma_center) / gamma_half_width
    f1 = _edg_chebyshev4(velocity_hat, f1_coefficients)
    log_f2 = _edg_chebyshev44(gamma_hat, velocity_hat, f2_coefficients)
    return f1, exp(log_f2) .* times ./ (2.0 * t_peri)
end

function _edg_closed_form_density_fit(planet)
    name = lowercase(String(planet.name))
    if name == "mars"
        return 175.0, 125.0, (-22.904887634115337, -11.841971600850178, 2.449273531108809,
            0.07122528331622979, -0.3143478095124491, 0.1734745257301972,
            -0.14396272246410283, 0.06104818577407287, 0.0738268254164405,
            -0.11983534975043726, 0.052506209635111954, 0.0022260853200571835,
            -0.0245022530854842)
    elseif name == "venus"
        return 150.0, 100.0, (-17.17603495187312, -15.971487996237625, 3.1916728583530163,
            0.8296122739608919, -0.891116626156514, -0.018054118127076735,
            0.2016004439406335, 0.0537350005662433, -0.06184791826147179,
            -0.002221184567720138, -0.03852447322489234, 0.11389766896205009,
            -0.023720859333601994)
    elseif name == "earth"
        return 275.5, 225.0, (-20.5897273251208688, -8.56737992269037107, 3.45019353462444567,
            -1.76485096058359936, 0.601801952591328626, 0.0724429484383713046,
            -0.330361927785729759, 0.307509524503914389, -0.158408358148278777,
            0.00777052279665210299, 0.0781913274877799047, -0.0903365513831131534,
            0.0560317850649417070, -0.0122200602431319916, -0.0147378888840894624,
            0.0175437183414442650, -0.00428819900201210177, -0.0110920129900025574,
            0.0182659646193858177, -0.0151388155198109841, 0.00619700872461139817)
    elseif name == "titan"
        return 1025.0, 975.0, (-19.4278696179288133, -14.4426058071846857, 2.69604850426823095,
            -0.336162464509102088, 0.0900778126111024535, -0.136262683746282226,
            0.232875124056659416, -0.226385345471503546, 0.118549197857537245,
            -0.0289659248379974324, 0.0165054791585107323, -0.0478606382329083910,
            0.0650600435627339962, -0.0503265704035375031, 0.0226198266668309125,
            -0.00645321729948887882, 0.00509967575456080676, -0.00943672721749253159,
            0.00948561230043162319, -0.00620859370479269826, 0.00120768474131884517)
    end
    throw(ArgumentError("No deprecated closed-form density fit is defined for planet $(planet.name)."))
end

function _edg_closed_form_density(planet, altitude_m::Float64)
    center_km, half_width_km, coefficients = _edg_closed_form_density_fit(planet)
    x = (altitude_m * 1e-3 - center_km) / half_width_km
    t_prev = 1.0
    exponent = coefficients[1]
    length(coefficients) == 1 && return exp(exponent)
    t_curr = x
    exponent += coefficients[2] * t_curr
    for j in 3:length(coefficients)
        t_next = 2.0 * x * t_curr - t_prev
        exponent += coefficients[j] * t_next
        t_prev, t_curr = t_curr, t_next
    end
    return exp(exponent)
end

function _edg_closed_form_heat_load_trajectory(
    config,
    p::ODEParams,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass::Float64,
    t::Float64,
    env,
    alpha_profile::Vector{Float64}=Float64[];
    prediction_step_s::Float64=0.1,
)
    planet = p.args.environment_model.planet
    duration = _edg_drag_passage_duration(config, p, pos, vel, mass)
    times = _edg_prediction_time_grid(duration, prediction_step_s)
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
    t_peri = max(0.5 * duration, eps(Float64))
    cost_3 = v0 * gamma0

    for j in eachindex(times)
        tau = times[j]
        altitude = altitude0 + cost_3 * (tau - tau^2 / (2.0 * t_peri))
        h[j] = altitude
        rho[j] = _edg_closed_form_density(planet, altitude)
        temperature[j] = planet.T
    end

    if isempty(alpha_profile)
        alpha_profile = fill(config.max_alpha_rad, n)
    elseif length(alpha_profile) > n
        alpha_profile = alpha_profile[1:n]
    elseif length(alpha_profile) < n
        alpha_profile = vcat(alpha_profile, fill(last(alpha_profile), n - length(alpha_profile)))
    end

    area_total = _edg_total_ref_area(spacecraft)
    area_panels = min(_edg_controlled_ref_area(spacecraft, config.controlled_panel_links), area_total)
    area_spacecraft = max(0.0, area_total - area_panels)
    T0 = max(Float64(planet.T), eps(Float64))
    S0 = v0 / sqrt(max(2.0 * planet.R * T0, eps(Float64)))
    cl_high, cd_high = _edg_weighted_aero_coefficients(spacecraft, T0, S0, config.max_alpha_rad)
    cl_low, cd_low = _edg_weighted_aero_coefficients(spacecraft, T0, S0, 0.0)
    cd_t = cd_low .+ alpha_profile .* (cd_high - cd_low) ./ (pi / 2)
    cl_t = cl_low .+ alpha_profile .* (cl_high - cl_low) ./ (pi / 2)
    cost_1 = rho .* cd_t .* area_total ./ (2.0 * mass) .* alpha_profile
    cost_2 = rho .* cl_t .* area_total ./ (2.0 * mass)
    f1, f2 = _edg_closed_form_f1_f2(planet, v0, gamma0, times, t_peri)
    epsilon = f1 .+
        f2 .* config.max_alpha_rad .* area_panels ./ area_total .+
        f2 .* (pi / 2) .* area_spacecraft ./ area_total
    k1 = cost_2 .+ (1.0 ./ (planet.Rp_e .+ h))
    k2 = (cost_1 .* cost_3) .* (1.0 .- times ./ t_peri)
    k3 = -planet.g_ref .- epsilon
    discriminant = (k2 ./ k1).^2 .- 4.0 .* (k3 ./ k1)
    speed_root = (k2 ./ k1 .- sqrt.(discriminant)) ./ 2.0
    speed_offset = v0 - speed_root[1]
    speed .= max.(speed_root .+ speed_offset, eps(Float64))
    gamma .= cost_3 .* (1.0 .- times ./ t_peri) ./ speed
    for j in eachindex(times)
        sound_speed = sqrt(max(0.0, planet.γ * planet.R * temperature[j]))
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
    _altitude_offset_m::Float64=0.0,
)
    planet = p.args.environment_model.planet
    area = _edg_total_ref_area(spacecraft)
    n = length(times)
    positions = Vector{SVector{3, Float64}}(undef, n)
    velocities = Vector{SVector{3, Float64}}(undef, n)
    positions[1] = pos0
    velocities[1] = vel0

    function acceleration(r, v, tau, alpha)
        gravity = _edg_prediction_gravity_acceleration(p, r, v, mass, t + tau)
        aero = _edg_targeting_aero_acceleration(config, p, spacecraft, r, v, mass, t + tau, alpha)
        return gravity + aero
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
        env = _edg_targeting_prediction_environment(p, r, v, t + times[j])
        radius = max(norm(r), eps(Float64))
        vel_norm = max(norm(v), eps(Float64))
        h[j] = env.altitude_m
        gamma[j] = asin(clamp(dot(r, v) / max(radius * vel_norm, eps(Float64)), -1.0, 1.0))
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
    )
end

function _edg_heat_load_costates(track, alpha_profile::Vector{Float64}, coeffs, mass::Float64, area::Float64, planet, scale_height::Float64, k::Float64)
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
    return (; lambda_switch, lambda_v, lambda_gamma, lambda_h)
end

function _edg_heat_load_lambdas(track, alpha_profile::Vector{Float64}, coeffs, mass::Float64, area::Float64, planet, scale_height::Float64, k::Float64)
    costates = _edg_heat_load_costates(track, alpha_profile, coeffs, mass, area, planet, scale_height, k)
    return costates.lambda_switch, costates.lambda_v
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
                min_alpha=_edg_constraint_min_alpha(config),
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

function _edg_heat_load_shooting_stage(
    config,
    p::ODEParams,
    spacecraft,
    mass::Float64,
    t::Float64,
    tau::Float64,
    y::SVector{9, Float64},
    coeffs,
    scale_height::Float64,
    k::Float64,
    alpha_past::Float64;
    heat_rate_control::Bool,
    structural_control::Bool,
    switch_smoothing::Float64=0.0,
)
    planet = p.args.environment_model.planet
    position = SVector{3, Float64}(y[1], y[2], y[3])
    velocity = SVector{3, Float64}(y[4], y[5], y[6])
    radius = max(norm(position), eps(Float64))
    speed = max(norm(velocity), eps(Float64))
    gamma = asin(clamp(dot(position, velocity) / (radius * speed), -1.0, 1.0))
    env = _edg_targeting_prediction_environment(p, position, velocity, t + tau)
    rho = env.rho
    area = _edg_total_ref_area(spacecraft)
    switch_denominator = max(abs(area * coeffs.cd_slope * pi), eps(Float64))
    lambda_switch = k * 2.0 * mass * speed / switch_denominator
    high_phase = y[7] >= lambda_switch
    raw_alpha = high_phase ? config.max_alpha_rad : config.min_alpha_rad
    high_weight = switch_smoothing > 0.0 ?
        0.5 * (1.0 + tanh((y[7] - lambda_switch) / switch_smoothing)) :
        (high_phase ? 1.0 : 0.0)
    base_alpha = config.min_alpha_rad + high_weight * (config.max_alpha_rad - config.min_alpha_rad)
    alpha = base_alpha > config.min_alpha_rad + 1e-8 ? _edg_targeting_constrained_alpha(
        config,
        p,
        spacecraft,
        env,
        base_alpha,
        alpha_past;
        heat_rate_control=heat_rate_control,
        structural_control=structural_control,
    ) : config.min_alpha_rad

    gravity_acceleration = _edg_prediction_gravity_acceleration(p, position, velocity, mass, t + tau)
    aerodynamic_acceleration = _edg_targeting_aero_acceleration(config, spacecraft, mass, alpha, env)
    acceleration = gravity_acceleration + aerodynamic_acceleration
    gravity = norm(gravity_acceleration)
    cd = coeffs.cd_low + alpha * coeffs.cd_slope
    lambda_v, lambda_gamma, lambda_h = y[7], y[8], y[9]

    lambda_v_dot =
        -3.0 * k * rho * speed^2 * alpha / pi +
        lambda_v * rho * area * cd * speed / mass -
        lambda_gamma * (rho * area * coeffs.cl_low / (2.0 * mass) + gravity / speed^2 + 1.0 / radius) -
        lambda_h * gamma
    lambda_gamma_dot = lambda_v * gravity - lambda_h * speed
    lambda_h_dot =
        k * rho * speed^3 * alpha / (pi * scale_height) -
        lambda_v * (rho * area * cd * speed^2 / (2.0 * mass * scale_height) + 2.0 * gravity * gamma / radius) +
        lambda_gamma * (rho * area * coeffs.cl_low * speed / (2.0 * mass * scale_height) -
                        2.0 * gravity / (radius * speed) + speed / radius^2)

    derivative = SVector{9, Float64}(
        velocity[1], velocity[2], velocity[3],
        acceleration[1], acceleration[2], acceleration[3],
        lambda_v_dot, lambda_gamma_dot, lambda_h_dot,
    )
    return derivative, raw_alpha, alpha
end

function _edg_integrate_heat_load_shooting(
    config,
    p::ODEParams,
    spacecraft,
    initial_position::SVector{3, Float64},
    initial_velocity::SVector{3, Float64},
    mass::Float64,
    t::Float64,
    times::Vector{Float64},
    coeffs,
    scale_height::Float64,
    k::Float64,
    initial_costates::SVector{3, Float64};
    heat_rate_control::Bool,
    structural_control::Bool,
    switch_smoothing::Float64=0.0,
)
    n = length(times)
    states = Vector{SVector{9, Float64}}(undef, n)
    raw_profile = Vector{Float64}(undef, n)
    alpha_profile = Vector{Float64}(undef, n)
    states[1] = SVector{9, Float64}(initial_position..., initial_velocity..., initial_costates...)
    alpha_past = config.max_alpha_rad

    for j in 1:(n - 1)
        dt = times[j + 1] - times[j]
        y = states[j]
        k1, raw_alpha, alpha = _edg_heat_load_shooting_stage(
            config, p, spacecraft, mass, t, times[j], y, coeffs, scale_height, k, alpha_past;
            heat_rate_control=heat_rate_control, structural_control=structural_control,
            switch_smoothing=switch_smoothing,
        )
        k2, _, _ = _edg_heat_load_shooting_stage(
            config, p, spacecraft, mass, t, times[j] + 0.5dt, y + 0.5dt * k1,
            coeffs, scale_height, k, alpha;
            heat_rate_control=heat_rate_control, structural_control=structural_control,
            switch_smoothing=switch_smoothing,
        )
        k3, _, _ = _edg_heat_load_shooting_stage(
            config, p, spacecraft, mass, t, times[j] + 0.5dt, y + 0.5dt * k2,
            coeffs, scale_height, k, alpha;
            heat_rate_control=heat_rate_control, structural_control=structural_control,
            switch_smoothing=switch_smoothing,
        )
        k4, _, _ = _edg_heat_load_shooting_stage(
            config, p, spacecraft, mass, t, times[j + 1], y + dt * k3,
            coeffs, scale_height, k, alpha;
            heat_rate_control=heat_rate_control, structural_control=structural_control,
            switch_smoothing=switch_smoothing,
        )
        states[j + 1] = y + dt * (k1 + 2.0k2 + 2.0k3 + k4) / 6.0
        raw_profile[j] = raw_alpha
        alpha_profile[j] = alpha
        alpha_past = alpha
    end
    _, raw_profile[end], alpha_profile[end] = _edg_heat_load_shooting_stage(
        config, p, spacecraft, mass, t, times[end], states[end], coeffs, scale_height, k, alpha_past;
        heat_rate_control=heat_rate_control, structural_control=structural_control,
        switch_smoothing=switch_smoothing,
    )
    return states, raw_profile, alpha_profile
end

function _edg_heat_load_shooting_profile(
    config,
    p::ODEParams,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass::Float64,
    t::Float64,
    seed_track,
    seed_profile::Vector{Float64},
    coeffs,
    scale_height::Float64,
    k::Float64;
    heat_rate_control::Bool,
    structural_control::Bool,
    initial_costates_guess::Union{Nothing, SVector{3, Float64}}=nothing,
)
    area = _edg_total_ref_area(spacecraft)
    initial_guess = _edg_heat_load_costates(
        seed_track, seed_profile, coeffs, mass, area,
        p.args.environment_model.planet, scale_height, k,
    )
    initial_position = seed_track.positions[1]
    initial_velocity = seed_track.velocities[1]
    initial_speed = norm(initial_velocity)
    costate_scale = SVector{3, Float64}(
        max(abs(initial_guess.lambda_v[1]), initial_speed, 500.0),
        max(abs(initial_guess.lambda_gamma[1]), 80_000.0),
        max(abs(initial_guess.lambda_h[1]), p.args.environment_model.planet.μ / norm(pos)^2, 9.0),
    )
    backward_guess = [
        initial_guess.lambda_v[1] / costate_scale[1],
        initial_guess.lambda_gamma[1] / costate_scale[2],
        initial_guess.lambda_h[1] / costate_scale[3],
    ]
    last_solution = Ref{Any}(nothing)
    switch_smoothing = Ref(0.0)

    function shooting_residual!(residual, z)
        initial_costates = SVector{3, Float64}(z[1], z[2], z[3]) .* costate_scale
        states, raw_profile, alpha_profile = _edg_integrate_heat_load_shooting(
            config, p, spacecraft, initial_position, initial_velocity,
            mass, t, seed_track.time, coeffs,
            scale_height, k, initial_costates;
            heat_rate_control=heat_rate_control, structural_control=structural_control,
            switch_smoothing=switch_smoothing[],
        )
        final_state = states[end]
        final_radius = norm(SVector{3, Float64}(final_state[1], final_state[2], final_state[3]))
        final_speed = norm(SVector{3, Float64}(final_state[4], final_state[5], final_state[6]))
        terminal_h = p.args.environment_model.planet.μ / max(final_radius^2, eps(Float64))
        residual[1] = (final_state[7] - final_speed) / costate_scale[1]
        residual[2] = final_state[8] / costate_scale[2]
        residual[3] = (final_state[9] - terminal_h) / costate_scale[3]
        last_solution[] = (; states, raw_profile, alpha_profile, initial_costates, residual=copy(residual))
        return nothing
    end

    function shooting_jacobian!(jacobian, z)
        for column in eachindex(z)
            lower = copy(z)
            upper = copy(z)
            relative_step = clamp(0.1 * switch_smoothing[] / costate_scale[1], 1e-5, 1e-3)
            step = relative_step * max(abs(z[column]), 1.0)
            lower[column] -= step
            upper[column] += step
            lower_residual = zeros(3)
            upper_residual = zeros(3)
            shooting_residual!(lower_residual, lower)
            shooting_residual!(upper_residual, upper)
            jacobian[:, column] .= (upper_residual .- lower_residual) ./ (2.0 * step)
        end
        return nothing
    end

    legacy_guess = [-500.0, 80_000.0, 9.0] ./ collect(costate_scale)
    best_result = nothing
    best_residual_norm = Inf
    smoothing_sequence = (
        costate_scale[1],
        0.1 * costate_scale[1],
        0.01 * costate_scale[1],
    )
    initial_guesses = Vector{Vector{Float64}}()
    if initial_costates_guess !== nothing
        push!(initial_guesses, collect(initial_costates_guess ./ costate_scale))
    end
    push!(initial_guesses, backward_guess, legacy_guess)
    for initial in initial_guesses
        guess = collect(initial)
        continuation_succeeded = true
        for smoothing in smoothing_sequence
            switch_smoothing[] = smoothing
            solution = NLsolve.nlsolve(
                shooting_residual!, shooting_jacobian!, guess;
                method=:trust_region, ftol=1e-9, xtol=1e-9, iterations=100, show_trace=false,
            )
            if !all(isfinite, solution.zero)
                continuation_succeeded = false
                break
            end
            guess = collect(solution.zero)
            shooting_residual!(zeros(3), guess)
            candidate = last_solution[]
            candidate_norm = norm(candidate.residual, Inf)
            if !isfinite(candidate_norm)
                continuation_succeeded = false
                break
            end
        end
        if continuation_succeeded
            candidate = last_solution[]
            candidate_norm = norm(candidate.residual, Inf)
            if candidate_norm < best_residual_norm
                best_result = candidate
                best_residual_norm = candidate_norm
            end
            candidate_norm <= _EDG_HEAT_LOAD_SHOOTING_RESIDUAL_TOL && break
        end
    end
    if best_result === nothing || best_residual_norm > _EDG_HEAT_LOAD_SHOOTING_RESIDUAL_TOL
        residual = best_result === nothing ? nothing : best_result.residual
        throw(ErrorException(
            "Heat-load TPBVP shooting failed to satisfy the terminal costates " *
            "(scaled infinity norm = $(best_residual_norm), residual = $(residual))."
        ))
    end
    return best_result.states, best_result.raw_profile, best_result.alpha_profile,
        best_residual_norm, best_result.initial_costates
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

function _edg_profile_heat_load(config, p::ODEParams, track, alpha_profile::Vector{Float64}; heat_rate_control::Bool)
    qdot = _edg_profile_heat_rates(config, p, track, alpha_profile; heat_rate_control=heat_rate_control)
    return sum(qdot) * last(track.time) / length(track.time)
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
    shooting_guess::Union{Nothing, SVector{3, Float64}}=nothing,
)
    planet = p.args.environment_model.planet
    area = _edg_total_ref_area(spacecraft)
    scale_height = _edg_heat_load_scale_height(p)
    prediction_step_s = config.heat_load_switch_solver == :tpbvp_integration ? 1.0 : 0.1
    base_track = _edg_closed_form_heat_load_trajectory(
        config,
        p,
        spacecraft,
        pos,
        vel,
        mass,
        t,
        env;
        prediction_step_s=prediction_step_s,
    )
    high_profile = _edg_constrained_heat_load_alpha_profile(
        config,
        p,
        spacecraft,
        controlled_panel_links,
        base_track,
        fill(config.max_alpha_rad, length(base_track.time));
        heat_rate_control=heat_rate_control,
        structural_control=structural_control,
    )

    if config.heat_load_switch_solver == :tpbvp_integration
        seed_track = _edg_integrated_heat_load_trajectory(
            config, p, spacecraft, pos, vel, mass, t, base_track.time, high_profile
        )
        _, raw_profile, _, _, solved_costates = _edg_heat_load_shooting_profile(
            config, p, spacecraft, pos, vel, mass, t, seed_track, high_profile,
            coeffs, scale_height, k;
            heat_rate_control=heat_rate_control,
            structural_control=structural_control,
            initial_costates_guess=shooting_guess,
        )
        applied_profile = _edg_first_two_switch_alpha_profile(config, raw_profile, high_profile)
        track = _edg_integrated_heat_load_trajectory(
            config, p, spacecraft, pos, vel, mass, t, base_track.time, applied_profile
        )
        high_profile = _edg_constrained_heat_load_alpha_profile(
            config, p, spacecraft, controlled_panel_links, track,
            fill(config.max_alpha_rad, length(track.time));
            heat_rate_control=heat_rate_control,
            structural_control=structural_control,
        )
        applied_profile = _edg_first_two_switch_alpha_profile(config, raw_profile, high_profile)
        track = _edg_integrated_heat_load_trajectory(
            config, p, spacecraft, pos, vel, mass, t, base_track.time, applied_profile
        )
        return track, raw_profile, applied_profile, solved_costates
    end

    track = base_track
    raw_profile = fill(config.max_alpha_rad, length(track.time))
    applied_profile = copy(high_profile)
    previous_heat_load = 0.0
    heat_load_change = Inf
    iterations = 0
    while heat_load_change > 1e-3 && iterations < 100
        raw_profile, _, _ = _edg_heat_load_alpha_profile(
            config, track, coeffs, mass, area, planet, scale_height, k, applied_profile
        )
        high_profile = _edg_constrained_heat_load_alpha_profile(
            config, p, spacecraft, controlled_panel_links, track,
            fill(config.max_alpha_rad, length(track.time));
            heat_rate_control=heat_rate_control,
            structural_control=structural_control,
        )
        applied_profile = _edg_first_two_switch_alpha_profile(config, raw_profile, high_profile)
        track = _edg_closed_form_heat_load_trajectory(
            config, p, spacecraft, pos, vel, mass, t, env, applied_profile
        )
        predicted_heat_load = _edg_profile_heat_load(config, p, track, applied_profile; heat_rate_control=false)
        heat_load_change = abs(predicted_heat_load - previous_heat_load)
        previous_heat_load = predicted_heat_load
        iterations += 1
    end
    return track, raw_profile, applied_profile, nothing
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
    target_load = config.heat_load_limit_j_cm2
    planet_name = lowercase(String(p.args.environment_model.planet.name))
    if config.heat_load_switch_solver == :tpbvp_integration
        low_k, high_k = planet_name == "mars" ? (3.3, 10.0) :
            planet_name == "venus" ? (0.1, 100.0) :
            planet_name == "earth" ? (1.0, 10.0) : (1e-5, 10.0)
    else
        low_k = 0.0
        high_k = planet_name == "venus" ? 100.0 : 10.0
    end

    last_track = Ref{Any}(nothing)
    last_profile = Ref{Vector{Float64}}(Float64[])
    shooting_solutions = Vector{Tuple{Float64, SVector{3, Float64}}}()

    function residual(k_root)
        k = config.heat_load_switch_solver == :tpbvp_integration ? k_root / 100.0 : k_root
        shooting_guess = if isempty(shooting_solutions)
            nothing
        else
            nearest = argmin(abs(item[1] - k_root) for item in shooting_solutions)
            shooting_solutions[nearest][2]
        end
        track, raw_profile, profile, solved_costates = _edg_heat_load_profile_for_k(
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
            shooting_guess,
        )
        solved_costates === nothing || push!(shooting_solutions, (Float64(k_root), solved_costates))
        last_track[] = track
        last_profile[] = profile
        predicted_load = heat_load_j_cm2 + _edg_profile_heat_load(
            config,
            p,
            track,
            profile;
            heat_rate_control=false,
        )
        return predicted_load - target_load
    end

    f_high = residual(high_k)
    f_low = residual(low_k)
    if f_high * f_low < 0.0
        bracket = if config.heat_load_switch_solver == :closed_form
            (0.0, 0.1)
        else
            (low_k, high_k)
        end
        root_tolerance = config.heat_load_switch_solver == :tpbvp_integration ? 1e-3 : 1e-5
        k_sol = Roots.find_zero(residual, bracket, Roots.Brent(); rtol=root_tolerance)
        residual(k_sol)
    elseif f_high < 0.0
        return (t, t)
    elseif f_low > 0.0
        return config.heat_load_switch_solver == :closed_form ?
            (t, t + last(last_track[].time) / 2.0) : (t, t + 1_000.0)
    else
        return (Inf, Inf)
    end

    track = last_track[]
    profile = copy(last_profile[])
    window = _edg_low_alpha_switch_window(config, t, track, profile)
    if all(isfinite, window) && window[1] < window[2]
        return window
    end
    return (Inf, Inf)
end

function _edg_recompute_second_heat_load_switch(
    config,
    p::ODEParams,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass_state::Float64,
    env,
    heat_load_j_cm2::Float64,
    t::Float64,
    current_switches::NTuple{2, Float64};
    passage_entry_time_s::Float64,
    reevaluation_mode::Int,
    heat_rate_control::Bool,
    structural_control::Bool,
)
    mass = _edg_predict_mass(spacecraft, mass_state)
    duration = _edg_drag_passage_duration(config, p, pos, vel, mass)
    times = _edg_prediction_time_grid(duration, config.heat_load_switch_solver == :tpbvp_integration ? 1.0 : 0.1)
    target_load = config.heat_load_limit_j_cm2

    function residual(second_switch::Float64)
        if config.heat_load_switch_solver == :closed_form
            base_track = _edg_closed_form_heat_load_trajectory(config, p, spacecraft, pos, vel, mass, t, env)
            high_profile = _edg_constrained_heat_load_alpha_profile(
                config,
                p,
                spacecraft,
                config.controlled_panel_links,
                base_track,
                fill(config.max_alpha_rad, length(base_track.time));
                heat_rate_control=heat_rate_control,
                structural_control=structural_control,
            )
            profile = [t + base_track.time[j] <= second_switch ? config.min_alpha_rad : high_profile[j]
                       for j in eachindex(base_track.time)]
            track = _edg_closed_form_heat_load_trajectory(config, p, spacecraft, pos, vel, mass, t, env, profile)
            future_load = _edg_profile_heat_load(config, p, track, profile; heat_rate_control=false)
        else
            base_track = _edg_integrated_heat_load_trajectory(
                config,
                p,
                spacecraft,
                pos,
                vel,
                mass,
                t,
                times,
                fill(config.max_alpha_rad, length(times)),
                env.altitude_m - (norm(pos) - p.args.environment_model.planet.Rp_e),
            )
            high_profile = _edg_constrained_heat_load_alpha_profile(
                config,
                p,
                spacecraft,
                config.controlled_panel_links,
                base_track,
                fill(config.max_alpha_rad, length(times));
                heat_rate_control=heat_rate_control,
                structural_control=structural_control,
            )
            profile = [t + times[j] <= second_switch ? config.min_alpha_rad : high_profile[j] for j in eachindex(times)]
            track = _edg_integrated_heat_load_trajectory(
                config,
                p,
                spacecraft,
                pos,
                vel,
                mass,
                t,
                times,
                profile,
                env.altitude_m - (norm(pos) - p.args.environment_model.planet.Rp_e),
            )
            future_load = _edg_profile_heat_load(config, p, track, profile; heat_rate_control=false)
        end
        return heat_load_j_cm2 + future_load - target_load
    end

    lower = t
    upper = t + 200.0
    lower_error = residual(lower)
    current_error = residual(current_switches[2])
    if config.heat_load_switch_solver == :closed_form
        x_tol = reevaluation_mode == 1 ? 0.1 : 0.01
        prediction_difference = abs(lower_error - current_error)
        if abs(current_error) <= x_tol &&
                !(prediction_difference < 0.5 && abs(t - current_switches[2]) > 20.0)
            return current_switches
        end
    else
        if abs(current_error) <= 0.01
            return current_switches
        elseif lower_error < 0.0
            return (current_switches[1], max(t, current_switches[1]))
        end
    end
    upper_error = residual(upper)
    if !(isfinite(lower_error) && isfinite(upper_error) && signbit(lower_error) != signbit(upper_error))
        return current_switches
    end
    second_switch = try
        Roots.find_zero(residual, (lower, upper), Roots.Brent())
    catch
        return current_switches
    end
    return (current_switches[1], max(second_switch, current_switches[1]))
end

function _edg_heat_load_security_required(
    config,
    p::ODEParams,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass_state::Float64,
    env,
    heat_load_j_cm2::Float64,
    t::Float64,
)
    mass = _edg_predict_mass(spacecraft, mass_state)
    track = _edg_closed_form_heat_load_trajectory(config, p, spacecraft, pos, vel, mass, t, env)
    low_profile = fill(config.min_alpha_rad, length(track.time))
    qdot = _edg_profile_heat_rates(config, p, track, low_profile; heat_rate_control=false)
    qdot[track.time .<= Float64(p.args.environment_model.planet.T)] .= 0.0
    remaining_load = sum(qdot) * last(track.time) / length(track.time)
    return heat_load_j_cm2 + remaining_load > config.heat_load_limit_j_cm2,
        t + last(track.time)
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
