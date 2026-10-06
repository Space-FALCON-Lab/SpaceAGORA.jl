function _spaceagora_rpo_modules()
    spaceagora = _load_spaceagora!(load_gramsuite=false)
    simulation_model = Base.invokelatest(getproperty, spaceagora, :SimulationModel)
    return (
        guidance=Base.invokelatest(getproperty, simulation_model, :GuidanceHooks),
        navigation=Base.invokelatest(getproperty, simulation_model, :NavigationHooks),
        control=Base.invokelatest(getproperty, simulation_model, :ControlHooks),
    )
end

"""Build a default Gateway/CubeSat scenario using SpaceAGORA's station asset."""
function build_rpo_hypr_rl_scenario(;
    start_rtn=[-8.0, -4.0, 2.0],
    goal_rtn=[5.0, 0.0, 0.0],
    station_asset::Symbol=:gateway,
    station_points::Int=2_000,
    station_seed::Int=740,
    station_keepout_radius_m::Real=0.25,
    initial_attitude_rtn_to_body=[0.0, 0.0, 0.0, 1.0],
    final_attitude_rtn_to_body=[0.0, 0.0, 0.0, 1.0],
    pso_config=nothing,
    tracking_settings=nothing,
    rrt_settings=nothing,
)
    spaceagora = _load_spaceagora!(load_gramsuite=false)
    modules = _spaceagora_rpo_modules()
    pointcloud_loader = Base.invokelatest(
        getproperty, spaceagora, :load_rpo_station_cad_pointcloud,
    )
    points = Base.invokelatest(
        pointcloud_loader, station_asset;
        n_points=station_points, rng=MersenneTwister(station_seed),
    )
    station_geometry = Base.invokelatest(
        getproperty, modules.navigation, :RPOStationGeometry,
    )
    station = Base.invokelatest(
        station_geometry, points;
        keepout_radius_m=Float64(station_keepout_radius_m),
        name=String(station_asset),
    )
    reference_geometry = Base.invokelatest(
        getproperty, modules.navigation, :RPOReferenceGeometry,
    )
    geometry = Base.invokelatest(
        reference_geometry, station,
    )
    return RPOHyPRRLScenario(
        start_rtn=Vector{Float64}(start_rtn),
        goal_rtn=Vector{Float64}(goal_rtn),
        geometry=geometry,
        pso_config=pso_config,
        tracking_settings=tracking_settings,
        rrt_settings=rrt_settings,
        initial_attitude_rtn_to_body=Vector{Float64}(initial_attitude_rtn_to_body),
        final_attitude_rtn_to_body=Vector{Float64}(final_attitude_rtn_to_body),
    )
end

function _rpo_spaceagora_settings(scenario::RPOHyPRRLScenario,
                                  config::RPOHyPRRLConfig)
    modules = _spaceagora_rpo_modules()
    pso_config = scenario.pso_config === nothing ?
        Base.invokelatest(
            getproperty(modules.guidance, :rpo_740_mpc_final_pso_config);
            safe_distance_m=config.safe_distance_m,
        ) : scenario.pso_config
    pso_config = Base.invokelatest(
        getproperty(modules.guidance, :rpo_pso_config), pso_config;
        safe_distance_m=config.safe_distance_m,
        w_len=0.0, w_fuel=0.0,
    )
    tracking = scenario.tracking_settings === nothing ?
        Base.invokelatest(getproperty(modules.guidance, :RPOLQMPCTrackingSettings)) :
        scenario.tracking_settings
    rrt = scenario.rrt_settings === nothing ?
        Base.invokelatest(getproperty(modules.guidance, :RPORRTConnectSettings)) :
        scenario.rrt_settings
    return modules, pso_config, tracking, rrt
end

function _rpo_fuel_wheel_objective(components, propellant_used_kg::Real,
                                   wheel_energy_j::Real,
                                   config::RPOHyPRRLConfig)
    fuel_reference_kg = max(Float64(get(components, :fuel_ref, 0.0)), 1.0e-12)
    normalized_fuel = Float64(propellant_used_kg) / fuel_reference_kg
    wheel_energy_reference_j = sum(
        config.reaction_wheel_max_momentum_nms^2 / (2.0 * inertia)
        for inertia in config.reaction_wheel_inertia_kgm2
    )
    normalized_wheel_energy = Float64(wheel_energy_j) /
        max(wheel_energy_reference_j, 1.0e-12)
    return (
        total=config.fuel_weight * normalized_fuel +
              config.wheel_weight * normalized_wheel_energy,
        normalized_fuel=normalized_fuel,
        normalized_wheel_energy=normalized_wheel_energy,
        fuel_reference_kg=fuel_reference_kg,
        wheel_energy_reference_j=wheel_energy_reference_j,
    )
end

function _rpo_original_clearance_penalty(components)
    # RPOHyPRRL evaluators zero the length and proxy-fuel weights, leaving the
    # original HyPR w_obs * J_obs sigmoid clearance term in .total.
    penalty = Float64(get(components, :total, Inf))
    return isfinite(penalty) ? max(0.0, penalty) : Inf
end

function _rpo_clearance_penalty_components(components, violation_count::Integer)
    penalty = _rpo_original_clearance_penalty(components)
    return merge(
        components,
        (
            total=penalty,
            violation_count=Int(violation_count),
            clearance_penalty=penalty,
            retimed_feasible=false,
        ),
    )
end

@inline function _rpo_clearance_feasibility_tol_m(pso_config)
    return pso_config !== nothing &&
        hasproperty(pso_config, :clearance_feasibility_tol_m) ?
        Float64(getproperty(pso_config, :clearance_feasibility_tol_m)) : 1.0e-3
end

@inline function _rpo_clearance_is_feasible(clearance::Real,
                                             safe_distance_m::Real,
                                             pso_config)
    return Float64(clearance) + _rpo_clearance_feasibility_tol_m(pso_config) >=
        Float64(safe_distance_m)
end

@inline function _rpo_pso_objective_components(
    components;
    total::Float64=Inf,
    violation_count::Int=Int(components.violation_count),
    clearance_penalty::Float64=0.0,
    propellant_used_kg::Float64=0.0,
    normalized_fuel::Float64=0.0,
    wheel_energy_j::Float64=0.0,
    normalized_wheel_energy::Float64=0.0,
    duration_s::Float64=0.0,
    feasible::Bool=false,
    cutoff_exceeded::Bool=false,
)
    return (
        total=total,
        J_len=Float64(components.J_len),
        J_len_norm=Float64(components.J_len_norm),
        J_obs=Float64(components.J_obs),
        J_fuel=propellant_used_kg,
        J_fuel_norm=normalized_fuel,
        min_clearance=Float64(components.min_clearance),
        violation_count=violation_count,
        cutoff_exceeded=cutoff_exceeded,
        len_ref=Float64(components.len_ref),
        fuel_ref=Float64(components.fuel_ref),
        clearance_penalty=clearance_penalty,
        retimed_propellant_used_kg=propellant_used_kg,
        retimed_wheel_energy_j=wheel_energy_j,
        retimed_wheel_energy_norm=normalized_wheel_energy,
        retimed_duration_s=duration_s,
        retimed_feasible=feasible,
    )
end

"""Typed PSO objective adapter for streamed retiming and actuator evaluation."""
struct RPOHyPRRLPSOObjectiveEvaluator{TTracking, TPrepare, TInterpolate, TAllocator}
    config::RPOHyPRRLConfig
    tracking_settings::TTracking
    initial_attitude_rtn_to_body::Vector{Float64}
    final_attitude_rtn_to_body::Vector{Float64}
    prepare_candidate::TPrepare
    interpolate_along_path!::TInterpolate
    allocator::TAllocator
end

function RPOHyPRRLPSOObjectiveEvaluator(config::RPOHyPRRLConfig,
                                        scenario::RPOHyPRRLScenario)
    modules, _, tracking, _ = _rpo_spaceagora_settings(scenario, config)
    prepare_candidate = getproperty(
        modules.guidance, :rpo_prepare_retimed_candidate,
    )
    interpolate_along_path! = getproperty(
        modules.guidance, :rpo_interpolate_along_path!,
    )
    return RPOHyPRRLPSOObjectiveEvaluator(
        config,
        tracking,
        copy(scenario.initial_attitude_rtn_to_body),
        copy(scenario.final_attitude_rtn_to_body),
        prepare_candidate,
        interpolate_along_path!,
        _rpo_thruster_allocator(config),
    )
end

function (evaluator::RPOHyPRRLPSOObjectiveEvaluator)(
    control_points_rtn,
    geometry,
    pso_config,
    safe_distance_m::Real,
    cost_cutoff::Real,
)
    prepared = evaluator.prepare_candidate(
        control_points_rtn,
        geometry,
        pso_config;
        safe_distance_m=Float64(safe_distance_m),
        retime_dt_s=evaluator.tracking_settings.dt_s,
        w_len=0.0,
        w_fuel=0.0,
        compute_fuel_proxy=false,
        retime_mean_motion_radps=evaluator.tracking_settings.mean_motion_radps,
        retime_command_limit_mps2=_rpo_guaranteed_acceleration_limit_mps2(
            evaluator.config,
            evaluator.tracking_settings.dt_s,
            evaluator.allocator,
        ),
    )
    components = prepared.components
    violation_count = !_rpo_clearance_is_feasible(
        components.min_clearance, safe_distance_m, pso_config,
    ) ?
        max(1, components.violation_count) : components.violation_count
    if violation_count > 0
        penalty = _rpo_original_clearance_penalty(components)
        return _rpo_pso_objective_components(
            components;
            total=penalty,
            violation_count=Int(violation_count),
            clearance_penalty=penalty,
        )
    end
    if !isfinite(components.total)
        return _rpo_pso_objective_components(
            components; violation_count=Int(violation_count),
        )
    end

    prepared.profile === nothing && return _rpo_pso_objective_components(
        components; violation_count=Int(violation_count),
    )
    accounting = _stream_retimed_feedforward_objective(
        prepared.profile,
        evaluator.tracking_settings,
        evaluator.initial_attitude_rtn_to_body,
        evaluator.final_attitude_rtn_to_body,
        evaluator.config,
        evaluator.allocator,
        components.fuel_ref,
        Float64(cost_cutoff),
        evaluator.interpolate_along_path!,
        optimize_pointing=true,
    )
    accounting.cutoff_exceeded && return _rpo_pso_objective_components(
        components;
        total=Inf,
        violation_count=Int(violation_count),
        clearance_penalty=0.0,
        propellant_used_kg=accounting.propellant_used_kg,
        normalized_fuel=accounting.normalized_fuel,
        wheel_energy_j=accounting.wheel_energy_j,
        normalized_wheel_energy=accounting.normalized_wheel_energy,
        duration_s=accounting.duration_s,
        cutoff_exceeded=true,
    )
    final_point = view(prepared.profile.samples, :, size(prepared.profile.samples, 2))
    goal_point = view(control_points_rtn, :, size(control_points_rtn, 2))
    final_error = norm(final_point - goal_point)
    wheel_feasible = accounting.wheel_peak_momentum_nms <=
        evaluator.config.reaction_wheel_max_momentum_nms + 1.0e-12
    feasible = final_error <= evaluator.tracking_settings.final_position_tol_m &&
        wheel_feasible && isfinite(accounting.total)
    return _rpo_pso_objective_components(
        components;
        total=feasible ? accounting.total : Inf,
        violation_count=Int(violation_count),
        clearance_penalty=0.0,
        propellant_used_kg=accounting.propellant_used_kg,
        normalized_fuel=accounting.normalized_fuel,
        wheel_energy_j=accounting.wheel_energy_j,
        normalized_wheel_energy=accounting.normalized_wheel_energy,
        duration_s=accounting.duration_s,
        feasible=feasible,
    )
end

"""Typed PSO adapter for the geometric objective with retimed HCW proxy fuel."""
struct RPOTimingAwareProxyPSOObjectiveEvaluator{TTracking, TPrepare, TInterpolate}
    config::RPOHyPRRLConfig
    tracking_settings::TTracking
    prepare_candidate::TPrepare
    interpolate_along_path!::TInterpolate
end

function RPOTimingAwareProxyPSOObjectiveEvaluator(
    config::RPOHyPRRLConfig,
    scenario::RPOHyPRRLScenario,
)
    modules, _, tracking, _ = _rpo_spaceagora_settings(scenario, config)
    return RPOTimingAwareProxyPSOObjectiveEvaluator(
        config,
        tracking,
        getproperty(modules.guidance, :rpo_prepare_retimed_candidate),
        getproperty(modules.guidance, :rpo_interpolate_along_path!),
    )
end

function (evaluator::RPOTimingAwareProxyPSOObjectiveEvaluator)(
    control_points_rtn,
    geometry,
    pso_config,
    safe_distance_m::Real,
    cost_cutoff::Real,
)
    prepared = evaluator.prepare_candidate(
        control_points_rtn,
        geometry,
        pso_config;
        safe_distance_m=Float64(safe_distance_m),
        retime_dt_s=evaluator.tracking_settings.dt_s,
        w_len=pso_config.w_len,
        w_obs=pso_config.w_obs,
        w_fuel=0.0,
        compute_fuel_proxy=false,
        retime_mean_motion_radps=evaluator.tracking_settings.mean_motion_radps,
        retime_command_limit_mps2=_rpo_guaranteed_acceleration_limit_mps2(
            evaluator.config,
            evaluator.tracking_settings.dt_s,
        ),
    )
    components = prepared.components
    violation_count = !_rpo_clearance_is_feasible(
        components.min_clearance, safe_distance_m, pso_config,
    ) ?
        max(1, components.violation_count) : components.violation_count
    if violation_count > 0
        penalty = _rpo_original_clearance_penalty(components)
        return _rpo_pso_objective_components(
            components;
            total=penalty,
            violation_count=Int(violation_count),
            clearance_penalty=penalty,
        )
    end
    if !isfinite(components.total)
        return _rpo_pso_objective_components(
            components; violation_count=Int(violation_count),
        )
    end
    if isfinite(cost_cutoff) && components.total > Float64(cost_cutoff)
        return _rpo_pso_objective_components(
            components;
            violation_count=Int(violation_count),
            cutoff_exceeded=true,
        )
    end
    prepared.profile === nothing && return _rpo_pso_objective_components(
        components; violation_count=Int(violation_count),
    )

    accounting = _stream_retimed_timing_aware_proxy_objective(
        prepared.profile,
        evaluator.tracking_settings.mean_motion_radps,
        pso_config.mass_kg,
        pso_config.isp_s,
        pso_config.g0_mps2,
        components.fuel_ref,
        pso_config.w_fuel,
        components.total,
        Float64(cost_cutoff),
        evaluator.interpolate_along_path!,
    )
    accounting.cutoff_exceeded && return _rpo_pso_objective_components(
        components;
        total=Inf,
        violation_count=Int(violation_count),
        propellant_used_kg=accounting.propellant_used_kg,
        normalized_fuel=accounting.normalized_fuel,
        duration_s=accounting.duration_s,
        cutoff_exceeded=true,
    )
    final_point = view(prepared.profile.samples, :, size(prepared.profile.samples, 2))
    goal_point = view(control_points_rtn, :, size(control_points_rtn, 2))
    final_error = norm(final_point - goal_point)
    feasible = final_error <= evaluator.tracking_settings.final_position_tol_m &&
        isfinite(accounting.total)
    return _rpo_pso_objective_components(
        components;
        total=feasible ? accounting.total : Inf,
        violation_count=Int(violation_count),
        propellant_used_kg=accounting.propellant_used_kg,
        normalized_fuel=accounting.normalized_fuel,
        duration_s=accounting.duration_s,
        feasible=feasible,
    )
end

function _rpo_endpoint_clearance_m(point, geometry)
    geometry === nothing && return Inf
    modules = _spaceagora_rpo_modules()
    return Float64(Base.invokelatest(
        getproperty(modules.navigation, :rpo_clearance_distance_to_station),
        point,
        geometry,
    ))
end

struct RPOHyPRRLEndpointSampler
    base::RPOHyPRRLScenario
    station_triangles::Matrix{Float64}
    triangle_cdf::Vector{Float64}
    endpoint_min_clearance_m::Float64
    endpoint_max_clearance_m::Float64
    min_separation_m::Float64
    surrounded_max_distance_m::Float64
    max_sampling_tries::Int
end

function _rpo_hypr_rl_triangle_cdf(triangles::AbstractMatrix)
    size(triangles, 1) == 3 ||
        throw(ArgumentError("station triangles must have three rows"))
    size(triangles, 2) % 3 == 0 ||
        throw(ArgumentError("station mesh contains an incomplete triangle"))
    n_triangles = size(triangles, 2) ÷ 3
    areas = zeros(Float64, n_triangles)
    for triangle_index in 1:n_triangles
        base = 3 * (triangle_index - 1)
        v1 = view(triangles, :, base + 1)
        v2 = view(triangles, :, base + 2)
        v3 = view(triangles, :, base + 3)
        areas[triangle_index] = 0.5 * norm(cross(v2 - v1, v3 - v1))
    end
    total_area = sum(areas)
    total_area > 0.0 || throw(ArgumentError("station mesh has zero triangle area"))
    return cumsum(areas) ./ total_area
end

function build_rpo_hypr_rl_endpoint_sampler(
    base::RPOHyPRRLScenario;
    station_asset::Symbol=:gateway,
    safe_distance_m::Real,
    endpoint_clearance_margin_m::Real=0.05,
    endpoint_max_clearance_m::Real=1.0,
    min_separation_m::Real=1.5,
    surrounded_max_distance_m::Real=2.0,
    max_sampling_tries::Integer=4_000,
)
    spaceagora = _load_spaceagora!(load_gramsuite=false)
    loader = Base.invokelatest(
        getproperty, spaceagora, :load_rpo_station_cad_triangles,
    )
    triangles = Matrix{Float64}(Base.invokelatest(loader, station_asset))
    minimum_clearance = Float64(safe_distance_m) +
        max(Float64(endpoint_clearance_margin_m), 0.0)
    maximum_clearance = max(Float64(endpoint_max_clearance_m), minimum_clearance)
    minimum_separation = Float64(min_separation_m)
    surrounded_distance = Float64(surrounded_max_distance_m)
    tries = Int(max_sampling_tries)
    minimum_separation >= 0.0 ||
        throw(ArgumentError("minimum endpoint separation must be nonnegative"))
    surrounded_distance > 0.0 ||
        throw(ArgumentError("surrounded endpoint distance must be positive"))
    tries > 0 || throw(ArgumentError("endpoint sampling attempts must be positive"))
    return RPOHyPRRLEndpointSampler(
        base, triangles, _rpo_hypr_rl_triangle_cdf(triangles),
        minimum_clearance, maximum_clearance, minimum_separation,
        surrounded_distance, tries,
    )
end

function _rpo_hypr_rl_sample_triangle_point(rng::AbstractRNG, v1, v2, v3)
    u = rand(rng)
    v = rand(rng)
    if u + v > 1.0
        u = 1.0 - u
        v = 1.0 - v
    end
    return v1 .+ u .* (v2 .- v1) .+ v .* (v3 .- v1)
end

function _rpo_hypr_rl_directional_nearest_distances(point, station_points)
    minimum_distance_squared = fill(Inf, 6)
    for column in axes(station_points, 2)
        dx = station_points[1, column] - point[1]
        dy = station_points[2, column] - point[2]
        dz = station_points[3, column] - point[3]
        distance_squared = dx * dx + dy * dy + dz * dz
        dx > 0.0 && (minimum_distance_squared[1] = min(
            minimum_distance_squared[1], distance_squared,
        ))
        dx < 0.0 && (minimum_distance_squared[2] = min(
            minimum_distance_squared[2], distance_squared,
        ))
        dy > 0.0 && (minimum_distance_squared[3] = min(
            minimum_distance_squared[3], distance_squared,
        ))
        dy < 0.0 && (minimum_distance_squared[4] = min(
            minimum_distance_squared[4], distance_squared,
        ))
        dz > 0.0 && (minimum_distance_squared[5] = min(
            minimum_distance_squared[5], distance_squared,
        ))
        dz < 0.0 && (minimum_distance_squared[6] = min(
            minimum_distance_squared[6], distance_squared,
        ))
    end
    return sqrt.(minimum_distance_squared)
end

function _rpo_hypr_rl_is_surrounded(point, station_points,
                                    maximum_distance_m::Real)
    distances = _rpo_hypr_rl_directional_nearest_distances(
        point, station_points,
    )
    return all(distance -> distance <= maximum_distance_m, distances)
end

function _rpo_hypr_rl_near_surface_endpoint(
    sampler::RPOHyPRRLEndpointSampler, rng::AbstractRNG,
)
    geometry = sampler.base.geometry
    station_points = getproperty(getproperty(geometry, :station), :points_body)
    half_extents = getproperty(getproperty(geometry, :chaser), :half_extents_body)
    buffer = getproperty(getproperty(geometry, :station), :keepout_radius_m) +
        Float64(Base.invokelatest(maximum, half_extents))
    mesh_center = vec(mean(sampler.station_triangles; dims=2))
    n_triangles = size(sampler.station_triangles, 2) ÷ 3
    for _ in 1:sampler.max_sampling_tries
        triangle_index = clamp(
            searchsortedfirst(sampler.triangle_cdf, rand(rng)), 1, n_triangles,
        )
        base = 3 * (triangle_index - 1)
        v1 = view(sampler.station_triangles, :, base + 1)
        v2 = view(sampler.station_triangles, :, base + 2)
        v3 = view(sampler.station_triangles, :, base + 3)
        surface_point = _rpo_hypr_rl_sample_triangle_point(rng, v1, v2, v3)
        normal = cross(v2 - v1, v3 - v1)
        normal_norm = norm(normal)
        normal_norm < 1.0e-9 && continue
        normal ./= normal_norm
        triangle_center = (v1 .+ v2 .+ v3) ./ 3.0
        dot(normal, triangle_center - mesh_center) < 0.0 && (normal .*= -1.0)
        clearance = sampler.endpoint_min_clearance_m + rand(rng) *
            (sampler.endpoint_max_clearance_m - sampler.endpoint_min_clearance_m)
        candidate = surface_point .+ normal .* (buffer + clearance)
        _rpo_endpoint_clearance_m(candidate, geometry) + 1.0e-9 >=
            sampler.endpoint_min_clearance_m && return candidate
    end
    throw(ErrorException(
        "failed to sample a near-surface RPO endpoint after " *
        "$(sampler.max_sampling_tries) attempts",
    ))
end

function sample_rpo_hypr_rl_scenario(sampler::RPOHyPRRLEndpointSampler,
                                     rng::AbstractRNG)
    station_points = getproperty(
        getproperty(sampler.base.geometry, :station), :points_body,
    )
    start = nothing
    for _ in 1:sampler.max_sampling_tries
        candidate = _rpo_hypr_rl_near_surface_endpoint(sampler, rng)
        if !_rpo_hypr_rl_is_surrounded(
            candidate, station_points, sampler.surrounded_max_distance_m,
        )
            start = candidate
            break
        end
    end
    start === nothing && throw(ErrorException(
        "failed to sample an RPO start endpoint that is not surrounded",
    ))

    goal = nothing
    for _ in 1:sampler.max_sampling_tries
        candidate = _rpo_hypr_rl_near_surface_endpoint(sampler, rng)
        if norm(candidate - start) >= sampler.min_separation_m &&
           !_rpo_hypr_rl_is_surrounded(
               candidate, station_points, sampler.surrounded_max_distance_m,
           )
            goal = candidate
            break
        end
    end
    goal === nothing && throw(ErrorException(
        "failed to sample an RPO goal endpoint with the required separation",
    ))
    base = sampler.base
    return RPOHyPRRLScenario(
        start_rtn=Vector{Float64}(start), goal_rtn=Vector{Float64}(goal),
        geometry=base.geometry, pso_config=base.pso_config,
        tracking_settings=base.tracking_settings, rrt_settings=base.rrt_settings,
        initial_attitude_rtn_to_body=base.initial_attitude_rtn_to_body,
        final_attitude_rtn_to_body=base.final_attitude_rtn_to_body,
    )
end

(sampler::RPOHyPRRLEndpointSampler)(rng::AbstractRNG, _) =
    sample_rpo_hypr_rl_scenario(sampler, rng)

"""Sample training/evaluation endpoints from the configured feasible scenario distribution."""
function sample_rpo_hypr_rl_scenario(base::RPOHyPRRLScenario, sigma_m::Real,
                                     rng::AbstractRNG;
                                     safe_distance_m::Real,
                                     max_sampling_tries::Int=10_000)
    sigma = Float64(sigma_m)
    safe_distance = Float64(safe_distance_m)
    for _ in 1:max_sampling_tries
        start = base.start_rtn .+ sigma .* randn(rng, 3)
        goal = base.goal_rtn .+ sigma .* randn(rng, 3)
        start_clearance = _rpo_endpoint_clearance_m(start, base.geometry)
        goal_clearance = _rpo_endpoint_clearance_m(goal, base.geometry)
        if _rpo_clearance_is_feasible(
               start_clearance, safe_distance, base.pso_config,
           ) && _rpo_clearance_is_feasible(
               goal_clearance, safe_distance, base.pso_config,
           )
            return RPOHyPRRLScenario(
                start_rtn=start,
                goal_rtn=goal,
                geometry=base.geometry,
                pso_config=base.pso_config,
                tracking_settings=base.tracking_settings,
                rrt_settings=base.rrt_settings,
                initial_attitude_rtn_to_body=base.initial_attitude_rtn_to_body,
                final_attitude_rtn_to_body=base.final_attitude_rtn_to_body,
            )
        end
    end
    throw(ErrorException(
        "failed to sample feasible RPO start and goal endpoints after " *
        "$max_sampling_tries attempts",
    ))
end

function _limit_seed_waypoints(path::AbstractMatrix{<:Real}, maximum_internal::Int)
    n = size(path, 2)
    n <= maximum_internal + 2 && return Matrix{Float64}(path)
    selected = unique(round.(Int, range(1, n; length=maximum_internal + 2)))
    return Matrix{Float64}(path[:, selected])
end

"""Construct the raw RRT topology seed edited and smoothed by HyPR-RL."""
function rpo_hypr_rl_seed_path(scenario::RPOHyPRRLScenario,
                               config::RPOHyPRRLConfig;
                               rng::AbstractRNG=Random.default_rng())
    scenario.seed_path_rtn === nothing ||
        return _limit_seed_waypoints(scenario.seed_path_rtn,
                                     config.max_translation_waypoints)
    modules, pso_config, _, rrt = _rpo_spaceagora_settings(scenario, config)
    for (label, point) in (("start", scenario.start_rtn), ("goal", scenario.goal_rtn))
        clearance = _rpo_endpoint_clearance_m(point, scenario.geometry)
        _rpo_clearance_is_feasible(
            clearance, config.safe_distance_m, pso_config,
        ) || throw(ArgumentError(
            "RPO $label endpoint clearance $clearance m is below the " *
            "required $(config.safe_distance_m) m",
        ))
    end

    fallback_path = hcat(scenario.start_rtn, scenario.goal_rtn)
    for _ in 1:3
        plan = Base.invokelatest(
            getproperty(modules.guidance, :rpo_rrt_connect_plan_path),
            scenario.start_rtn, scenario.goal_rtn, scenario.geometry, pso_config;
            safe_distance_m=config.safe_distance_m, settings=rrt, rng=rng,
        )
        fallback_path = plan.raw_path
        plan.path_found && return _limit_seed_waypoints(
            plan.raw_path, config.max_translation_waypoints,
        )
    end
    return _limit_seed_waypoints(
        fallback_path, config.max_translation_waypoints,
    )
end

struct _RPOThrusterAllocator
    directions::Matrix{Float64}
    maximum_thrust_n::Vector{Float64}
    force_to_pair_command::Matrix{Float64}
    fallback_projection::Matrix{Float64}
    impulse_to_propellant::Vector{Float64}
    wheel_energy_reference_j::Float64
    paired_axes::Bool
end

"""Per-thruster pulse-residual state carried between MPC updates."""
mutable struct _RPOThrusterPulseScheduler
    residual_impulse_ns::Vector{Float64}
    continuously_on::BitVector
end

function _rpo_thruster_pulse_scheduler(config::RPOHyPRRLConfig)
    n_thrusters = length(config.thruster_max_thrust_n)
    return _RPOThrusterPulseScheduler(
        zeros(n_thrusters), falses(n_thrusters),
    )
end

"""Delivered impulse for a completed linearly opening and closing valve pulse."""
@inline function _rpo_completed_pulse_impulse(maximum_thrust::Float64,
                                               commanded_on_time_s::Float64,
                                               config::RPOHyPRRLConfig)
    return maximum_thrust * (
        commanded_on_time_s - 0.5 * config.thruster_opening_time_s +
        0.5 * config.thruster_closing_time_s
    )
end

"""Commanded valve-on duration that realizes a completed requested impulse."""
@inline function _rpo_pulse_on_time_for_impulse(requested_impulse_ns::Float64,
                                                 maximum_thrust::Float64,
                                                 config::RPOHyPRRLConfig)
    return requested_impulse_ns / maximum_thrust +
        0.5 * (
            config.thruster_opening_time_s - config.thruster_closing_time_s
        )
end

function _rpo_thruster_allocator(config::RPOHyPRRLConfig)
    directions = config.thruster_directions_body
    axes = directions[:, [1, 3, 5]]
    paired_axes = all(isapprox(
        directions[:, 2 * axis], -directions[:, 2 * axis - 1];
        atol=1.0e-12, rtol=1.0e-12,
    ) for axis in 1:3) && isapprox(
        transpose(axes) * axes, Matrix{Float64}(I, 3, 3);
        atol=1.0e-12, rtol=1.0e-12,
    )
    pair_projection = paired_axes ? inv(axes) : zeros(3, 3)
    gram = directions * transpose(directions) + 1.0e-12I
    fallback_projection = transpose(directions) / gram
    impulse_to_propellant = 1.0 ./
        (config.thruster_isp_s .* config.g0_mps2)
    wheel_energy_reference = sum(
        config.reaction_wheel_max_momentum_nms^2 / (2.0 * inertia)
        for inertia in config.reaction_wheel_inertia_kgm2
    )
    return _RPOThrusterAllocator(
        directions,
        config.thruster_max_thrust_n,
        pair_projection,
        fallback_projection,
        impulse_to_propellant,
        wheel_energy_reference,
        paired_axes,
    )
end

"""Conservative acceleration available in every direction for one complete pulse."""
function _rpo_guaranteed_acceleration_limit_mps2(
    config::RPOHyPRRLConfig,
    dt_s::Real,
    allocator::_RPOThrusterAllocator=_rpo_thruster_allocator(config),
)
    allocator.paired_axes || throw(ArgumentError(
        "hard retiming acceleration constraints require opposing orthonormal thruster pairs",
    ))
    dt = Float64(dt_s)
    dt > 0.0 || throw(ArgumentError("retiming time step must be positive"))
    valve_fraction = clamp(
        1.0 - (config.thruster_opening_time_s +
               config.thruster_closing_time_s) / (2.0 * dt),
        0.0,
        1.0,
    )
    paired_force_limit = minimum(
        min(
            allocator.maximum_thrust_n[2 * axis - 1],
            allocator.maximum_thrust_n[2 * axis],
        ) for axis in 1:3
    )
    limit = valve_fraction * paired_force_limit / config.mass_kg
    limit > 0.0 || throw(ArgumentError(
        "thruster and valve settings provide no guaranteed acceleration",
    ))
    return limit
end

function _allocate_thrusters!(thrust::Vector{Float64},
                              desired_force_body::AbstractVector{<:Real},
    allocator::_RPOThrusterAllocator)
    if allocator.paired_axes
        @inbounds for axis in 1:3
            positive_index = 2 * axis - 1
            negative_index = 2 * axis
            value = allocator.force_to_pair_command[axis, 1] * desired_force_body[1] +
                allocator.force_to_pair_command[axis, 2] * desired_force_body[2] +
                allocator.force_to_pair_command[axis, 3] * desired_force_body[3]
            thrust[positive_index] = clamp(
                max(value, 0.0), 0.0, allocator.maximum_thrust_n[positive_index],
            )
            thrust[negative_index] = clamp(
                max(-value, 0.0), 0.0, allocator.maximum_thrust_n[negative_index],
            )
        end
    else
        mul!(thrust, allocator.fallback_projection, desired_force_body)
        @inbounds for index in eachindex(thrust)
            thrust[index] = clamp(
                thrust[index], 0.0, allocator.maximum_thrust_n[index],
            )
        end
    end
    return thrust
end

function _allocate_thrusters(desired_force_body::AbstractVector{<:Real},
                             config::RPOHyPRRLConfig)
    allocator = _rpo_thruster_allocator(config)
    return _allocate_thrusters!(zeros(6), desired_force_body, allocator)
end

function _thruster_translation_step(command_rtn, attitude_rtn_to_body, dt_s,
                                    config::RPOHyPRRLConfig,
                                    allocator::_RPOThrusterAllocator,
                                    scheduler::_RPOThrusterPulseScheduler)
    dt = Float64(dt_s)
    rotation_rtn_to_body = _quat_rotation_matrix(attitude_rtn_to_body)
    desired_acceleration = Vector{Float64}(command_rtn)
    desired_force_body = config.mass_kg .* rotation_rtn_to_body * desired_acceleration
    commanded_thrust = _allocate_thrusters!(
        zeros(6), desired_force_body, allocator,
    )
    allocated_force_body = allocator.directions * commanded_thrust
    hard_saturated = norm(allocated_force_body - desired_force_body) /
        config.mass_kg > 1.0e-9
    impulse = zeros(6)
    for thruster_index in 1:6
        maximum_thrust = config.thruster_max_thrust_n[thruster_index]
        if maximum_thrust <= 0.0
            continue
        end
        requested_impulse = scheduler.residual_impulse_ns[thruster_index] +
            commanded_thrust[thruster_index] * dt
        maximum_impulse = maximum_thrust * dt
        minimum_on_time = min(config.thruster_min_firing_time_s, dt)
        minimum_impulse = _rpo_completed_pulse_impulse(
            maximum_thrust, minimum_on_time, config,
        )

        if scheduler.continuously_on[thruster_index]
            # A valve that was already open either holds nominal thrust or
            # closes during this update.  The latter includes its finite
            # closing-tail impulse; no off-time dwell is introduced.
            maximum_closing_impulse = maximum_thrust * (
                dt - 0.5 * config.thruster_closing_time_s
            )
            if requested_impulse > maximum_closing_impulse + 1.0e-12
                impulse[thruster_index] = maximum_impulse
                scheduler.residual_impulse_ns[thruster_index] =
                    requested_impulse - maximum_impulse
                continue
            end
            held_on_time = clamp(
                requested_impulse / maximum_thrust -
                0.5 * config.thruster_closing_time_s,
                0.0, dt - config.thruster_closing_time_s,
            )
            impulse[thruster_index] = maximum_thrust * (
                held_on_time + 0.5 * config.thruster_closing_time_s
            )
            scheduler.residual_impulse_ns[thruster_index] =
                requested_impulse - impulse[thruster_index]
            scheduler.continuously_on[thruster_index] = false
            continue
        end

        available_on_time = max(0.0, dt - config.thruster_closing_time_s)
        maximum_completed_impulse = _rpo_completed_pulse_impulse(
            maximum_thrust, available_on_time, config,
        )
        if requested_impulse > maximum_completed_impulse + 1.0e-12
            # Leave the valve open at the update boundary rather than
            # reserving a closing tail in every saturated update.
            impulse[thruster_index] = maximum_thrust * (
                dt - 0.5 * config.thruster_opening_time_s
            )
            scheduler.residual_impulse_ns[thruster_index] =
                requested_impulse - impulse[thruster_index]
            scheduler.continuously_on[thruster_index] = true
            continue
        end
        requested_on_time = _rpo_pulse_on_time_for_impulse(
            requested_impulse, maximum_thrust, config,
        )
        can_fire = requested_impulse + 1.0e-12 >= minimum_impulse &&
                   requested_on_time + 1.0e-12 >= minimum_on_time &&
                   available_on_time + 1.0e-12 >= minimum_on_time
        if can_fire
            on_time = min(requested_on_time, available_on_time)
            impulse[thruster_index] = _rpo_completed_pulse_impulse(
                maximum_thrust, on_time, config,
            )
            scheduler.residual_impulse_ns[thruster_index] =
                requested_impulse - impulse[thruster_index]
        else
            scheduler.residual_impulse_ns[thruster_index] = requested_impulse
        end
    end
    average_thrust = impulse ./ dt
    realized_force_body = allocator.directions * average_thrust
    realized_acceleration = transpose(rotation_rtn_to_body) * realized_force_body /
                            config.mass_kg
    propellant = dot(impulse, allocator.impulse_to_propellant)
    error = norm(realized_acceleration - desired_acceleration)
    return (
        realized_acceleration_rtn=realized_acceleration,
        propellant_used_kg=propellant,
        thruster_impulse_ns=impulse,
        allocation_error_mps2=error,
        saturated=hard_saturated,
    )
end

function _thruster_translation_step(command_rtn, attitude_rtn_to_body, dt_s,
                                    config::RPOHyPRRLConfig,
                                    allocator::_RPOThrusterAllocator)
    return _thruster_translation_step(
        command_rtn, attitude_rtn_to_body, dt_s, config, allocator,
        _rpo_thruster_pulse_scheduler(config),
    )
end

function _thruster_translation_step(command_rtn, attitude_rtn_to_body, dt_s,
                                    config::RPOHyPRRLConfig)
    return _thruster_translation_step(
        command_rtn, attitude_rtn_to_body, dt_s, config,
        _rpo_thruster_allocator(config),
    )
end

function _account_thruster_and_attitude_history(u_hist::AbstractMatrix{<:Real},
                                                q_reference::AbstractMatrix{<:Real},
                                                dt_s::Real,
                                                scenario::RPOHyPRRLScenario,
                                                config::RPOHyPRRLConfig)
    n_steps = size(u_hist, 2)
    impulse = zeros(6)
    propellant = 0.0
    allocation_error = 0.0
    saturated_steps = 0
    wheel_energy = 0.0
    wheel_momentum = zeros(3)
    wheel_peak = 0.0
    attitude = _quat_normalize(scenario.initial_attitude_rtn_to_body)
    angular_rate = zeros(3)
    dt = Float64(dt_s)
    allocator = _rpo_thruster_allocator(config)
    scheduler = _rpo_thruster_pulse_scheduler(config)
    for step_index in 1:n_steps
        desired_attitude = view(q_reference, :, min(step_index, size(q_reference, 2)))
        attitude, angular_rate, torque = _integrate_attitude(
            attitude, angular_rate, desired_attitude, dt, config,
        )
        wheel_momentum .-= torque .* dt
        wheel_peak = max(wheel_peak, maximum(abs, wheel_momentum))
        wheel_energy += abs(dot(torque, angular_rate)) * dt

        desired_acceleration_rtn = Vector{Float64}(view(u_hist, :, step_index))
        actuation = _thruster_translation_step(
            desired_acceleration_rtn, attitude, dt, config, allocator, scheduler,
        )
        impulse .+= actuation.thruster_impulse_ns
        propellant += actuation.propellant_used_kg
        allocation_error += actuation.allocation_error_mps2 * dt
        saturated_steps += actuation.saturated ? 1 : 0
    end
    return (
        propellant_used_kg=propellant,
        allocation_error_impulse_mps=allocation_error,
        saturation_fraction=n_steps > 0 ? saturated_steps / n_steps : 0.0,
        thruster_impulse_ns=impulse,
        wheel_energy_j=wheel_energy,
        wheel_peak_momentum_nms=wheel_peak,
        final_attitude_rtn_to_body=attitude,
        final_angular_rate_radps=angular_rate,
    )
end

"""
Integrate ideal continuous-thrust fuel from a geometry-retimed HCW reference.

The geometric cost is retained from baseline HyPR, while the proxy-fuel term
uses the actual retimed positions, velocities, duration, and control timestep.
Endpoint velocities are zero, so the first and last intervals include departure
acceleration and terminal braking, respectively.
"""
function _stream_retimed_timing_aware_proxy_objective(
    profile,
    mean_motion_radps::Real,
    mass_kg::Real,
    isp_s::Real,
    g0_mps2::Real,
    fuel_reference_kg::Real,
    fuel_weight::Real,
    geometric_cost::Real,
    cost_cutoff::Float64,
    interpolate_along_path!,
)
    n_reference = length(profile.s_ref)
    n_steps = max(n_reference - 1, 0)
    dt = Float64(profile.dt_s)
    mean_motion = Float64(mean_motion_radps)
    fuel_scale = Float64(mass_kg) * dt /
        max(Float64(isp_s) * Float64(g0_mps2), 1.0e-12)
    fuel_reference = max(Float64(fuel_reference_kg), 1.0e-12)
    local_fuel_weight = Float64(fuel_weight)
    base_cost = Float64(geometric_cost)

    current_position = zeros(3)
    next_position = zeros(3)
    following_position = zeros(3)
    current_velocity = zeros(3)
    next_velocity = zeros(3)
    propellant = 0.0

    if n_steps > 0
        interpolate_along_path!(
            current_position, profile.samples, profile.s_samples, profile.s_ref[1],
        )
        interpolate_along_path!(
            next_position, profile.samples, profile.s_samples, profile.s_ref[2],
        )
        # Initial relative velocity is zero, matching the materialized reference.
    end

    cutoff_exceeded = false
    completed_steps = 0
    for step_index in 1:n_steps
        if step_index < n_steps
            interpolate_along_path!(
                following_position,
                profile.samples,
                profile.s_samples,
                profile.s_ref[step_index + 2],
            )
            @inbounds for axis in 1:3
                next_velocity[axis] =
                    (following_position[axis] - next_position[axis]) / dt
            end
        else
            next_velocity .= 0.0
        end

        @inbounds begin
            ax = (next_velocity[1] - current_velocity[1]) / dt -
                (3.0 * mean_motion^2 * current_position[1] +
                 2.0 * mean_motion * current_velocity[2])
            ay = (next_velocity[2] - current_velocity[2]) / dt +
                2.0 * mean_motion * current_velocity[1]
            az = (next_velocity[3] - current_velocity[3]) / dt +
                mean_motion^2 * current_position[3]
            propellant += fuel_scale * sqrt(ax * ax + ay * ay + az * az)
        end
        completed_steps = step_index

        normalized_fuel = propellant / fuel_reference
        partial_cost = base_cost + local_fuel_weight * normalized_fuel
        if isfinite(cost_cutoff) && partial_cost > cost_cutoff
            cutoff_exceeded = true
            break
        end

        if step_index < n_steps
            current_position .= next_position
            next_position .= following_position
            current_velocity .= next_velocity
        end
    end

    normalized_fuel = propellant / fuel_reference
    return (
        total=base_cost + local_fuel_weight * normalized_fuel,
        normalized_fuel=normalized_fuel,
        propellant_used_kg=propellant,
        duration_s=completed_steps * dt,
        cutoff_exceeded=cutoff_exceeded,
    )
end


function _stream_retimed_feedforward_objective(
    profile,
    tracking,
    initial_attitude_rtn_to_body::Vector{Float64},
    final_attitude_rtn_to_body::Vector{Float64},
    config::RPOHyPRRLConfig,
    allocator::_RPOThrusterAllocator,
    fuel_reference_kg::Real,
    cost_cutoff::Float64,
    interpolate_along_path!,
    ;
    optimize_pointing::Bool=false,
)
    n_reference = length(profile.s_ref)
    n_steps = max(n_reference - 1, 0)
    dt = Float64(profile.dt_s)
    mean_motion = Float64(tracking.mean_motion_radps)
    fuel_reference = max(Float64(fuel_reference_kg), 1.0e-12)
    wheel_reference = max(allocator.wheel_energy_reference_j, 1.0e-12)

    current_position = zeros(3)
    next_position = zeros(3)
    following_position = zeros(3)
    current_velocity = zeros(3)
    next_velocity = zeros(3)
    command = zeros(3)
    impulse = zeros(6)
    attitude = _quat_normalize(initial_attitude_rtn_to_body)
    angular_rate = zeros(3)
    wheel_momentum = zeros(3)
    propellant = 0.0
    allocation_error = 0.0
    saturated_steps = 0
    wheel_energy = 0.0
    wheel_peak = 0.0
    scheduler = _rpo_thruster_pulse_scheduler(config)

    if n_steps > 0
        interpolate_along_path!(
            current_position, profile.samples, profile.s_samples, profile.s_ref[1],
        )
        interpolate_along_path!(
            next_position, profile.samples, profile.s_samples, profile.s_ref[2],
        )
        # Initial relative velocity is zero, matching the materialized reference.
    end

    cutoff_exceeded = false
    completed_steps = 0
    for step_index in 1:n_steps
        if step_index < n_steps
            interpolate_along_path!(
                following_position,
                profile.samples,
                profile.s_samples,
                profile.s_ref[step_index + 2],
            )
            @inbounds for axis in 1:3
                next_velocity[axis] =
                    (following_position[axis] - next_position[axis]) / dt
            end
        else
            next_velocity .= 0.0
        end

        @inbounds begin
            command[1] = (next_velocity[1] - current_velocity[1]) / dt -
                (3.0 * mean_motion^2 * current_position[1] +
                 2.0 * mean_motion * current_velocity[2])
            command[2] = (next_velocity[2] - current_velocity[2]) / dt +
                2.0 * mean_motion * current_velocity[1]
            command[3] = (next_velocity[3] - current_velocity[3]) / dt +
                mean_motion^2 * current_position[3]
        end

        attitude_denominator = optimize_pointing ? n_steps - 1 : n_steps
        attitude_fraction = (step_index - 1) / max(attitude_denominator, 1)
        baseline_attitude = _quat_slerp(
            initial_attitude_rtn_to_body,
            final_attitude_rtn_to_body,
            attitude_fraction,
        )
        desired_attitude = optimize_pointing ?
            _rpo_fuel_optimal_pointing_attitude(
                attitude,
                baseline_attitude,
                command,
                attitude_fraction,
                config,
            ) : baseline_attitude
        attitude, angular_rate, torque = _integrate_attitude(
            attitude, angular_rate, desired_attitude, dt, config,
        )
        wheel_momentum .-= torque .* dt
        wheel_peak = max(wheel_peak, maximum(abs, wheel_momentum))
        wheel_energy += abs(dot(torque, angular_rate)) * dt

        actuation = _thruster_translation_step(
            command, attitude, dt, config, allocator, scheduler,
        )
        impulse .+= actuation.thruster_impulse_ns
        propellant += actuation.propellant_used_kg
        allocation_error += actuation.allocation_error_mps2 * dt
        saturated_steps += actuation.saturated ? 1 : 0
        completed_steps = step_index

        normalized_fuel = propellant / fuel_reference
        normalized_wheel_energy = wheel_energy / wheel_reference
        partial_cost = config.fuel_weight * normalized_fuel +
            config.wheel_weight * normalized_wheel_energy
        if isfinite(cost_cutoff) && partial_cost > cost_cutoff
            cutoff_exceeded = true
            break
        end

        if step_index < n_steps
            current_position .= next_position
            next_position .= following_position
            current_velocity .= next_velocity
        end
    end

    normalized_fuel = propellant / fuel_reference
    normalized_wheel_energy = wheel_energy / wheel_reference
    total = config.fuel_weight * normalized_fuel +
        config.wheel_weight * normalized_wheel_energy
    return (
        total=total,
        normalized_fuel=normalized_fuel,
        normalized_wheel_energy=normalized_wheel_energy,
        propellant_used_kg=propellant,
        allocation_error_impulse_mps=allocation_error,
        saturation_fraction=completed_steps > 0 ?
            saturated_steps / completed_steps : 0.0,
        thruster_impulse_ns=impulse,
        wheel_energy_j=wheel_energy,
        wheel_peak_momentum_nms=wheel_peak,
        duration_s=completed_steps * dt,
        cutoff_exceeded=cutoff_exceeded,
    )
end

function _rpo_feedforward_acceleration_history(r_ref::AbstractMatrix{<:Real},
                                               v_ref::AbstractMatrix{<:Real},
                                               mean_motion_radps::Real,
                                               dt_s::Real)
    size(r_ref, 1) == 3 || throw(DimensionMismatch("r_ref must have three rows"))
    size(v_ref) == size(r_ref) ||
        throw(DimensionMismatch("v_ref must match r_ref"))
    n_steps = max(size(r_ref, 2) - 1, 0)
    commands = zeros(3, n_steps)
    n = Float64(mean_motion_radps)
    dt = Float64(dt_s)
    dt > 0.0 || throw(ArgumentError("feedforward time step must be positive"))
    for step_index in 1:n_steps
        position = view(r_ref, :, step_index)
        velocity = view(v_ref, :, step_index)
        acceleration = (view(v_ref, :, step_index + 1) .- velocity) ./ dt
        drift = (
            3.0 * n^2 * position[1] + 2.0 * n * velocity[2],
            -2.0 * n * velocity[1],
            -n^2 * position[3],
        )
        commands[:, step_index] .= acceleration .- drift
    end
    return commands
end

function _rpo_infeasible_path_evaluation(
    scenario::RPOHyPRRLScenario,
    control_points_rtn::AbstractMatrix{<:Real},
    components;
    reason::Symbol=:path_infeasible,
)
    clearance_penalty = _rpo_original_clearance_penalty(components)
    final_error = norm(
        view(control_points_rtn, :, size(control_points_rtn, 2)) -
        scenario.goal_rtn,
    )
    return RPOHyPRRLEvaluation(
        false, clearance_penalty, clearance_penalty,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        Float64(components.min_clearance), final_error,
        Float64[], zeros(3, 0), zeros(3, 0), zeros(4, 0), zeros(6),
        (
            reason=reason,
            evaluator_mode=:retimed_feedforward,
            path_components=components,
            clearance_penalty=clearance_penalty,
        ),
    )
end

"""
Evaluate an RPO editor candidate using the retimed reference and nominal HCW
feedforward acceleration. This is the inexpensive training fidelity: it keeps
six-thruster pulse/Isp accounting and reaction-wheel attitude propagation but
does not run the closed-loop LQ-MPC tracker.
"""
function evaluate_rpo_training_candidate(
    scenario::RPOHyPRRLScenario,
    config::RPOHyPRRLConfig,
    control_points_rtn::AbstractMatrix{<:Real},
    attitude_progress::AbstractVector{<:Real},
    attitude_quaternions::AbstractMatrix{<:Real},
)
    try
        modules, pso_config, tracking, _ = _rpo_spaceagora_settings(scenario, config)
        prepared = Base.invokelatest(
            getproperty(modules.guidance, :rpo_prepare_retimed_candidate),
            control_points_rtn, scenario.geometry, pso_config;
            safe_distance_m=config.safe_distance_m,
            retime_dt_s=tracking.dt_s,
            w_len=0.0,
            w_fuel=0.0,
            compute_fuel_proxy=false,
            retime_mean_motion_radps=tracking.mean_motion_radps,
            retime_command_limit_mps2=
                _rpo_guaranteed_acceleration_limit_mps2(config, tracking.dt_s),
        )
        components = prepared.components
        if !isfinite(components.total) || !isfinite(components.min_clearance)
            return _failed_rpo_evaluation(reason=:nonfinite_path_cost)
        end
        if components.violation_count > 0 || !_rpo_clearance_is_feasible(
               components.min_clearance, config.safe_distance_m, pso_config,
           )
            return _rpo_infeasible_path_evaluation(
                scenario, control_points_rtn, components,
            )
        end
        prepared.profile === nothing &&
            return _failed_rpo_evaluation(reason=:missing_retiming_profile)
        t_ref, r_ref, v_ref = Base.invokelatest(
            getproperty(modules.guidance, :rpo_reference_from_profile),
            prepared.profile,
        )
        commands = _rpo_feedforward_acceleration_history(
            r_ref, v_ref, tracking.mean_motion_radps, tracking.dt_s,
        )
        reference_progress = isempty(t_ref) ? Float64[] :
            collect(range(0.0, 1.0; length=length(t_ref)))
        q_reference = _attitude_reference(
            attitude_progress, attitude_quaternions, reference_progress,
        )
        accounting = _account_thruster_and_attitude_history(
            commands, q_reference, tracking.dt_s, scenario, config,
        )
        final_error = isempty(r_ref) ? Inf :
            norm(view(r_ref, :, size(r_ref, 2)) - scenario.goal_rtn)
        wheel_feasible = accounting.wheel_peak_momentum_nms <=
            config.reaction_wheel_max_momentum_nms + 1.0e-12
        feasible = final_error <= tracking.final_position_tol_m && wheel_feasible
        duration_s = isempty(t_ref) ? 0.0 : t_ref[end]
        objective = _rpo_fuel_wheel_objective(
            components, accounting.propellant_used_kg,
            accounting.wheel_energy_j, config,
        )
        return RPOHyPRRLEvaluation(
            feasible, objective.total, objective.total, accounting.propellant_used_kg,
            duration_s, accounting.allocation_error_impulse_mps,
            accounting.saturation_fraction, accounting.wheel_energy_j,
            accounting.wheel_peak_momentum_nms, components.min_clearance,
            final_error, Vector{Float64}(t_ref), Matrix{Float64}(r_ref),
            Matrix{Float64}(v_ref), q_reference, accounting.thruster_impulse_ns,
            (
                evaluator_mode=:retimed_feedforward,
                path_components=components,
                objective_components=objective,
                coupled_tracking_success=missing,
                keepout_violations=0,
                wheel_feasible=wheel_feasible,
                feedforward_command_history=commands,
                final_attitude_rtn_to_body=accounting.final_attitude_rtn_to_body,
                final_angular_rate_radps=accounting.final_angular_rate_radps,
            ),
        )
    catch error
        return _failed_rpo_evaluation(
            reason=:evaluation_exception,
            error_message=sprint(showerror, error),
        )
    end
end

function _coupled_rpo_tracking(modules, t_ref, r_ref, v_ref, goal_rtn,
                               attitude_progress, attitude_quaternions,
                               tracking, scenario::RPOHyPRRLScenario,
                               config::RPOHyPRRLConfig;
                               optimize_pointing::Bool=false)
    Q = Diagonal(vcat(fill(tracking.q_pos, 3), fill(tracking.q_vel, 3)))
    R = Diagonal(fill(tracking.r_accel, 3))
    Qf = Diagonal(vcat(fill(tracking.qf_pos, 3), fill(tracking.qf_vel, 3)))
    u_max = Vector{Float64}(tracking.u_max_mps2)
    controller = Base.invokelatest(
        getproperty(modules.control, :init_rpo_lqmpc),
        tracking.mean_motion_radps, tracking.dt_s, Q, R, Qf, tracking.horizon;
        u_min=-u_max, u_max=u_max,
    )
    n_plan_steps = max(size(r_ref, 2) - 1, 0)
    settle_steps = max(0, Int(ceil(tracking.settle_time_s / tracking.dt_s)))
    total_steps = n_plan_steps + settle_steps
    control_progress = total_steps == 0 ? Float64[] :
        collect(range(0.0, 1.0; length=total_steps))
    baseline_attitudes = _attitude_reference(
        attitude_progress, attitude_quaternions, control_progress,
    )
    desired_attitude_history = zeros(4, total_steps)

    state = vcat(Vector{Float64}(r_ref[:, 1]), Vector{Float64}(v_ref[:, 1]))
    plant = Base.invokelatest(
        getproperty(modules.guidance, :rpo_init_two_body_plant),
        state, tracking.mean_motion_radps,
    )
    state_history = zeros(6, total_steps + 1)
    command_history = zeros(3, total_steps)
    realized_history = zeros(3, total_steps)
    actual_attitude_history = zeros(4, total_steps + 1)
    state_history[:, 1] .= state
    attitude = _quat_normalize(scenario.initial_attitude_rtn_to_body)
    angular_rate = zeros(3)
    actual_attitude_history[:, 1] .= attitude
    wheel_momentum = zeros(3)
    wheel_peak = 0.0
    wheel_energy = 0.0
    impulse = zeros(6)
    propellant = 0.0
    allocation_error = 0.0
    saturated_steps = 0
    min_clearance = Inf
    keepout_violations = 0
    dt = Float64(tracking.dt_s)
    allocator = _rpo_thruster_allocator(config)
    scheduler = _rpo_thruster_pulse_scheduler(config)

    for step_index in 1:total_steps
        reference_index = min(step_index, size(r_ref, 2))
        preview = Base.invokelatest(
            getproperty(modules.guidance, :rpo_lqmpc_reference_preview),
            r_ref, v_ref, reference_index, tracking.horizon,
        )
        command = Vector{Float64}(Base.invokelatest(
            getproperty(modules.control, :rpo_lqmpc_control),
            controller, state, preview,
        ))
        baseline_attitude = view(baseline_attitudes, :, step_index)
        desired_attitude = optimize_pointing ?
            _rpo_fuel_optimal_pointing_attitude(
                attitude,
                baseline_attitude,
                command,
                control_progress[step_index],
                config,
            ) : baseline_attitude
        desired_attitude_history[:, step_index] .= desired_attitude
        attitude, angular_rate, torque = _integrate_attitude(
            attitude, angular_rate, desired_attitude, dt, config,
        )
        wheel_momentum .-= torque .* dt
        wheel_peak = max(wheel_peak, maximum(abs, wheel_momentum))
        wheel_energy += abs(dot(torque, angular_rate)) * dt

        actuation = _thruster_translation_step(
            command, attitude, dt, config, allocator, scheduler,
        )
        realized_acceleration = actuation.realized_acceleration_rtn
        impulse .+= actuation.thruster_impulse_ns
        propellant += actuation.propellant_used_kg
        allocation_error += actuation.allocation_error_mps2 * dt
        saturated_steps += actuation.saturated ? 1 : 0
        state .= Base.invokelatest(
            getproperty(modules.guidance, :rpo_step_two_body!),
            plant, realized_acceleration, dt,
        )
        clearance = Float64(Base.invokelatest(
            getproperty(modules.navigation, :rpo_capsule_clearance_to_station),
            view(state_history, 1:3, step_index), view(state, 1:3), scenario.geometry,
        ))
        min_clearance = min(min_clearance, clearance)
        keepout_violations += _rpo_clearance_is_feasible(
            clearance, config.safe_distance_m, scenario.pso_config,
        ) ? 0 : 1
        state_history[:, step_index + 1] .= state
        command_history[:, step_index] .= command
        realized_history[:, step_index] .= realized_acceleration
        actual_attitude_history[:, step_index + 1] .= attitude
    end
    final_error = norm(state[1:3] - scenario.goal_rtn)
    wheel_feasible = wheel_peak <= config.reaction_wheel_max_momentum_nms + 1.0e-12
    success = final_error <= tracking.final_position_tol_m &&
              keepout_violations == 0 && wheel_feasible
    return (
        success=success,
        final_position_error_m=final_error,
        min_clearance_m=min_clearance,
        keepout_violations=keepout_violations,
        propellant_used_kg=propellant,
        allocation_error_impulse_mps=allocation_error,
        saturation_fraction=total_steps > 0 ? saturated_steps / total_steps : 0.0,
        thruster_impulse_ns=impulse,
        wheel_energy_j=wheel_energy,
        wheel_peak_momentum_nms=wheel_peak,
        wheel_feasible=wheel_feasible,
        state_history=state_history,
        command_history=command_history,
        realized_acceleration_history=realized_history,
        actual_attitude_history=actual_attitude_history,
        desired_attitude_history=desired_attitude_history,
        final_attitude_rtn_to_body=attitude,
        final_angular_rate_radps=angular_rate,
    )
end

"""
Evaluate a translation/attitude candidate after path retiming and LQ-MPC
tracking with nonlinear two-body orbital propagation. Translation is allocated
through all six fixed body thrusters; the
attitude is tracked using reaction-wheel torque and momentum limits.
"""
function evaluate_rpo_candidate(scenario::RPOHyPRRLScenario,
                                config::RPOHyPRRLConfig,
                                control_points_rtn::AbstractMatrix{<:Real},
                                attitude_progress::AbstractVector{<:Real},
                                attitude_quaternions::AbstractMatrix{<:Real};
                                optimize_pointing::Bool=false)
    try
        modules, pso_config, tracking, _ = _rpo_spaceagora_settings(scenario, config)
        prepared = Base.invokelatest(
            getproperty(modules.guidance, :rpo_prepare_retimed_candidate),
            control_points_rtn, scenario.geometry, pso_config;
            safe_distance_m=config.safe_distance_m,
            retime_dt_s=tracking.dt_s,
            w_len=0.0,
            w_fuel=0.0,
            compute_fuel_proxy=false,
            retime_mean_motion_radps=tracking.mean_motion_radps,
            retime_command_limit_mps2=
                _rpo_guaranteed_acceleration_limit_mps2(config, tracking.dt_s),
        )
        components = prepared.components
        if components.violation_count > 0 || !isfinite(components.total)
            return _failed_rpo_evaluation(reason=:path_infeasible)
        end
        prepared.profile === nothing &&
            return _failed_rpo_evaluation(reason=:missing_retiming_profile)
        t_ref, r_ref, v_ref = Base.invokelatest(
            getproperty(modules.guidance, :rpo_reference_from_profile),
            prepared.profile,
        )
        coupled = _coupled_rpo_tracking(
            modules, t_ref, r_ref, v_ref, scenario.goal_rtn,
            attitude_progress, attitude_quaternions, tracking, scenario, config,
            optimize_pointing=optimize_pointing,
        )
        q_reference = if optimize_pointing && !isempty(t_ref)
            reference = zeros(4, length(t_ref))
            reference[:, 1] .= _quat_normalize(
                scenario.initial_attitude_rtn_to_body,
            )
            @inbounds for column in 2:length(t_ref)
                step_index = min(
                    column - 1, size(coupled.desired_attitude_history, 2),
                )
                reference[:, column] .= view(
                    coupled.desired_attitude_history, :, step_index,
                )
            end
            reference
        else
            reference_progress = isempty(t_ref) ? Float64[] :
                collect(range(0.0, 1.0; length=length(t_ref)))
            _attitude_reference(
                attitude_progress, attitude_quaternions, reference_progress,
            )
        end
        objective = _rpo_fuel_wheel_objective(
            components, coupled.propellant_used_kg,
            coupled.wheel_energy_j, config,
        )
        return RPOHyPRRLEvaluation(
            coupled.success, objective.total, objective.total, coupled.propellant_used_kg,
            isempty(t_ref) ? 0.0 : t_ref[end],
            coupled.allocation_error_impulse_mps,
            coupled.saturation_fraction,
            coupled.wheel_energy_j,
            coupled.wheel_peak_momentum_nms,
            coupled.min_clearance_m,
            coupled.final_position_error_m,
            Vector{Float64}(t_ref),
            Matrix{Float64}(r_ref),
            Matrix{Float64}(v_ref),
            q_reference,
            coupled.thruster_impulse_ns,
            (
                evaluator_mode=:full_lqmpc,
                propagation_model=:nonlinear_two_body,
                path_components=components,
                objective_components=objective,
                coupled_tracking_success=coupled.success,
                keepout_violations=coupled.keepout_violations,
                wheel_feasible=coupled.wheel_feasible,
                state_history=coupled.state_history,
                command_history=coupled.command_history,
                realized_acceleration_history=coupled.realized_acceleration_history,
                actual_attitude_history=coupled.actual_attitude_history,
                desired_attitude_history=coupled.desired_attitude_history,
                pointing_optimized=optimize_pointing,
                final_attitude_rtn_to_body=coupled.final_attitude_rtn_to_body,
                final_angular_rate_radps=coupled.final_angular_rate_radps,
            ),
        )
    catch error
        return _failed_rpo_evaluation(
            reason=:evaluation_exception,
            error_message=sprint(showerror, error),
        )
    end
end
