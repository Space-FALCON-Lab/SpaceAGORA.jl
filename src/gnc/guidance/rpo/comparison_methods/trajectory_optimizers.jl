const _RPO_PLOTLYJS_MODULE = Ref{Union{Nothing,Module}}(nothing)
const _RPO_PLOTLYJS_LOAD_LOCK = ReentrantLock()

"""
Load PlotlyJS only when an RPO plotting function is called.

This keeps the optional WebIO/PlotlyJS visualization stack out of headless
simulation processes, where another dependency may already have loaded a newer
and WebIO-incompatible JSON.jl version.
"""
function _rpo_plotlyjs()
    return lock(_RPO_PLOTLYJS_LOAD_LOCK) do
        if _RPO_PLOTLYJS_MODULE[] === nothing
            _RPO_PLOTLYJS_MODULE[] = Base.require(
                Base.PkgId(
                    Base.UUID("f0f68f2c-4968-5e81-91da-67840de0976a"),
                    "PlotlyJS",
                ),
            )
        end
        return _RPO_PLOTLYJS_MODULE[]::Module
    end
end

"""CHOMP-like trajectory optimizer settings for RPO comparison planning."""
Base.@kwdef struct RPOCHOMPSettings
    n_iters::Int = 100
    learning_rate::Float64 = 0.06
    gradient_eps::Float64 = 1.0e-3
    w_smooth::Float64 = 1.0
end

"""STOMP-like trajectory optimizer settings for RPO comparison planning."""
Base.@kwdef struct RPOSTOMPSettings
    n_iters::Int = 100
    n_rollouts::Int = 20
    noise_std::Float64 = 0.25
    lambda::Float64 = 10.0
    update_step::Float64 = 1.0
    w_smooth::Float64 = 1.0
end

"""Shared trajectory optimizer settings for CHOMP and STOMP comparisons."""
Base.@kwdef struct RPOTrajectoryOptimizerSettings
    w_obs_scale::Float64 = 10.0
    obstacle_margin_m::Float64 = 0.5
    no_change_iters::Int = 8
    no_change_tol::Float64 = 1.0e-7
    match_hypr_iters::Bool = false
    match_hypr_runtime::Bool = false
    runtime_limit_s::Float64 = 30.0
    runtime_max_iters::Int = 100_000
    force_full_iters::Bool = true
end

"""Normalize planner identifiers to canonical symbols."""
function normalize_rpo_comparison_planner_type(planner_type)
    planner = planner_type isa Symbol ?
        Symbol(replace(lowercase(string(planner_type)), "-" => "_")) :
        Symbol(replace(lowercase(String(planner_type)), "-" => "_"))
    planner in (:hypr, :pso) && return :hypr
    planner in (:pso_unrefined, :unrefined_pso, :hypr_unrefined, :unrefined_hypr) && return :pso_unrefined
    planner in (:rrt_connect, :rrtconnect) && return :rrt_connect
    planner in (:rrt_connect_bezier, :rrtconnect_bezier, :rrt_bezier, :rrt_connect_smooth, :rrt_connect_smoothed) && return :rrt_connect_bezier
    planner in (:rrt_star, :rrtstar) && return :rrt_star
    planner in (:chomp, :stomp) && return planner
    throw(ArgumentError("Unsupported RPO comparison planner $(planner_type). Expected :hypr, :pso_unrefined, :rrt_connect, :rrt_connect_bezier, :rrt_star, :chomp, or :stomp."))
end

"""Return a display label for an RPO comparison planner."""
function rpo_comparison_planner_label(planner_type)
    planner = normalize_rpo_comparison_planner_type(planner_type)
    planner == :hypr && return "HYPR"
    planner == :pso_unrefined && return "PSO (unrefined)"
    planner == :rrt_connect && return "RRT-Connect"
    planner == :rrt_connect_bezier && return "RRT-Connect + Bezier"
    planner == :rrt_star && return "RRT*"
    planner == :chomp && return "CHOMP"
    planner == :stomp && return "STOMP"
    return uppercase(string(planner))
end

"""Create linearly spaced internal waypoints between start and goal."""
function rpo_trajectory_internal_points(start, goal, n_waypoints::Integer)
    theta = zeros(3, n_waypoints)
    for i in 1:n_waypoints
        α = i / (n_waypoints + 1)
        theta[:, i] .= (1.0 - α) .* start .+ α .* goal
    end
    return theta
end

"""Resample an existing path into internal optimizer waypoints."""
function rpo_trajectory_internal_points_from_seed(seed_path, start, goal, n_waypoints::Integer)
    n_waypoints <= 0 && return zeros(3, 0)
    seed_path === nothing && return rpo_trajectory_internal_points(start, goal, n_waypoints)

    pts = Matrix{Float64}(seed_path)
    size(pts, 1) == 3 || return rpo_trajectory_internal_points(start, goal, n_waypoints)
    size(pts, 2) >= 2 || return rpo_trajectory_internal_points(start, goal, n_waypoints)
    pts[:, 1] .= start
    pts[:, end] .= goal

    cumulative = zeros(size(pts, 2))
    for i in 2:size(pts, 2)
        cumulative[i] = cumulative[i - 1] + norm(pts[:, i] - pts[:, i - 1])
    end
    total_len = cumulative[end]
    total_len > 1.0e-12 || return rpo_trajectory_internal_points(start, goal, n_waypoints)

    theta = zeros(3, n_waypoints)
    for i in 1:n_waypoints
        s = total_len * i / (n_waypoints + 1)
        seg = clamp(searchsortedfirst(cumulative, s), 2, length(cumulative))
        span = cumulative[seg] - cumulative[seg - 1]
        α = span <= 1.0e-12 ? 0.0 : (s - cumulative[seg - 1]) / span
        theta[:, i] .= (1.0 - α) .* pts[:, seg - 1] .+ α .* pts[:, seg]
    end
    return theta
end

"""Combine start, internal waypoints, and goal into a full trajectory matrix."""
function rpo_trajectory_points_from_internal(theta, start, goal)
    n_waypoints = size(theta, 2)
    points = zeros(3, n_waypoints + 2)
    points[:, 1] .= start
    n_waypoints > 0 && (points[:, 2:(end - 1)] .= theta)
    points[:, end] .= goal
    return points
end

"""Create optimizer bounds around the start-goal corridor."""
function rpo_trajectory_search_bounds(start, goal, search_margin, n_waypoints::Integer)
    margin = max(0.0, Float64(search_margin))
    lo = min.(start, goal) .- margin
    hi = max.(start, goal) .+ margin
    return repeat(lo, 1, n_waypoints), repeat(hi, 1, n_waypoints)
end

"""Clamp internal waypoints to optimizer search bounds."""
rpo_clamp_internal_waypoints(theta, lo, hi) = size(theta, 2) == 0 ? copy(theta) : clamp.(theta, lo, hi)

"""Build the finite-difference smoothness matrix for internal waypoints."""
function rpo_second_difference_metric(n_waypoints::Integer)
    n_waypoints <= 0 && return zeros(0, 0), zeros(0, 0), zeros(0, 0)

    A = zeros(n_waypoints, n_waypoints)
    for i in 1:n_waypoints
        A[i, i] = -2.0
        i > 1 && (A[i, i - 1] = 1.0)
        i < n_waypoints && (A[i, i + 1] = 1.0)
    end

    R = A' * A + 1.0e-8 .* Matrix{Float64}(I, n_waypoints, n_waypoints)
    R_inv = Matrix(inv(Symmetric(R)))
    M = copy(R_inv)
    for j in 1:n_waypoints
        max_col = maximum(abs.(M[:, j]))
        if max_col > 1.0e-12
            M[:, j] .*= (1.0 / n_waypoints) / max_col
        end
    end
    return R, R_inv, M
end

"""Measure squared second-difference smoothness for a waypoint trajectory."""
function rpo_trajectory_smoothness_cost(points)
    pts = Matrix{Float64}(points)
    size(pts, 2) < 3 && return 0.0
    acc = 0.0
    for i in 2:(size(pts, 2) - 1)
        d2 = pts[:, i - 1] .- 2.0 .* pts[:, i] .+ pts[:, i + 1]
        acc += dot(d2, d2)
    end
    return acc
end

"""Evaluate CHOMP-style soft obstacle potential from clearance."""
function rpo_chomp_obstacle_potential(clearance, safe_distance, margin)
    d = clearance - safe_distance
    eps = max(Float64(margin), 1.0e-6)
    if d < 0.0
        return -d + 0.5 * eps
    elseif d <= eps
        return 0.5 * (d - eps)^2 / eps
    end
    return 0.0
end

"""Accumulate soft obstacle cost over sampled trajectory points."""
function rpo_soft_obstacle_cost_from_samples(samples, geometry; safe_distance_m, obstacle_margin_m)
    n_samples = size(samples, 2)
    n_samples == 0 && return 0.0

    acc = 0.0
    for i in 1:n_samples
        clearance = rpo_clearance_distance_to_station(samples[:, i], geometry)
        acc += rpo_chomp_obstacle_potential(clearance, safe_distance_m, obstacle_margin_m)
    end
    return acc / n_samples
end

"""Compute the smoothness-plus-obstacle objective for trajectory optimizers."""
function rpo_trajectory_soft_objective(points, geometry, cfg::RPOPSOConfig; safe_distance_m, obstacle_margin_m, w_smooth)
    samples = rpo_sample_path(
        points,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=cfg.sample_ds_m,
        curve_type=cfg.curve_type,
    )
    refs = rpo_path_cost_normalization_refs(points, cfg)
    J_len = rpo_path_length(samples)
    J_fuel = cfg.w_fuel > 0.0 ? rpo_fuel_proxy_from_samples(samples, cfg) : 0.0
    J_obs_soft = rpo_soft_obstacle_cost_from_samples(
        samples,
        geometry;
        safe_distance_m=safe_distance_m,
        obstacle_margin_m=obstacle_margin_m,
    )
    J_smooth = rpo_trajectory_smoothness_cost(points) / max(refs.len_ref^2, 1.0e-9)
    return cfg.w_len * (J_len / refs.len_ref)^2 +
        cfg.w_fuel * (J_fuel / refs.fuel_ref)^2 +
        cfg.w_obs * J_obs_soft +
        Float64(w_smooth) * J_smooth
end

"""Estimate the CHOMP objective gradient with centered finite differences."""
function rpo_chomp_numeric_gradient(theta, start, goal, objective_fn, gradient_eps)
    grad = zeros(size(theta))
    for i in 1:size(theta, 2)
        for d in 1:3
            h = max(Float64(gradient_eps), 1.0e-6 * max(abs(theta[d, i]), 1.0))
            theta[d, i] += h
            c_plus = objective_fn(rpo_trajectory_points_from_internal(theta, start, goal))
            theta[d, i] -= 2.0 * h
            c_minus = objective_fn(rpo_trajectory_points_from_internal(theta, start, goal))
            theta[d, i] += h
            grad[d, i] = (c_plus - c_minus) / (2.0 * h)
        end
    end
    return grad
end

"""Plan an RPO path using the CHOMP-like comparison optimizer."""
function rpo_chomp_plan_path(
    start_rtn,
    goal_rtn,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    settings::RPOCHOMPSettings=RPOCHOMPSettings(n_iters=cfg.n_iters),
    optimizer::RPOTrajectoryOptimizerSettings=RPOTrajectoryOptimizerSettings(),
    max_runtime_s::Real=Inf,
    initial_path=nothing,
)
    local_cfg = rpo_pso_config(cfg; w_obs=cfg.w_obs * optimizer.w_obs_scale)
    n_internal = max(1, Int(local_cfg.n_waypoints))
    start = Vector{Float64}(start_rtn)
    goal = Vector{Float64}(goal_rtn)
    theta = rpo_trajectory_internal_points_from_seed(initial_path, start, goal, n_internal)
    lo, hi = rpo_trajectory_search_bounds(start, goal, local_cfg.search_margin_m, n_internal)
    lo .= min.(lo, theta)
    hi .= max.(hi, theta)
    _, R_inv, _ = rpo_second_difference_metric(n_internal)
    obstacle_margin = max(Float64(optimizer.obstacle_margin_m), Float64(safe_distance_m))

    objective_fn = points -> rpo_trajectory_soft_objective(
        points,
        geometry,
        local_cfg;
        safe_distance_m=safe_distance_m,
        obstacle_margin_m=obstacle_margin,
        w_smooth=settings.w_smooth,
    )

    best_points = rpo_trajectory_points_from_internal(theta, start, goal)
    best_obj = objective_fn(best_points)
    best_components = rpo_normalized_path_cost_components(best_points, geometry, local_cfg; safe_distance_m=safe_distance_m)
    best_cost = best_components.total
    history = Matrix{Float64}[]
    cost_history = Float64[]
    no_change_count = 0
    previous_obj = best_obj
    iterations_completed = 0
    start_ns = time_ns()
    max_runtime = Float64(max_runtime_s)

    for iter in 1:max(1, Int(settings.n_iters))
        if iter > 1 && isfinite(max_runtime) && (time_ns() - start_ns) / 1.0e9 >= max_runtime
            break
        end

        grad = rpo_chomp_numeric_gradient(theta, start, goal, objective_fn, settings.gradient_eps)
        cov_grad = zeros(size(grad))
        for d in 1:3
            cov_grad[d, :] .= R_inv * vec(grad[d, :])
        end

        accepted = false
        candidate_obj = best_obj
        candidate_theta = theta
        for scale in (1.0, 0.5, 0.25, 0.125, 0.0625)
            trial_theta = rpo_clamp_internal_waypoints(theta .- settings.learning_rate * scale .* cov_grad, lo, hi)
            trial_obj = objective_fn(rpo_trajectory_points_from_internal(trial_theta, start, goal))
            if isfinite(trial_obj) && trial_obj < candidate_obj
                candidate_obj = trial_obj
                candidate_theta = trial_theta
                accepted = true
                break
            end
        end

        if accepted
            theta = candidate_theta
            best_obj = candidate_obj
            best_points = rpo_trajectory_points_from_internal(theta, start, goal)
            best_components = rpo_normalized_path_cost_components(best_points, geometry, local_cfg; safe_distance_m=safe_distance_m)
            best_cost = best_components.total
        end

        iterations_completed = iter
        push!(history, copy(best_points))
        push!(cost_history, best_cost)

        unchanged = abs(previous_obj - best_obj) <= optimizer.no_change_tol * max(1.0, abs(previous_obj))
        no_change_count = unchanged ? no_change_count + 1 : 0
        previous_obj = best_obj
        if isfinite(max_runtime) && (time_ns() - start_ns) / 1.0e9 >= max_runtime
            break
        end
        if !optimizer.force_full_iters && (!accepted || no_change_count >= max(0, optimizer.no_change_iters))
            break
        end
    end

    refined, refined_cost, improved = rpo_post_refine_path(best_points, geometry, local_cfg; safe_distance_m=safe_distance_m)
    refined_components = rpo_normalized_path_cost_components(refined, geometry, local_cfg; safe_distance_m=safe_distance_m)
    return (
        path=refined,
        raw_path=best_points,
        cost=refined_cost,
        raw_cost=best_cost,
        components=refined_components,
        raw_components=best_components,
        config=local_cfg,
        adaptive=(enabled=false,),
        refinement_improved=improved,
        cost_history=cost_history,
        history=history,
        iterations=iterations_completed,
        objective=best_obj,
    )
end

"""Evaluate waypoint obstacle state cost used by STOMP rollouts."""
function rpo_stomp_waypoint_state_cost(
    points,
    internal_idx::Integer,
    geometry;
    safe_distance_m,
    w_obs,
    w_len,
    w_smooth,
    obstacle_margin_m,
)
    point_idx = internal_idx + 1
    p = points[:, point_idx]
    clearance = rpo_clearance_distance_to_station(p, geometry)
    obs = rpo_chomp_obstacle_potential(clearance, safe_distance_m, obstacle_margin_m)
    local_len = norm(p - points[:, point_idx - 1]) + norm(points[:, point_idx + 1] - p)
    d2 = points[:, point_idx - 1] .- 2.0 .* p .+ points[:, point_idx + 1]
    return w_obs * obs + w_len * local_len + w_smooth * dot(d2, d2)
end

"""Plan an RPO path using the STOMP-like comparison optimizer."""
function rpo_stomp_plan_path(
    start_rtn,
    goal_rtn,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    settings::RPOSTOMPSettings=RPOSTOMPSettings(n_iters=cfg.n_iters),
    optimizer::RPOTrajectoryOptimizerSettings=RPOTrajectoryOptimizerSettings(),
    max_runtime_s::Real=Inf,
    rng=Random.default_rng(),
    initial_path=nothing,
)
    local_cfg = rpo_pso_config(cfg; w_obs=cfg.w_obs * optimizer.w_obs_scale)
    n_internal = max(1, Int(local_cfg.n_waypoints))
    start = Vector{Float64}(start_rtn)
    goal = Vector{Float64}(goal_rtn)
    theta = rpo_trajectory_internal_points_from_seed(initial_path, start, goal, n_internal)
    lo, hi = rpo_trajectory_search_bounds(start, goal, local_cfg.search_margin_m, n_internal)
    lo .= min.(lo, theta)
    hi .= max.(hi, theta)
    _, R_inv, M = rpo_second_difference_metric(n_internal)
    cov = max(settings.noise_std, 1.0e-9)^2 .* R_inv .+ 1.0e-10 .* Matrix{Float64}(I, n_internal, n_internal)
    noise_chol = cholesky(Symmetric(cov)).L
    obstacle_margin = max(Float64(optimizer.obstacle_margin_m), Float64(safe_distance_m))

    objective_fn = points -> rpo_trajectory_soft_objective(
        points,
        geometry,
        local_cfg;
        safe_distance_m=safe_distance_m,
        obstacle_margin_m=obstacle_margin,
        w_smooth=settings.w_smooth,
    )

    best_points = rpo_trajectory_points_from_internal(theta, start, goal)
    best_obj = objective_fn(best_points)
    best_components = rpo_normalized_path_cost_components(best_points, geometry, local_cfg; safe_distance_m=safe_distance_m)
    best_cost = best_components.total
    history = Matrix{Float64}[]
    cost_history = Float64[]
    no_change_count = 0
    previous_obj = best_obj
    iterations_completed = 0
    K = max(1, Int(settings.n_rollouts))
    start_ns = time_ns()
    max_runtime = Float64(max_runtime_s)

    for iter in 1:max(1, Int(settings.n_iters))
        if iter > 1 && isfinite(max_runtime) && (time_ns() - start_ns) / 1.0e9 >= max_runtime
            break
        end

        epsilons = zeros(3, n_internal, K)
        local_costs = zeros(n_internal, K)
        rollout_thetas = Vector{Matrix{Float64}}(undef, K)
        rollout_objs = fill(Inf, K)

        for k in 1:K
            eps = zeros(3, n_internal)
            for d in 1:3
                eps[d, :] .= noise_chol * randn(rng, n_internal)
            end
            trial_theta = rpo_clamp_internal_waypoints(theta .+ eps, lo, hi)
            eps .= trial_theta .- theta
            trial_points = rpo_trajectory_points_from_internal(trial_theta, start, goal)
            rollout_thetas[k] = trial_theta
            epsilons[:, :, k] .= eps
            rollout_objs[k] = objective_fn(trial_points)
            for i in 1:n_internal
                local_costs[i, k] = rpo_stomp_waypoint_state_cost(
                    trial_points,
                    i,
                    geometry;
                    safe_distance_m=safe_distance_m,
                    w_obs=local_cfg.w_obs,
                    w_len=local_cfg.w_len,
                    w_smooth=settings.w_smooth,
                    obstacle_margin_m=obstacle_margin,
                )
            end
        end

        delta = zeros(3, n_internal)
        temp = max(settings.lambda, 1.0e-9)
        for i in 1:n_internal
            costs_i = local_costs[i, :]
            min_cost = minimum(costs_i)
            weights = exp.(-(costs_i .- min_cost) ./ temp)
            denom = sum(weights)
            if denom <= 1.0e-12 || !isfinite(denom)
                weights .= 1.0 / K
            else
                weights ./= denom
            end
            for k in 1:K
                delta[:, i] .+= weights[k] .* epsilons[:, i, k]
            end
        end

        smooth_delta = zeros(3, n_internal)
        for d in 1:3
            smooth_delta[d, :] .= M * vec(delta[d, :])
        end

        candidate_theta = theta
        candidate_obj = best_obj
        for scale in (1.0, 0.5, 0.25, 0.125)
            trial_theta = rpo_clamp_internal_waypoints(theta .+ settings.update_step * scale .* smooth_delta, lo, hi)
            trial_obj = objective_fn(rpo_trajectory_points_from_internal(trial_theta, start, goal))
            if isfinite(trial_obj) && trial_obj < candidate_obj
                candidate_theta = trial_theta
                candidate_obj = trial_obj
                break
            end
        end

        best_rollout_idx = argmin(rollout_objs)
        if rollout_objs[best_rollout_idx] < candidate_obj
            candidate_theta = rollout_thetas[best_rollout_idx]
            candidate_obj = rollout_objs[best_rollout_idx]
        end

        accepted = candidate_obj < best_obj
        if accepted
            theta = candidate_theta
            best_obj = candidate_obj
            best_points = rpo_trajectory_points_from_internal(theta, start, goal)
            best_components = rpo_normalized_path_cost_components(best_points, geometry, local_cfg; safe_distance_m=safe_distance_m)
            best_cost = best_components.total
        end

        iterations_completed = iter
        push!(history, copy(best_points))
        push!(cost_history, best_cost)

        unchanged = abs(previous_obj - best_obj) <= optimizer.no_change_tol * max(1.0, abs(previous_obj))
        no_change_count = unchanged ? no_change_count + 1 : 0
        previous_obj = best_obj
        if isfinite(max_runtime) && (time_ns() - start_ns) / 1.0e9 >= max_runtime
            break
        end
        if !optimizer.force_full_iters && no_change_count >= max(0, optimizer.no_change_iters)
            break
        end
    end

    refined, refined_cost, improved = rpo_post_refine_path(best_points, geometry, local_cfg; safe_distance_m=safe_distance_m)
    refined_components = rpo_normalized_path_cost_components(refined, geometry, local_cfg; safe_distance_m=safe_distance_m)
    return (
        path=refined,
        raw_path=best_points,
        cost=refined_cost,
        raw_cost=best_cost,
        components=refined_components,
        raw_components=best_components,
        config=local_cfg,
        adaptive=(enabled=false,),
        refinement_improved=improved,
        cost_history=cost_history,
        history=history,
        iterations=iterations_completed,
        objective=best_obj,
    )
end
