"""RRT-Connect planner settings for RPO warm starts and comparison planners."""
Base.@kwdef struct RPORRTConnectSettings
    n_iters::Int = 1000
    step_size_m::Float64 = 0.75
    goal_sample_rate::Float64 = 0.05
    collision_sample_ds_m::Float64 = 0.10
    adaptive_collision_sampling_enable::Bool = true
    collision_max_sample_ds_m::Float64 = 0.50
    collision_far_clearance_m::Float64 = 1.0
    collision_sampling_power::Float64 = 1.0
    collision_safe_distance_fraction::Float64 = 0.5
    collision_obstacle_guard_fraction::Float64 = 0.5
    connect_max_steps::Int = 10_000
    shortcut_iters::Int = 80
end

"""RRT* planner settings for RPO comparison planning."""
Base.@kwdef struct RPORRTStarSettings
    n_iters::Int = 1000
    step_size_m::Float64 = 0.75
    goal_sample_rate::Float64 = 0.05
    collision_sample_ds_m::Float64 = 0.10
    adaptive_collision_sampling_enable::Bool = true
    collision_max_sample_ds_m::Float64 = 0.50
    collision_far_clearance_m::Float64 = 1.0
    collision_sampling_power::Float64 = 1.0
    collision_safe_distance_fraction::Float64 = 0.5
    collision_obstacle_guard_fraction::Float64 = 0.5
    neighbor_radius_m::Float64 = 2.0
    shortcut_iters::Int = 80
end

"""Tree storage for RPO RRT nodes, parent links, and accumulated costs."""
struct RPORRTConnectTree
    nodes::Vector{SVector{3, Float64}}
    parents::Vector{Int}
    costs::Vector{Float64}
end

"""Tree storage for RPO RRT nodes, parent links, and accumulated costs."""
RPORRTConnectTree(root::SVector{3, Float64}) = RPORRTConnectTree([root], [0], [0.0])

"""Return the nearest node in an RPO RRT tree."""
function rpo_rrt_nearest_index(tree::RPORRTConnectTree, q::SVector{3, Float64})
    return hypr_rrt_nearest_index(tree, q)
end

"""Return RPO RRT tree nodes within a radius of the query point."""
function rpo_rrt_near_indices(tree::RPORRTConnectTree, q::SVector{3, Float64}, radius_m::Real)
    return hypr_rrt_near_indices(tree, q, max(Float64(radius_m), 0.0))
end

"""Sample a random RPO state inside axis-aligned bounds."""
function rpo_rrt_random_state(rng, lo, hi)
    return SVector{3, Float64}(
        lo[1] + rand(rng) * (hi[1] - lo[1]),
        lo[2] + rand(rng) * (hi[2] - lo[2]),
        lo[3] + rand(rng) * (hi[3] - lo[3]),
    )
end

"""Advance an RPO RRT state toward a target by one step."""
function rpo_rrt_steer(q_near::SVector{3, Float64}, q_target::SVector{3, Float64}, step_size_m::Real)
    return hypr_rrt_steer(q_near, q_target, step_size_m)
end

"""Choose collision-check spacing for RPO RRT segments."""
function rpo_rrt_collision_min_ds_m(
    sample_ds_m::Real,
    geometry,
    safe_distance_m::Real,
    safe_distance_fraction::Real,
    guard_fraction::Real,
)
    min_ds = max(Float64(sample_ds_m), 1.0e-9)
    if Float64(safe_distance_m) > 0.0
        min_ds = min(min_ds, Float64(safe_distance_fraction) * Float64(safe_distance_m))
    end
    inflated_radius = rpo_inflated_obstacle_radius_m(geometry, safe_distance_m)
    inflated_radius > 0.0 || return min_ds
    return max(min(min_ds, Float64(guard_fraction) * inflated_radius), 1.0e-9)
end

"""Check whether an RPO RRT edge is collision-free under the active safety margin."""
function rpo_rrt_segment_is_safe(
    q_from,
    q_to,
    geometry;
    safe_distance_m::Real,
    sample_ds_m::Real,
    adaptive_enable::Bool=true,
    max_sample_ds_m::Real=0.50,
    far_clearance_m::Real=1.0,
    sampling_power::Real=1.0,
    safe_distance_fraction::Real=0.5,
    obstacle_guard_fraction::Real=0.5,
)
    q0 = SVector{3, Float64}(q_from)
    q1 = SVector{3, Float64}(q_to)
    safe = Float64(safe_distance_m)
    samples = if adaptive_enable
        min_ds = rpo_rrt_collision_min_ds_m(
            sample_ds_m,
            geometry,
            safe,
            safe_distance_fraction,
            obstacle_guard_fraction,
        )
        rpo_adaptive_segment_samples(
            q0,
            q1,
            geometry;
            safe_distance_m=safe,
            min_ds_m=min_ds,
            max_ds_m=max(Float64(max_sample_ds_m), min_ds),
            far_clearance_m=far_clearance_m,
            power=sampling_power,
        )
    else
        dist = norm(q1 - q0)
        n = max(1, Int(ceil(dist / max(Float64(sample_ds_m), 1.0e-6))))
        out = zeros(3, n + 1)
        @inbounds for i in 0:n
            α = i / n
            out[:, i + 1] .= (1.0 - α) * q0 + α * q1
        end
        out
    end
    @inbounds for i in 1:size(samples, 2)
        q = SVector{3, Float64}(samples[:, i])
        rpo_clearance_distance_to_station(q, geometry) + 1.0e-9 >= safe || return false
    end
    return true
end

"""Check whether an RPO RRT edge is collision-free under the active safety margin."""
function rpo_rrt_segment_is_safe(q_from, q_to, geometry, settings; safe_distance_m::Real)
    return rpo_rrt_segment_is_safe(
        q_from,
        q_to,
        geometry;
        safe_distance_m=safe_distance_m,
        sample_ds_m=settings.collision_sample_ds_m,
        adaptive_enable=settings.adaptive_collision_sampling_enable,
        max_sample_ds_m=settings.collision_max_sample_ds_m,
        far_clearance_m=settings.collision_far_clearance_m,
        sampling_power=settings.collision_sampling_power,
        safe_distance_fraction=settings.collision_safe_distance_fraction,
        obstacle_guard_fraction=settings.collision_obstacle_guard_fraction,
    )
end

"""Extend an RPO RRT tree one step toward a target state."""
function rpo_rrt_extend!(
    tree::RPORRTConnectTree,
    q_target::SVector{3, Float64},
    geometry,
    settings::RPORRTConnectSettings;
    safe_distance_m::Real,
)
    nearest_idx = rpo_rrt_nearest_index(tree, q_target)
    q_new, status = rpo_rrt_steer(tree.nodes[nearest_idx], q_target, settings.step_size_m)
    status == :trapped && return :trapped, nearest_idx
    if !rpo_rrt_segment_is_safe(
        tree.nodes[nearest_idx],
        q_new,
        geometry,
        settings;
        safe_distance_m=safe_distance_m,
    )
        return :trapped, nearest_idx
    end
    push!(tree.nodes, q_new)
    push!(tree.parents, nearest_idx)
    push!(tree.costs, tree.costs[nearest_idx] + norm(q_new - tree.nodes[nearest_idx]))
    return status, length(tree.nodes)
end

"""Repeatedly extend an RPO RRT tree toward a target until blocked or connected."""
function rpo_rrt_connect!(
    tree::RPORRTConnectTree,
    q_target::SVector{3, Float64},
    geometry,
    settings::RPORRTConnectSettings;
    safe_distance_m::Real,
)
    status = :advanced
    idx = 1
    steps = 0
    while status == :advanced && steps < max(1, settings.connect_max_steps)
        status, idx = rpo_rrt_extend!(
            tree,
            q_target,
            geometry,
            settings;
            safe_distance_m=safe_distance_m,
        )
        steps += 1
    end
    return status, idx
end

"""Reconstruct an RPO RRT path from a tree node back to the root."""
function rpo_rrt_tree_path(tree::RPORRTConnectTree, idx::Integer)
    return hypr_rrt_tree_path(tree, idx)
end

"""Refresh RPO RRT accumulated costs below a rewired parent node."""
function rpo_rrt_refresh_subtree_costs!(tree::RPORRTConnectTree, parent_idx::Integer)
    return hypr_rrt_refresh_subtree_costs!(tree, parent_idx)
end

"""Join the two sides of an RPO bidirectional RRT connection."""
function rpo_rrt_connect_join_paths(start_tree::RPORRTConnectTree, start_idx::Integer, goal_tree::RPORRTConnectTree, goal_idx::Integer)
    return hypr_rrt_join_paths(start_tree, start_idx, goal_tree, goal_idx)
end

"""Randomly shortcut an RPO RRT path while preserving safety."""
function rpo_rrt_shortcut_path(path, geometry, settings::RPORRTConnectSettings; safe_distance_m::Real, rng=Random.default_rng())
    pts = Matrix{Float64}(path)
    size(pts, 2) <= 2 && return pts
    for _ in 1:max(0, settings.shortcut_iters)
        n_pts = size(pts, 2)
        n_pts <= 2 && break
        i = rand(rng, 1:(n_pts - 2))
        j = rand(rng, (i + 2):n_pts)
        if rpo_rrt_segment_is_safe(
            pts[:, i],
            pts[:, j],
            geometry,
            settings;
            safe_distance_m=safe_distance_m,
        )
            keep = vcat(1:i, j:n_pts)
            candidate = pts[:, keep]
            rpo_path_length(candidate) <= rpo_path_length(pts) + 1.0e-9 && (pts = candidate)
        end
    end
    return pts
end

"""Add an RRT* node and optionally rewire nearby nodes through the cheaper parent."""
function rpo_rrt_star_add_node!(
    tree::RPORRTConnectTree,
    q_new::SVector{3, Float64},
    geometry,
    settings::RPORRTStarSettings;
    safe_distance_m::Real,
)
    nearest_idx = rpo_rrt_nearest_index(tree, q_new)
    near_idxs = rpo_rrt_near_indices(tree, q_new, settings.neighbor_radius_m)
    isempty(near_idxs) && push!(near_idxs, nearest_idx)

    best_parent = nearest_idx
    best_cost = tree.costs[nearest_idx] + norm(q_new - tree.nodes[nearest_idx])
    @inbounds for idx in near_idxs
        edge_cost = norm(q_new - tree.nodes[idx])
        candidate_cost = tree.costs[idx] + edge_cost
        if candidate_cost + 1.0e-9 < best_cost &&
            rpo_rrt_segment_is_safe(
                tree.nodes[idx],
                q_new,
                geometry,
                settings;
                safe_distance_m=safe_distance_m,
            )
            best_parent = idx
            best_cost = candidate_cost
        end
    end

    if !rpo_rrt_segment_is_safe(
        tree.nodes[best_parent],
        q_new,
        geometry,
        settings;
        safe_distance_m=safe_distance_m,
    )
        return 0
    end

    push!(tree.nodes, q_new)
    push!(tree.parents, best_parent)
    push!(tree.costs, best_cost)
    new_idx = length(tree.nodes)

    @inbounds for idx in near_idxs
        idx == best_parent && continue
        idx == new_idx && continue
        candidate_cost = best_cost + norm(tree.nodes[idx] - q_new)
        if candidate_cost + 1.0e-9 < tree.costs[idx] &&
            rpo_rrt_segment_is_safe(
                q_new,
                tree.nodes[idx],
                geometry,
                settings;
                safe_distance_m=safe_distance_m,
            )
            tree.parents[idx] = new_idx
            tree.costs[idx] = candidate_cost
            rpo_rrt_refresh_subtree_costs!(tree, idx)
        end
    end
    return new_idx
end

"""Plan an RPO path with bidirectional RRT-Connect."""
function rpo_rrt_connect_plan_path(
    start_rtn,
    goal_rtn,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    settings::RPORRTConnectSettings=RPORRTConnectSettings(),
    max_runtime_s::Real=Inf,
    rng=Random.default_rng(),
)
    local_cfg = rpo_pso_config(cfg; curve_type=:polyline)
    start = SVector{3, Float64}(start_rtn)
    goal = SVector{3, Float64}(goal_rtn)
    lo, hi = rpo_pso_bounds(start, goal, local_cfg)
    direct_path = hcat(collect(start), collect(goal))
    iterations_completed = 0
    start_ns = time_ns()
    max_runtime = Float64(max_runtime_s)

    if rpo_rrt_segment_is_safe(
        start,
        goal,
        geometry,
        settings;
        safe_distance_m=safe_distance_m,
    )
        components = rpo_normalized_path_cost_components(direct_path, geometry, local_cfg; safe_distance_m=safe_distance_m)
        return (
            path=direct_path,
            raw_path=direct_path,
            cost=components.total,
            raw_cost=components.total,
            components=components,
            raw_components=components,
            config=local_cfg,
            adaptive=(enabled=false,),
            refinement_improved=false,
            cost_history=[components.total],
            history=[direct_path],
            iterations=iterations_completed,
            objective=components.total,
            path_found=true,
        )
    end

    start_tree = RPORRTConnectTree(start)
    goal_tree = RPORRTConnectTree(goal)
    found_path = nothing

    for iter in 1:max(1, settings.n_iters)
        if iter > 1 && isfinite(max_runtime) && (time_ns() - start_ns) / 1.0e9 >= max_runtime
            break
        end

        grow_from_start = isodd(iter)
        q_bias = grow_from_start ? goal : start
        q_rand = rand(rng) < clamp(settings.goal_sample_rate, 0.0, 1.0) ?
            q_bias :
            rpo_rrt_random_state(rng, lo, hi)

        tree_a = grow_from_start ? start_tree : goal_tree
        tree_b = grow_from_start ? goal_tree : start_tree
        status, idx_a = rpo_rrt_extend!(
            tree_a,
            q_rand,
            geometry,
            settings;
            safe_distance_m=safe_distance_m,
        )
        iterations_completed = iter
        status == :trapped && continue

        connect_status, idx_b = rpo_rrt_connect!(
            tree_b,
            tree_a.nodes[idx_a],
            geometry,
            settings;
            safe_distance_m=safe_distance_m,
        )
        if connect_status == :reached
            start_idx = grow_from_start ? idx_a : idx_b
            goal_idx = grow_from_start ? idx_b : idx_a
            found_path = rpo_rrt_connect_join_paths(start_tree, start_idx, goal_tree, goal_idx)
            break
        end
    end

    raw_path = found_path === nothing ? direct_path : Matrix{Float64}(found_path)
    shortcut_path = found_path === nothing ?
        raw_path :
        rpo_rrt_shortcut_path(raw_path, geometry, settings; safe_distance_m=safe_distance_m, rng=rng)
    refined, refined_cost, improved = rpo_post_refine_path(shortcut_path, geometry, local_cfg; safe_distance_m=safe_distance_m)
    raw_components = rpo_normalized_path_cost_components(raw_path, geometry, local_cfg; safe_distance_m=safe_distance_m)
    refined_components = rpo_normalized_path_cost_components(refined, geometry, local_cfg; safe_distance_m=safe_distance_m)
    return (
        path=refined,
        raw_path=raw_path,
        cost=refined_cost,
        raw_cost=raw_components.total,
        components=refined_components,
        raw_components=raw_components,
        config=local_cfg,
        adaptive=(enabled=false,),
        refinement_improved=improved,
        cost_history=[refined_components.total],
        history=[raw_path],
        iterations=iterations_completed,
        objective=refined_components.total,
        path_found=found_path !== nothing,
    )
end

"""Plan an RPO RRT-Connect path and refit it to Bezier control points."""
function rpo_rrt_connect_bezier_plan_path(
    start_rtn,
    goal_rtn,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    settings::RPORRTConnectSettings=RPORRTConnectSettings(),
    max_runtime_s::Real=Inf,
    rng=Random.default_rng(),
)
    base_plan = rpo_rrt_connect_plan_path(
        start_rtn,
        goal_rtn,
        geometry,
        cfg;
        safe_distance_m=safe_distance_m,
        settings=settings,
        max_runtime_s=max_runtime_s,
        rng=rng,
    )
    bezier_cfg = rpo_pso_config(cfg; curve_type=:bezier)
    samples = rpo_sample_path(
        base_plan.path,
        bezier_cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=bezier_cfg.sample_ds_m,
        curve_type=:polyline,
    )
    base_controls = max(2, bezier_cfg.n_waypoints + 2)
    max_controls = max(base_controls, min(size(samples, 2), max(base_controls + 6, size(base_plan.path, 2))))

    best_path = rpo_fit_bezier_fixed_endpoints(samples, base_controls, bezier_cfg)
    best_components = rpo_normalized_path_cost_components(best_path, geometry, bezier_cfg; safe_distance_m=safe_distance_m)
    for n_control in (base_controls + 1):max_controls
        candidate = rpo_fit_bezier_fixed_endpoints(samples, n_control, bezier_cfg)
        comps = rpo_normalized_path_cost_components(candidate, geometry, bezier_cfg; safe_distance_m=safe_distance_m)
        if comps.J_obs < best_components.J_obs - 1.0e-9 ||
                (comps.J_obs <= best_components.J_obs + 1.0e-9 && comps.total < best_components.total)
            best_path = candidate
            best_components = comps
        end
        best_components.violation_count == 0 && break
    end

    refined, refined_cost, improved = rpo_post_refine_path(best_path, geometry, bezier_cfg; safe_distance_m=safe_distance_m)
    refined_components = rpo_normalized_path_cost_components(refined, geometry, bezier_cfg; safe_distance_m=safe_distance_m)
    return (
        path=refined,
        raw_path=base_plan.raw_path,
        smoothed_path=best_path,
        cost=refined_cost,
        raw_cost=base_plan.raw_cost,
        components=refined_components,
        raw_components=base_plan.raw_components,
        config=bezier_cfg,
        adaptive=(enabled=false,),
        refinement_improved=improved,
        cost_history=[refined_components.total],
        history=[base_plan.raw_path],
        iterations=base_plan.iterations,
        objective=refined_components.total,
        path_found=base_plan.path_found,
    )
end

"""Plan an RPO path with RRT* and optional shortcut smoothing."""
function rpo_rrt_star_plan_path(
    start_rtn,
    goal_rtn,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    settings::RPORRTStarSettings=RPORRTStarSettings(),
    max_runtime_s::Real=Inf,
    rng=Random.default_rng(),
)
    local_cfg = rpo_pso_config(cfg; curve_type=:polyline)
    start = SVector{3, Float64}(start_rtn)
    goal = SVector{3, Float64}(goal_rtn)
    lo, hi = rpo_pso_bounds(start, goal, local_cfg)
    direct_path = hcat(collect(start), collect(goal))
    iterations_completed = 0
    start_ns = time_ns()
    max_runtime = Float64(max_runtime_s)

    if rpo_rrt_segment_is_safe(
        start,
        goal,
        geometry,
        settings;
        safe_distance_m=safe_distance_m,
    )
        components = rpo_normalized_path_cost_components(direct_path, geometry, local_cfg; safe_distance_m=safe_distance_m)
        return (
            path=direct_path,
            raw_path=direct_path,
            cost=components.total,
            raw_cost=components.total,
            components=components,
            raw_components=components,
            config=local_cfg,
            adaptive=(enabled=false,),
            refinement_improved=false,
            cost_history=[components.total],
            history=[direct_path],
            iterations=iterations_completed,
            objective=components.total,
            path_found=true,
        )
    end

    tree = RPORRTConnectTree(start)
    best_goal_idx = 0
    best_goal_cost = Inf
    cost_history = Float64[]
    history = Matrix{Float64}[]

    for iter in 1:max(1, settings.n_iters)
        if iter > 1 && isfinite(max_runtime) && (time_ns() - start_ns) / 1.0e9 >= max_runtime
            break
        end

        q_target = rand(rng) < clamp(settings.goal_sample_rate, 0.0, 1.0) ?
            goal :
            rpo_rrt_random_state(rng, lo, hi)
        nearest_idx = rpo_rrt_nearest_index(tree, q_target)
        q_new, status = rpo_rrt_steer(tree.nodes[nearest_idx], q_target, settings.step_size_m)
        iterations_completed = iter
        status == :trapped && continue

        new_idx = rpo_rrt_star_add_node!(
            tree,
            q_new,
            geometry,
            settings;
            safe_distance_m=safe_distance_m,
        )
        new_idx == 0 && continue

        goal_edge = norm(goal - tree.nodes[new_idx])
        if goal_edge <= settings.step_size_m &&
            tree.costs[new_idx] + goal_edge + 1.0e-9 < best_goal_cost &&
            rpo_rrt_segment_is_safe(
                tree.nodes[new_idx],
                goal,
                geometry,
                settings;
                safe_distance_m=safe_distance_m,
            )
            best_goal_idx = new_idx
            best_goal_cost = tree.costs[new_idx] + goal_edge
            raw_candidate = hcat(rpo_rrt_tree_path(tree, best_goal_idx), collect(goal))
            push!(history, raw_candidate)
            comps = rpo_normalized_path_cost_components(raw_candidate, geometry, local_cfg; safe_distance_m=safe_distance_m)
            push!(cost_history, comps.total)
        end
    end

    raw_path = best_goal_idx == 0 ? direct_path : hcat(rpo_rrt_tree_path(tree, best_goal_idx), collect(goal))
    shortcut_settings = RPORRTConnectSettings(
        n_iters=settings.n_iters,
        step_size_m=settings.step_size_m,
        goal_sample_rate=settings.goal_sample_rate,
        collision_sample_ds_m=settings.collision_sample_ds_m,
        adaptive_collision_sampling_enable=settings.adaptive_collision_sampling_enable,
        collision_max_sample_ds_m=settings.collision_max_sample_ds_m,
        collision_far_clearance_m=settings.collision_far_clearance_m,
        collision_sampling_power=settings.collision_sampling_power,
        collision_safe_distance_fraction=settings.collision_safe_distance_fraction,
        collision_obstacle_guard_fraction=settings.collision_obstacle_guard_fraction,
        shortcut_iters=settings.shortcut_iters,
    )
    shortcut_path = best_goal_idx == 0 ?
        raw_path :
        rpo_rrt_shortcut_path(raw_path, geometry, shortcut_settings; safe_distance_m=safe_distance_m, rng=rng)
    refined, refined_cost, improved = rpo_post_refine_path(shortcut_path, geometry, local_cfg; safe_distance_m=safe_distance_m)
    raw_components = rpo_normalized_path_cost_components(raw_path, geometry, local_cfg; safe_distance_m=safe_distance_m)
    refined_components = rpo_normalized_path_cost_components(refined, geometry, local_cfg; safe_distance_m=safe_distance_m)
    isempty(cost_history) && push!(cost_history, refined_components.total)
    isempty(history) && push!(history, raw_path)
    return (
        path=refined,
        raw_path=raw_path,
        cost=refined_cost,
        raw_cost=raw_components.total,
        components=refined_components,
        raw_components=raw_components,
        config=local_cfg,
        adaptive=(enabled=false,),
        refinement_improved=improved,
        cost_history=cost_history,
        history=history,
        iterations=iterations_completed,
        objective=refined_components.total,
        path_found=best_goal_idx != 0,
    )
end
