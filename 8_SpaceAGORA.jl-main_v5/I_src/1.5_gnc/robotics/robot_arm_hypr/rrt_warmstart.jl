"""RRT-Connect tree storage for robot-arm joint-space warm starts."""
struct RobotArmRRTConnectTree
    nodes::Vector{Vector{Float64}}
    parents::Vector{Int}
    costs::Vector{Float64}
end

"""RRT-Connect tree storage for robot-arm joint-space warm starts."""
RobotArmRRTConnectTree(root) = RobotArmRRTConnectTree([Float64.(collect(root))], [0], [0.0])

"""Return the nearest robot-arm RRT node to a joint-space query."""
function _robot_arm_rrt_nearest_index(tree::RobotArmRRTConnectTree, q::AbstractVector{<:Real})
    return hypr_rrt_nearest_index(tree, Float64.(collect(q)))
end

"""Sample a random joint-space state inside limits."""
function _robot_arm_rrt_random_state(rng, lo::AbstractVector{<:Real}, hi::AbstractVector{<:Real})
    q = zeros(length(lo))
    @inbounds for i in eachindex(q)
        q[i] = lo[i] + rand(rng) * (hi[i] - lo[i])
    end
    return q
end

"""Move one robot-arm RRT state toward another by at most one step."""
function _robot_arm_rrt_steer(q_near::AbstractVector{<:Real}, q_target::AbstractVector{<:Real}, step_size_rad::Real)
    return hypr_rrt_steer(Float64.(collect(q_near)), Float64.(collect(q_target)), step_size_rad)
end

"""Sample a robot-arm RRT edge for collision and limit checks."""
function _robot_arm_rrt_segment_samples(q_from, q_to, sample_ds_rad::Real)
    q0 = Float64.(collect(q_from))
    q1 = Float64.(collect(q_to))
    n = length(q0)
    n_steps = max(1, Int(ceil(norm(q1 - q0) / max(Float64(sample_ds_rad), 1.0e-9))))
    samples = zeros(n, n_steps + 1)
    @inbounds for i in 0:n_steps
        α = i / n_steps
        samples[:, i + 1] .= (1.0 - α) .* q0 .+ α .* q1
    end
    return samples
end

"""Check whether a robot-arm RRT edge remains collision-free."""
function _robot_arm_rrt_segment_is_safe(
    q_from,
    q_to,
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    obstacles::AbstractVector{RobotArmSphereObstacle},
    cfg::RobotArmHYPRConfig,
)
    isempty(obstacles) && return true
    samples = _robot_arm_rrt_segment_samples(q_from, q_to, cfg.rrt_warmstart_collision_sample_ds_rad)
    stats = robot_arm_clearance_stats_from_samples(model, base_pose, samples, obstacles, cfg.safe_distance_m)
    return stats.violation_count == 0 && stats.min_clearance + 1.0e-9 >= cfg.safe_distance_m
end

"""Extend a robot-arm RRT tree toward a target joint state."""
function _robot_arm_rrt_extend!(
    tree::RobotArmRRTConnectTree,
    q_target::AbstractVector{<:Real},
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    obstacles::AbstractVector{RobotArmSphereObstacle},
    cfg::RobotArmHYPRConfig,
)
    nearest_idx = _robot_arm_rrt_nearest_index(tree, q_target)
    q_new, status = _robot_arm_rrt_steer(tree.nodes[nearest_idx], q_target, cfg.rrt_warmstart_step_size_rad)
    status == :trapped && return :trapped, nearest_idx
    if !_robot_arm_rrt_segment_is_safe(tree.nodes[nearest_idx], q_new, model, base_pose, obstacles, cfg)
        return :trapped, nearest_idx
    end
    push!(tree.nodes, q_new)
    push!(tree.parents, nearest_idx)
    push!(tree.costs, tree.costs[nearest_idx] + norm(q_new - tree.nodes[nearest_idx]))
    return status, length(tree.nodes)
end

"""Connect a robot-arm RRT tree toward a target until blocked or reached."""
function _robot_arm_rrt_connect!(
    tree::RobotArmRRTConnectTree,
    q_target::AbstractVector{<:Real},
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    obstacles::AbstractVector{RobotArmSphereObstacle},
    cfg::RobotArmHYPRConfig,
)
    status = :advanced
    idx = 1
    steps = 0
    while status == :advanced && steps < cfg.rrt_warmstart_connect_max_steps
        status, idx = _robot_arm_rrt_extend!(tree, q_target, model, base_pose, obstacles, cfg)
        steps += 1
    end
    return status, idx
end

"""Reconstruct a robot-arm RRT path from a tree node."""
function _robot_arm_rrt_tree_path(tree::RobotArmRRTConnectTree, idx::Integer)
    return hypr_rrt_tree_path(tree, idx)
end

"""Join start and goal robot-arm RRT trees into one path."""
function _robot_arm_rrt_join_paths(
    start_tree::RobotArmRRTConnectTree,
    start_idx::Integer,
    goal_tree::RobotArmRRTConnectTree,
    goal_idx::Integer,
)
    return hypr_rrt_join_paths(start_tree, start_idx, goal_tree, goal_idx)
end

"""Shortcut a robot-arm RRT path while preserving clearance."""
function _robot_arm_rrt_shortcut_path(
    path,
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    obstacles::AbstractVector{RobotArmSphereObstacle},
    cfg::RobotArmHYPRConfig,
    rng,
)
    pts = Matrix{Float64}(path)
    size(pts, 2) <= 2 && return pts
    for _ in 1:cfg.rrt_warmstart_shortcut_iters
        n_pts = size(pts, 2)
        n_pts <= 2 && break
        i = rand(rng, 1:(n_pts - 2))
        j = rand(rng, (i + 2):n_pts)
        if _robot_arm_rrt_segment_is_safe(pts[:, i], pts[:, j], model, base_pose, obstacles, cfg)
            candidate = pts[:, vcat(1:i, j:n_pts)]
            _robot_arm_path_length(candidate) <= _robot_arm_path_length(pts) + 1.0e-9 && (pts = candidate)
        end
    end
    return pts
end

"""Resample a robot-arm joint polyline to a fixed number of points."""
function _robot_arm_resample_polyline_points(path, n_points::Int)
    n_points >= 2 || throw(ArgumentError("n_points must be at least 2."))
    pts = Matrix{Float64}(path)
    size(pts, 2) == n_points && return copy(pts)
    result = zeros(size(pts, 1), n_points)
    result[:, 1] .= pts[:, 1]
    result[:, end] .= pts[:, end]
    size(pts, 2) == 1 && return result

    cumulative = zeros(size(pts, 2))
    @inbounds for j in 2:size(pts, 2)
        cumulative[j] = cumulative[j - 1] + norm(pts[:, j] - pts[:, j - 1])
    end
    total = cumulative[end]
    if total <= 1.0e-12
        for j in 2:(n_points - 1)
            result[:, j] .= pts[:, 1]
        end
        return result
    end

    seg = 1
    @inbounds for j in 2:(n_points - 1)
        target_s = total * (j - 1) / (n_points - 1)
        while seg < length(cumulative) - 1 && cumulative[seg + 1] < target_s
            seg += 1
        end
        denom = max(cumulative[seg + 1] - cumulative[seg], 1.0e-12)
        α = clamp((target_s - cumulative[seg]) / denom, 0.0, 1.0)
        result[:, j] .= (1.0 - α) .* pts[:, seg] .+ α .* pts[:, seg + 1]
    end
    return result
end

"""Score a robot-arm RRT path by HYPR objective components."""
function _robot_arm_rrt_path_score(
    path,
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    obstacles::AbstractVector{RobotArmSphereObstacle},
    cfg::RobotArmHYPRConfig,
)
    samples = robot_arm_sample_hypr_path(path, cfg.n_samples; curve_type=:polyline)
    len_ref = max(norm(path[:, end] - path[:, 1]), 1.0e-6)
    stats = robot_arm_clearance_stats_from_samples(model, base_pose, samples, obstacles, cfg.safe_distance_m)
    penalty_ref = max(cfg.safe_distance_m^2, 1.0e-8)
    return (_robot_arm_path_length(samples) / len_ref)^2 +
        cfg.w_obs * (stats.violation_count + stats.clearance_penalty / penalty_ref)
end

"""Create diagnostics for a skipped or failed robot-arm RRT warm start."""
function _robot_arm_empty_rrt_warmstart_diagnostics(cfg::RobotArmHYPRConfig)
    return (
        enabled=cfg.rrt_warmstart_enable,
        attempted=false,
        path_found=false,
        iterations=0,
        cost=Inf,
        raw_cost=Inf,
        n_points=0,
    )
end

"""Normalize robot-arm RRT warm-start diagnostics into named fields."""
function _robot_arm_rrt_warmstart_fields(warmstart)
    return (
        rrt_warmstart_enabled=warmstart.enabled,
        rrt_warmstart_attempted=warmstart.attempted,
        rrt_warmstart_path_found=warmstart.path_found,
        rrt_warmstart_iterations=warmstart.iterations,
        rrt_warmstart_cost=warmstart.cost,
        rrt_warmstart_raw_cost=warmstart.raw_cost,
        rrt_warmstart_n_points=warmstart.n_points,
        rrt_warmstart=warmstart,
    )
end

"""Generate a robot-arm RRT-Connect warm-start path for HYPR."""
function _robot_arm_rrt_connect_warmstart_path(
    q_start,
    q_goal,
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    obstacles::AbstractVector{RobotArmSphereObstacle},
    cfg::RobotArmHYPRConfig,
    rng,
)
    cfg.rrt_warmstart_enable || return nothing, _robot_arm_empty_rrt_warmstart_diagnostics(cfg)
    q0 = Float64.(collect(q_start))
    qf = Float64.(collect(q_goal))
    direct_path = hcat(q0, qf)
    iterations_completed = 0
    start_ns = time_ns()
    max_runtime = Float64(cfg.rrt_warmstart_runtime_limit_s)

    if _robot_arm_rrt_segment_is_safe(q0, qf, model, base_pose, obstacles, cfg)
        cost = _robot_arm_rrt_path_score(direct_path, model, base_pose, obstacles, cfg)
        return direct_path, (
            enabled=true,
            attempted=true,
            path_found=true,
            iterations=iterations_completed,
            cost=cost,
            raw_cost=cost,
            n_points=size(direct_path, 2),
        )
    end

    lo = [joint.lower_rad for joint in model.joints]
    hi = [joint.upper_rad for joint in model.joints]
    start_tree = RobotArmRRTConnectTree(q0)
    goal_tree = RobotArmRRTConnectTree(qf)
    found_path = nothing

    for iter in 1:max(1, cfg.rrt_warmstart_iters)
        if iter > 1 && isfinite(max_runtime) && (time_ns() - start_ns) / 1.0e9 >= max_runtime
            break
        end

        grow_from_start = isodd(iter)
        q_bias = grow_from_start ? qf : q0
        q_rand = rand(rng) < cfg.rrt_warmstart_goal_sample_rate ?
            q_bias :
            _robot_arm_rrt_random_state(rng, lo, hi)

        tree_a = grow_from_start ? start_tree : goal_tree
        tree_b = grow_from_start ? goal_tree : start_tree
        status, idx_a = _robot_arm_rrt_extend!(tree_a, q_rand, model, base_pose, obstacles, cfg)
        iterations_completed = iter
        status == :trapped && continue

        connect_status, idx_b = _robot_arm_rrt_connect!(tree_b, tree_a.nodes[idx_a], model, base_pose, obstacles, cfg)
        if connect_status == :reached
            start_idx = grow_from_start ? idx_a : idx_b
            goal_idx = grow_from_start ? idx_b : idx_a
            found_path = _robot_arm_rrt_join_paths(start_tree, start_idx, goal_tree, goal_idx)
            break
        end
    end

    raw_path = found_path === nothing ? direct_path : Matrix{Float64}(found_path)
    shortcut_path = found_path === nothing ?
        raw_path :
        _robot_arm_rrt_shortcut_path(raw_path, model, base_pose, obstacles, cfg, rng)
    cost = _robot_arm_rrt_path_score(shortcut_path, model, base_pose, obstacles, cfg)
    raw_cost = _robot_arm_rrt_path_score(raw_path, model, base_pose, obstacles, cfg)
    diagnostics = (
        enabled=true,
        attempted=true,
        path_found=found_path !== nothing,
        iterations=iterations_completed,
        cost=cost,
        raw_cost=raw_cost,
        n_points=size(shortcut_path, 2),
    )
    return found_path === nothing ? nothing : shortcut_path, diagnostics
end
