"""Evaluate robot-arm HYPR path cost components for a candidate joint path."""
function robot_arm_hypr_path_cost_components(
    points,
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    obstacles::AbstractVector{RobotArmSphereObstacle},
    cfg::RobotArmHYPRConfig;
    cost_cutoff::Real=Inf,
)
    samples = robot_arm_sample_hypr_path(points, cfg.n_samples; curve_type=cfg.curve_type)
    len_ref = max(norm(points[:, end] - points[:, 1]), 1.0e-6)
    J_len = _robot_arm_path_length(samples)
    J_len_norm = J_len / len_ref
    J_smooth = _robot_arm_path_smoothness(samples)
    stats = robot_arm_clearance_stats_from_samples(model, base_pose, samples, obstacles, cfg.safe_distance_m)
    penalty_ref = max(cfg.safe_distance_m^2, 1.0e-8)
    J_obs = stats.violation_count + stats.clearance_penalty / penalty_ref
    partial = cfg.w_len * J_len_norm^2 + cfg.w_smooth * J_smooth + cfg.w_obs * J_obs
    if partial > Float64(cost_cutoff)
        return (
            total=Inf,
            J_len=J_len,
            J_len_norm=J_len_norm,
            J_smooth=J_smooth,
            J_obs=J_obs,
            min_clearance=stats.min_clearance,
            violation_count=stats.violation_count,
            violation_fraction=stats.violation_fraction,
            clearance_penalty=stats.clearance_penalty,
            len_ref=len_ref,
        )
    end
    return (
        total=partial,
        J_len=J_len,
        J_len_norm=J_len_norm,
        J_smooth=J_smooth,
        J_obs=J_obs,
        min_clearance=stats.min_clearance,
        violation_count=stats.violation_count,
        violation_fraction=stats.violation_fraction,
        clearance_penalty=stats.clearance_penalty,
        len_ref=len_ref,
    )
end

"""Compare robot-arm HYPR refinement candidates under configured tolerance."""
function _robot_arm_hypr_refinement_better(candidate, current, cfg::RobotArmHYPRConfig)
    candidate.J_obs < current.J_obs && return true
    candidate.J_obs > current.J_obs + 1.0e-9 && return false
    if candidate.min_clearance > current.min_clearance + 1.0e-6 && current.J_obs > 0
        return true
    end
    abs_improvement = current.total - candidate.total
    rel_improvement = abs_improvement / max(abs(current.total), 1.0e-12)
    return abs_improvement > cfg.refinement_min_abs_cost_improvement &&
        rel_improvement > cfg.refinement_min_rel_cost_improvement
end

"""Run local shortcut-style refinement on robot-arm HYPR path points."""
function _robot_arm_hypr_post_refine_points(
    points,
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    obstacles::AbstractVector{RobotArmSphereObstacle},
    cfg::RobotArmHYPRConfig,
)
    current = Matrix{Float64}(points)
    current_components = robot_arm_hypr_path_cost_components(current, model, base_pose, obstacles, cfg)
    if !cfg.refinement_enable || cfg.refinement_rounds == 0 || size(current, 2) <= 2
        return current, current_components, false, 0
    end

    spans = max.([joint.upper_rad - joint.lower_rad for joint in model.joints], 1.0e-9)
    lo = [joint.lower_rad for joint in model.joints]
    hi = [joint.upper_rad for joint in model.joints]
    steps = cfg.refinement_step_fraction .* spans
    improved = false
    rounds_run = 0

    for round in 1:cfg.refinement_rounds
        changed = false
        @inbounds for j in 2:(size(current, 2) - 1)
            for d in axes(current, 1)
                best_value = current[d, j]
                best_components = current_components
                for direction in (-1.0, 1.0)
                    candidate = copy(current)
                    candidate[d, j] = clamp(current[d, j] + direction * steps[d], lo[d], hi[d])
                    candidate[d, j] == current[d, j] && continue
                    comps = robot_arm_hypr_path_cost_components(
                        candidate,
                        model,
                        base_pose,
                        obstacles,
                        cfg;
                        cost_cutoff=isfinite(best_components.total) ? best_components.total : Inf,
                    )
                    if _robot_arm_hypr_refinement_better(comps, best_components, cfg)
                        best_value = candidate[d, j]
                        best_components = comps
                    end
                end
                if best_value != current[d, j]
                    current[d, j] = best_value
                    current_components = best_components
                    improved = true
                    changed = true
                end
            end
        end
        rounds_run = round
        steps .*= cfg.refinement_shrink
        changed || break
    end

    return current, current_components, improved, rounds_run
end

"""Build a RobotArmPlan from a joint reference and associated timing."""
function _robot_arm_plan_from_q_reference(
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    q_start,
    q_goal,
    target,
    t_ref,
    q_ref;
    planner::Symbol=:hypr,
)
    n, nt = size(q_ref)
    dq_ref = zeros(n, nt)
    ddq_ref = zeros(n, nt)
    if nt >= 2
        dq_ref[:, 1] .= (q_ref[:, 2] - q_ref[:, 1]) ./ (t_ref[2] - t_ref[1])
        dq_ref[:, end] .= (q_ref[:, end] - q_ref[:, end - 1]) ./ (t_ref[end] - t_ref[end - 1])
    end
    if nt >= 3
        @inbounds for k in 2:(nt - 1)
            dt_prev = t_ref[k] - t_ref[k - 1]
            dt_next = t_ref[k + 1] - t_ref[k]
            dq_ref[:, k] .= (q_ref[:, k + 1] - q_ref[:, k - 1]) ./ (dt_prev + dt_next)
            ddq_ref[:, k] .= 2.0 .* (
                (q_ref[:, k + 1] - q_ref[:, k]) ./ dt_next .-
                (q_ref[:, k] - q_ref[:, k - 1]) ./ dt_prev
            ) ./ (dt_prev + dt_next)
        end
        ddq_ref[:, 1] .= ddq_ref[:, 2]
        ddq_ref[:, end] .= ddq_ref[:, end - 1]
    end
    ee_ref = zeros(3, nt)
    @inbounds for k in 1:nt
        ee_ref[:, k] .= cloth_fk(model, base_pose, q_ref[:, k]).end_effector_position
    end
    target_v = SVector{3, Float64}(target)
    final_error = norm(SVector{3, Float64}(ee_ref[:, end]) - target_v)
    return RobotArmPlan(
        model,
        base_pose,
        collect(Float64.(t_ref)),
        q_ref,
        dq_ref,
        ddq_ref,
        ee_ref,
        Float64.(collect(q_start)),
        Float64.(collect(q_goal)),
        target_v,
        final_error,
        planner,
    )
end

"""Plan a robot-arm trajectory with HYPR sampling, PSO search, optional RRT warm start, and retiming."""
function plan_robot_arm_motion_hypr(
    model::ClothArmModel,
    base_pose::ClothArmBasePose,
    q_start,
    target;
    planner_config::RobotArmPlannerConfig=RobotArmPlannerConfig(),
    hypr_config::RobotArmHYPRConfig=RobotArmHYPRConfig(),
    obstacles::AbstractVector{RobotArmSphereObstacle}=RobotArmSphereObstacle[],
    rng=Random.default_rng(),
)
    cfg = _validate_robot_arm_hypr_config(hypr_config)
    q0 = Float64.(collect(q_start))
    n = length(q0)
    length(model.joints) == n || throw(ArgumentError("q_start does not match the arm model."))
    q_goal = cloth_ik(
        model,
        base_pose,
        SVector{3, Float64}(target);
        q_seed=q0,
        max_iters=planner_config.ik_max_iters,
        position_tol_m=planner_config.ik_tol_m,
        damping=planner_config.ik_damping,
    )
    warmstart = _robot_arm_empty_rrt_warmstart_diagnostics(cfg)
    if cfg.n_waypoints == 0
        points = hcat(q0, q_goal)
        sampled = robot_arm_sample_hypr_path(points, cfg.n_samples; curve_type=:polyline)
        t_nom = _reference_times(planner_config.dt_s, planner_config.duration_s)
        q_ref = robot_arm_sample_hypr_path(points, length(t_nom); curve_type=:polyline)
        t_ref = _robot_arm_hypr_retime_reference(model, base_pose, q_ref, planner_config, cfg)
        plan = _robot_arm_plan_from_q_reference(model, base_pose, q0, q_goal, target, t_ref, q_ref; planner=:hypr)
        wrench_ratios = _robot_arm_hypr_base_wrench_ratios(model, base_pose, q_ref, t_ref, cfg)
        comps_raw = robot_arm_hypr_path_cost_components(points, model, base_pose, obstacles, cfg)
        comps = merge(comps_raw, (
            refinement_enabled=cfg.refinement_enable,
            refinement_improved=false,
            refinement_rounds=0,
            retime_enabled=cfg.retime_enable,
            retime_duration_s=t_ref[end],
            retime_nominal_duration_s=planner_config.duration_s,
            retime_scale_max=t_ref[end] / max(planner_config.duration_s, 1.0e-12),
            retime_base_force_ratio=wrench_ratios.force_ratio,
            retime_base_torque_ratio=wrench_ratios.torque_ratio,
        ), _robot_arm_rrt_warmstart_fields(warmstart))
        return RobotArmHYPRResult(plan, points, sampled, comps.total, comps, [comps.total], cfg, false, 0, 0)
    end

    lo_base = [joint.lower_rad for joint in model.joints]
    hi_base = [joint.upper_rad for joint in model.joints]
    lo_rep = repeat(lo_base, cfg.n_waypoints)
    hi_rep = repeat(hi_base, cfg.n_waypoints)
    span_rep = max.(hi_rep .- lo_rep, 1.0e-9)
    dim = n * cfg.n_waypoints

    positions = zeros(dim, cfg.n_particles)
    velocities = zeros(dim, cfg.n_particles)
    pbest = zeros(dim, cfg.n_particles)
    pbest_cost = fill(Inf, cfg.n_particles)
    pbest_components = Vector{Any}(undef, cfg.n_particles)
    gbest = zeros(dim)
    gbest_cost = Inf
    gbest_components = nothing
    cost_history = Float64[]

    warmstart_path, warmstart =
        _robot_arm_rrt_connect_warmstart_path(q0, q_goal, model, base_pose, obstacles, cfg, rng)
    seeded = warmstart_path === nothing ?
        _robot_arm_seed_control_points(q0, q_goal, cfg.n_waypoints) :
        _robot_arm_resample_polyline_points(warmstart_path, cfg.n_waypoints + 2)
    seeded[:, 1] .= q0
    seeded[:, end] .= q_goal
    base_flat = _robot_arm_flatten_internal_points(seeded)
    @inbounds for pidx in 1:cfg.n_particles
        for d in 1:dim
            if pidx == 1
                positions[d, pidx] = clamp(base_flat[d], lo_rep[d], hi_rep[d])
            else
                positions[d, pidx] = clamp(
                    base_flat[d] + cfg.spread_scale * span_rep[d] * randn(rng),
                    lo_rep[d],
                    hi_rep[d],
                )
            end
            velocities[d, pidx] = 0.1 * cfg.max_velocity_fraction * span_rep[d] * randn(rng)
            pbest[d, pidx] = positions[d, pidx]
        end
    end

    early_stop_best_cost = Inf
    early_stop_stale_iters = 0
    early_stop_iter = 0
    cull_replacements = 0

    # Evaluate one flattened robot-arm particle as HYPR path-cost components.
    function evaluate_position(pos; cost_cutoff=Inf)
        points = _robot_arm_control_points(q0, q_goal, pos, cfg.n_waypoints)
        return robot_arm_hypr_path_cost_components(points, model, base_pose, obstacles, cfg; cost_cutoff=cost_cutoff)
    end

    for iter in 1:cfg.n_iters
        @inbounds for pidx in 1:cfg.n_particles
            comps = evaluate_position(@view positions[:, pidx]; cost_cutoff=pbest_cost[pidx])
            curr_cost = comps.total
            if curr_cost < pbest_cost[pidx]
                pbest[:, pidx] .= positions[:, pidx]
                pbest_cost[pidx] = curr_cost
                pbest_components[pidx] = comps
            end
            if curr_cost < gbest_cost
                gbest .= positions[:, pidx]
                gbest_cost = curr_cost
                gbest_components = comps
            end
        end

        push!(cost_history, gbest_cost)
        if cfg.early_stopping_enable && iter >= cfg.early_stopping_min_iters &&
                _robot_arm_hypr_early_stopping_feasible(gbest_components, cfg)
            if _robot_arm_hypr_material_improvement(gbest_cost, early_stop_best_cost, cfg)
                early_stop_best_cost = gbest_cost
                early_stop_stale_iters = 0
            else
                early_stop_stale_iters += 1
            end
            if early_stop_stale_iters >= cfg.early_stopping_patience
                early_stop_iter = iter
                break
            end
        elseif iter == 1
            early_stop_best_cost = gbest_cost
        end

        cull_replacements += _robot_arm_hypr_cull_swarm!(
            positions,
            velocities,
            pbest,
            pbest_cost,
            lo_rep,
            hi_rep,
            gbest,
            cfg,
            iter,
            rng,
        )

        weights = _robot_arm_hypr_iteration_weights(cfg, iter)
        @inbounds for pidx in 1:cfg.n_particles
            for d in 1:dim
                v = weights.w_inertia * velocities[d, pidx] +
                    weights.c1 * rand(rng) * (pbest[d, pidx] - positions[d, pidx]) +
                    weights.c2 * rand(rng) * (gbest[d] - positions[d, pidx])
                vmax = cfg.max_velocity_fraction * span_rep[d]
                velocities[d, pidx] = clamp(v, -vmax, vmax)
                positions[d, pidx] = clamp(positions[d, pidx] + velocities[d, pidx], lo_rep[d], hi_rep[d])
            end
        end
    end

    best_points_raw = _robot_arm_control_points(q0, q_goal, gbest, cfg.n_waypoints)
    best_points, refined_components_raw, refinement_improved, refinement_rounds =
        _robot_arm_hypr_post_refine_points(best_points_raw, model, base_pose, obstacles, cfg)
    sampled = robot_arm_sample_hypr_path(best_points, cfg.n_samples; curve_type=cfg.curve_type)
    t_nom = _reference_times(planner_config.dt_s, planner_config.duration_s)
    q_ref = robot_arm_sample_hypr_path(best_points, length(t_nom); curve_type=cfg.curve_type)
    t_ref = _robot_arm_hypr_retime_reference(model, base_pose, q_ref, planner_config, cfg)
    plan = _robot_arm_plan_from_q_reference(model, base_pose, q0, q_goal, target, t_ref, q_ref; planner=:hypr)
    wrench_ratios = _robot_arm_hypr_base_wrench_ratios(model, base_pose, q_ref, t_ref, cfg)
    final_components = merge(refined_components_raw, (
        refinement_enabled=cfg.refinement_enable,
        refinement_improved=refinement_improved,
        refinement_rounds=refinement_rounds,
        retime_enabled=cfg.retime_enable,
        retime_duration_s=t_ref[end],
        retime_nominal_duration_s=planner_config.duration_s,
        retime_scale_max=t_ref[end] / max(planner_config.duration_s, 1.0e-12),
        retime_base_force_ratio=wrench_ratios.force_ratio,
        retime_base_torque_ratio=wrench_ratios.torque_ratio,
    ), _robot_arm_rrt_warmstart_fields(warmstart))
    return RobotArmHYPRResult(
        plan,
        best_points,
        sampled,
        final_components.total,
        final_components,
        cost_history,
        cfg,
        early_stop_iter > 0,
        early_stop_iter,
        cull_replacements,
    )
end
