"""Return the PSO inertia and acceleration weights for the current iteration."""
function rpo_pso_iteration_weights(cfg::RPOPSOConfig, iter::Int)
    return hypr_iteration_weights(
        cfg.schedule_enable,
        cfg.n_iters,
        iter,
        cfg.w_inertia,
        cfg.c1,
        cfg.c2,
        cfg.schedule_transition_fraction,
        cfg.schedule_w_min,
        cfg.schedule_w_end_fraction,
        cfg.schedule_c1_end_fraction,
        cfg.schedule_c2_end_fraction,
        cfg.schedule_c_min,
        cfg.schedule_c_max,
    )
end

"""Project a point onto a line segment for path-shaping heuristics."""
function rpo_pso_project_to_segment(q, a, b)
    av = SVector{3, Float64}(a)
    bv = SVector{3, Float64}(b)
    qv = SVector{3, Float64}(q)
    ab = bv - av
    denom = dot(ab, ab)
    denom <= eps(Float64) && return av
    α = clamp(dot(qv - av, ab) / denom, 0.0, 1.0)
    return av + α * ab
end

"""Scale waypoint reseeding noise so endpoints move less than middle waypoints."""
function rpo_pso_tapered_noise_scale(j::Int, n_waypoints::Int)
    n_waypoints <= 1 && return 1.0
    edge_distance = min(j, n_waypoints + 1 - j)
    return clamp(edge_distance / max((n_waypoints + 1) / 2, 1.0), 0.25, 1.0)
end

"""Replace weak PSO particles with corridor-biased random samples while preserving elites."""
function rpo_pso_cull_swarm!(
    positions,
    velocities,
    pbest,
    pbest_cost,
    lo_rep,
    hi_rep,
    gbest,
    start_rtn,
    goal_rtn,
    n_waypoints::Int,
    cfg::RPOPSOConfig,
    iter::Int,
    rng,
)
    if !cfg.cull_enable || iter < cfg.cull_start_iter || cfg.cull_fraction_max <= 0.0 || !isfinite(minimum(pbest_cost))
        return 0
    end
    n_replace = min(cfg.n_particles - 1, floor(Int, cfg.cull_fraction_max * cfg.n_particles))
    n_replace <= 0 && return 0
    worst = sortperm(pbest_cost; rev=true)[1:n_replace]
    if cfg.curve_type != :bezier || n_waypoints <= 0
        @inbounds for pidx in worst
            for d in axes(positions, 1)
                span = hi_rep[d] - lo_rep[d]
                jitter = cfg.cull_noise_scale * span * (rand(rng) - 0.5)
                positions[d, pidx] = clamp(gbest[d] + jitter, lo_rep[d], hi_rep[d])
                velocities[d, pidx] = 0.0
                pbest[d, pidx] = positions[d, pidx]
            end
            pbest_cost[pidx] = Inf
        end
        return n_replace
    end

    start = SVector{3, Float64}(start_rtn)
    goal = SVector{3, Float64}(goal_rtn)
    gbest_path = rpo_position_to_path(gbest, start, goal, n_waypoints)
    n_control = n_waypoints + 2
    @inbounds for pidx in worst
        for j in 1:n_waypoints
            q_star = SVector{3, Float64}(gbest_path[:, j + 1])
            α = j / (n_control - 1)
            q_line = (1.0 - α) * start + α * goal
            q_loc = rpo_pso_project_to_segment(q_star, gbest_path[:, j], gbest_path[:, j + 2])
            d_bez = 0.6 * (q_line - q_star) + 0.4 * (q_loc - q_star)
            target = q_star + d_bez
            taper = rpo_pso_tapered_noise_scale(j, n_waypoints)
            for axis in 1:3
                d = 3 * (j - 1) + axis
                span = max(hi_rep[d] - lo_rep[d], 1.0e-12)
                noise = cfg.cull_noise_scale * taper * span * randn(rng)
                old_position = positions[d, pidx]
                positions[d, pidx] = clamp(target[axis] + noise, lo_rep[d], hi_rep[d])
                velocities[d, pidx] = cfg.cull_arc_velocity_scale * (positions[d, pidx] - old_position)
                pbest[d, pidx] = positions[d, pidx]
            end
        end
        pbest_cost[pidx] = Inf
    end
    return n_replace
end

"""Resolve the effective safety distance from planner input and config defaults."""
function rpo_pso_effective_safe_distance(cfg::RPOPSOConfig, safe_distance_m)
    safe_distance_m === nothing && return cfg.safe_distance_m
    safe = Float64(safe_distance_m)
    return safe > 0.0 || cfg.safe_distance_m <= 0.0 ? safe : cfg.safe_distance_m
end

"""Create default diagnostics for a skipped or unavailable RRT warm start."""
function rpo_pso_empty_warmstart_diagnostics(cfg::RPOPSOConfig)
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

"""Generate an optional RRT-Connect seed path for the RPO PSO swarm."""
function rpo_pso_rrt_warmstart_path(start, goal, geometry, cfg::RPOPSOConfig, safe_distance_m::Real, rng)
    cfg.rrt_warmstart_enable || return nothing, rpo_pso_empty_warmstart_diagnostics(cfg)
    settings = RPORRTConnectSettings(
        n_iters=cfg.rrt_warmstart_iters,
        step_size_m=cfg.rrt_warmstart_step_size_m,
        goal_sample_rate=cfg.rrt_warmstart_goal_sample_rate,
        collision_sample_ds_m=cfg.rrt_warmstart_collision_sample_ds_m,
        connect_max_steps=cfg.rrt_warmstart_connect_max_steps,
        shortcut_iters=cfg.rrt_warmstart_shortcut_iters,
    )
    plan = rpo_rrt_connect_plan_path(
        start,
        goal,
        geometry,
        cfg;
        safe_distance_m=safe_distance_m,
        settings=settings,
        max_runtime_s=cfg.rrt_warmstart_runtime_limit_s,
        rng=rng,
    )
    diagnostics = (
        enabled=true,
        attempted=true,
        path_found=plan.path_found,
        iterations=plan.iterations,
        cost=plan.cost,
        raw_cost=plan.raw_cost,
        n_points=size(plan.path, 2),
    )
    return plan.path_found ? plan.path : nothing, diagnostics
end

"""Mark elite RPO PSO particles that culling should leave untouched."""
function rpo_pso_protected_particle_mask(costs, elite_fraction)
    return hypr_protected_particle_mask(costs, elite_fraction)
end

"""Check whether a new RPO PSO best cost materially improves over a reference cost."""
function rpo_pso_material_improvement(new_cost::Real, reference_cost::Real, cfg::RPOPSOConfig)::Bool
    return hypr_material_improvement(
        new_cost,
        reference_cost,
        cfg.early_stopping_min_abs_improvement,
        cfg.early_stopping_min_rel_improvement,
    )
end

"""Return whether path-cost components are feasible enough for early-stopping checks."""
function rpo_pso_early_stopping_feasible(components, cfg::RPOPSOConfig)::Bool
    !cfg.early_stopping_require_feasible && return true
    components === nothing && return false
    return getproperty(components, :violation_count) == 0
end

"""Return the post-attempt stagnation count; only an accepted learning attempt clears it."""
@inline function rpo_pso_stagnation_count_after_learning(count::Integer, accepted::Bool)
    return accepted ? 0 : Int(count)
end

"""Run the RPO HYPR PSO planner and return the best path, cost, and diagnostics."""
function rpo_pso_plan_path(
    start_rtn,
    goal_rtn,
    geometry,
    base_cfg::RPOPSOConfig;
    safe_distance_m=nothing,
    rng=Random.default_rng(),
    iteration_callback=nothing,
    objective_evaluator=nothing,
)
    adaptive_safe_distance = safe_distance_m === nothing ||
        (Float64(safe_distance_m) == 0.0 && base_cfg.safe_distance_m > 0.0) ?
        base_cfg.safe_distance_m :
        Float64(safe_distance_m)
    base_cfg = rpo_pso_config(base_cfg; safe_distance_m=adaptive_safe_distance)
    cfg, adaptive = rpo_adaptive_pso_config(base_cfg, start_rtn, goal_rtn, geometry; safe_distance_m=adaptive_safe_distance)
    effective_safe_distance = rpo_pso_effective_safe_distance(cfg, safe_distance_m)
    cfg = rpo_pso_config(cfg; safe_distance_m=effective_safe_distance)

    start = SVector{3, Float64}(start_rtn)
    goal = SVector{3, Float64}(goal_rtn)
    current_n_waypoints = max(0, cfg.n_waypoints)
    warmstart_path, warmstart = current_n_waypoints > 0 ?
        rpo_pso_rrt_warmstart_path(start, goal, geometry, cfg, effective_safe_distance, rng) :
        (nothing, rpo_pso_empty_warmstart_diagnostics(cfg))
    if warmstart_path !== nothing && cfg.curve_type == :bezier
        warmstart_waypoints = max(0, size(warmstart_path, 2) - 2)
        warmstart_cap = cfg.reexplore_max_waypoints > 0 ?
            max(current_n_waypoints, cfg.reexplore_max_waypoints) :
            max(current_n_waypoints, warmstart_waypoints)
        current_n_waypoints = min(max(current_n_waypoints, warmstart_waypoints), warmstart_cap)
        cfg = rpo_pso_config(cfg; n_waypoints=current_n_waypoints)
    end
    if current_n_waypoints == 0
        path = hcat(start, goal)
        comps = rpo_path_objective_components(
            path,
            geometry,
            cfg;
            safe_distance_m=effective_safe_distance,
            objective_evaluator=objective_evaluator,
        )
        iteration_callback !== nothing && iteration_callback(0, comps.total, comps)
        return (
            path=path,
            cost=comps.total,
            components=comps,
            config=cfg,
            adaptive=adaptive,
            refinement_improved=false,
            early_stopped=false,
            early_stop_iter=0,
            iteration_timed_out=false,
            iteration_timeout_iter=0,
            iteration_timeout_phase=:none,
            iteration_timeout_events=NamedTuple[],
            warmstart=warmstart,
            cost_history=[comps.total],
        )
    end

    current_search_margin = cfg.search_margin_m
    dim = 3 * current_n_waypoints
    lo_rep = zeros(dim)
    hi_rep = zeros(dim)
    positions = zeros(dim, cfg.n_particles)
    velocities = zeros(dim, cfg.n_particles)
    pbest = zeros(dim, cfg.n_particles)
    pbest_cost = fill(Inf, cfg.n_particles)
    pbest_obs = fill(Inf, cfg.n_particles)
    stagnation_count = zeros(Int, cfg.n_particles)
    gbest = zeros(dim)
    gbest_cost = Inf
    gbest_components = nothing
    cost_history = Float64[]
    early_stop_best_cost = Inf
    early_stop_stale_iters = 0
    early_stop_iter = 0
    optimization_timed_out = false
    iteration_timeout_iter = 0
    iteration_timeout_phase = :none
    iteration_timeout_events = NamedTuple[]
    thread_rngs = [MersenneTwister(rand(rng, UInt)) for _ in 1:max(1, maxthreadid())]

    # Return elapsed wall-clock seconds for the current PSO iteration.
    function iteration_elapsed_s(iter_start_ns::UInt64)
        return (time_ns() - iter_start_ns) / 1.0e9
    end

    # Check whether the per-iteration runtime budget has been exceeded.
    function iteration_timed_out(iter_start_ns::UInt64)
        limit = cfg.iteration_runtime_limit_s
        return isfinite(limit) && limit > 0.0 && iteration_elapsed_s(iter_start_ns) >= limit
    end

    # Record the global-best cost and invoke the optional iteration callback.
    function record_iteration!(iter)
        push!(cost_history, gbest_cost)
        iteration_callback !== nothing && iteration_callback(iter, gbest_cost, gbest_components)
        return nothing
    end

    # Record a timeout event and mark the current optimization as timed out.
    function record_iteration_timeout!(iter, phase::Symbol, iter_start_ns::UInt64)
        optimization_timed_out = true
        iteration_timeout_iter = iter
        iteration_timeout_phase = phase
        push!(iteration_timeout_events, (
            iter=iter,
            phase=phase,
            elapsed_s=iteration_elapsed_s(iter_start_ns),
        ))
        return nothing
    end

    # Rebuild the local PSO config for the currently active waypoint count and search margin.
    function cfg_for_current()
        return rpo_pso_config(cfg; n_waypoints=current_n_waypoints, search_margin_m=current_search_margin)
    end

    # Build seed control points from a supplied warm start or the straight start-goal chord.
    function seed_control_points(seed_points, n_waypoints, local_cfg::RPOPSOConfig; seed_curve_type::Symbol=local_cfg.curve_type)
        if seed_points === nothing
            points = zeros(3, n_waypoints + 2)
            points[:, 1] .= start
            points[:, end] .= goal
            for j in 1:n_waypoints
                alpha = j / (n_waypoints + 1)
                points[:, j + 1] .= (1.0 - alpha) .* start .+ alpha .* goal
            end
            return points
        end
        if local_cfg.curve_type == :bezier
            samples = seed_curve_type == :bezier ?
                rpo_sample_path(seed_points, local_cfg.sample_ds_m; curve_type=:bezier) :
                rpo_sample_path_polyline(seed_points, local_cfg.sample_ds_m)
            seeded = rpo_fit_bezier_fixed_endpoints(samples, n_waypoints + 2, local_cfg)
        else
            seeded = rpo_resample_polyline_points(seed_points, n_waypoints + 2)
        end
        seeded[:, 1] .= start
        seeded[:, end] .= goal
        return seeded
    end

    # Reinitialize PSO particle state after adaptive sizing or reexploration changes.
    function reset_swarm!(
        new_n_waypoints,
        new_search_margin;
        seed_points=nothing,
        seed_curve_type::Symbol=cfg.curve_type,
        use_warmstart_bounds::Bool=false,
    )
        current_n_waypoints = max(1, Int(new_n_waypoints))
        current_search_margin = max(0.0, Float64(new_search_margin))
        local_cfg = cfg_for_current()
        lo, hi = use_warmstart_bounds && seed_points !== nothing ?
            rpo_pso_warmstart_bounds(seed_points, local_cfg) :
            rpo_pso_bounds(start, goal, local_cfg)
        dim = 3 * current_n_waypoints
        lo_rep = repeat(collect(lo), current_n_waypoints)
        hi_rep = repeat(collect(hi), current_n_waypoints)
        span_rep = hi_rep .- lo_rep
        positions = zeros(dim, cfg.n_particles)
        velocities = zeros(dim, cfg.n_particles)
        pbest = zeros(dim, cfg.n_particles)
        pbest_cost = fill(Inf, cfg.n_particles)
        pbest_obs = fill(Inf, cfg.n_particles)
        stagnation_count = zeros(Int, cfg.n_particles)
        gbest = zeros(dim)
        gbest_cost = Inf
        gbest_components = nothing

        seeded = seed_control_points(seed_points, current_n_waypoints, local_cfg; seed_curve_type=seed_curve_type)
        base = zeros(dim)
        for j in 1:current_n_waypoints
            offset = 3 * (j - 1)
            base[offset + 1] = seeded[1, j + 1]
            base[offset + 2] = seeded[2, j + 1]
            base[offset + 3] = seeded[3, j + 1]
        end
        @inbounds for pidx in 1:cfg.n_particles
            for d in 1:dim
                if pidx == 1
                    positions[d, pidx] = clamp(base[d], lo_rep[d], hi_rep[d])
                else
                    positions[d, pidx] = clamp(
                        base[d] + cfg.spread_scale * span_rep[d] * randn(rng),
                        lo_rep[d],
                        hi_rep[d],
                    )
                end
                velocities[d, pidx] = 0.1 * (rand(rng) - 0.5) * span_rep[d]
            end
        end
        return nothing
    end

    # Evaluate one flattened particle position as an RPO path-cost component bundle.
    function evaluate_position(pos; cost_cutoff=Inf)
        path = rpo_position_to_path(pos, start, goal, current_n_waypoints)
        return rpo_path_objective_components(
            path,
            geometry,
            cfg;
            safe_distance_m=effective_safe_distance,
            cost_cutoff=cost_cutoff,
            objective_evaluator=objective_evaluator,
        )
    end

    # Update particle-best and global-best state from one candidate evaluation.
    function update_swarm_best_from_components!(pidx, comps)
        if comps.total < pbest_cost[pidx]
            for d in 1:dim
                pbest[d, pidx] = positions[d, pidx]
            end
            pbest_cost[pidx] = comps.total
            pbest_obs[pidx] = comps.J_obs
            stagnation_count[pidx] = 0
        else
            stagnation_count[pidx] += 1
        end
        if comps.total < gbest_cost
            for d in 1:dim
                gbest[d] = positions[d, pidx]
            end
            gbest_cost = comps.total
            gbest_components = comps
        end
        return nothing
    end

    # Evaluate all particles and update PSO best-state bookkeeping.
    function evaluate_swarm!(; iter_start_ns=nothing)
        curr_cost = fill(Inf, cfg.n_particles)
        curr_obs = fill(Inf, cfg.n_particles)
        curr_components = Vector{Any}(undef, cfg.n_particles)
        if iter_start_ns !== nothing && isfinite(cfg.iteration_runtime_limit_s) && cfg.iteration_runtime_limit_s > 0.0
            @inbounds for pidx in 1:cfg.n_particles
                iteration_timed_out(iter_start_ns) && return curr_cost, curr_obs, true
                comps = evaluate_position(view(positions, :, pidx); cost_cutoff=pbest_cost[pidx])
                curr_cost[pidx] = comps.total
                curr_obs[pidx] = comps.J_obs
                curr_components[pidx] = comps
                update_swarm_best_from_components!(pidx, comps)
            end
            return curr_cost, curr_obs, iteration_timed_out(iter_start_ns)
        end
        @threads for pidx in 1:cfg.n_particles
            comps = evaluate_position(view(positions, :, pidx); cost_cutoff=pbest_cost[pidx])
            curr_cost[pidx] = comps.total
            curr_obs[pidx] = comps.J_obs
            curr_components[pidx] = comps
        end
        @inbounds for pidx in 1:cfg.n_particles
            comps = curr_components[pidx]
            update_swarm_best_from_components!(pidx, comps)
        end
        return curr_cost, curr_obs, false
    end

    # Pick the larger waypoint count used when obstacle stagnation triggers reexploration.
    function reexplore_waypoint_count()
        scaled = Int(ceil(max(1.0, cfg.reexplore_waypoint_scale) * current_n_waypoints))
        incremented = current_n_waypoints + max(1, Int(cfg.reexplore_waypoint_increment))
        target = max(current_n_waypoints + 1, scaled, incremented)
        cfg.reexplore_max_waypoints > 0 && (target = min(target, cfg.reexplore_max_waypoints))
        return target
    end

    # Locally perturb stagnant particles around promising waypoint blocks.
    function apply_stagnation_learning!(curr_cost, curr_obs)
        if !cfg.stagnation_learning_enable || cfg.stagnation_learning_threshold <= 0
            return nothing
        end
        protected = rpo_pso_protected_particle_mask(curr_cost, cfg.stagnation_learning_elite_fraction)
        max_blocks = clamp(cfg.stagnation_learning_max_blocks, 1, max(current_n_waypoints, 1))
        candidate = zeros(dim)
        @inbounds for pidx in 1:cfg.n_particles
            if protected[pidx] || stagnation_count[pidx] < cfg.stagnation_learning_threshold
                continue
            end
            accepted = false
            n_trials = 0
            for wp_idx in randperm(rng, current_n_waypoints)
                for d in 1:dim
                    candidate[d] = positions[d, pidx]
                end
                block = (3 * (wp_idx - 1) + 1):(3 * wp_idx)
                for d in block
                    candidate[d] = gbest[d]
                end
                comps = evaluate_position(candidate; cost_cutoff=curr_cost[pidx])
                improves_cost = comps.total + 1.0e-9 < curr_cost[pidx]
                obstacle_nonworse = comps.J_obs <= curr_obs[pidx] + 1.0e-9
                if improves_cost && obstacle_nonworse
                    for d in 1:dim
                        positions[d, pidx] = candidate[d]
                    end
                    for d in block
                        velocities[d, pidx] = 0.0
                    end
                    curr_cost[pidx] = comps.total
                    curr_obs[pidx] = comps.J_obs
                    if comps.total < pbest_cost[pidx]
                        for d in 1:dim
                            pbest[d, pidx] = candidate[d]
                        end
                        pbest_cost[pidx] = comps.total
                        pbest_obs[pidx] = comps.J_obs
                    end
                    if comps.total < gbest_cost
                        for d in 1:dim
                            gbest[d] = candidate[d]
                        end
                        gbest_cost = comps.total
                        gbest_components = comps
                    end
                    accepted = true
                    break
                end
                n_trials += 1
                n_trials >= max_blocks && break
            end
            stagnation_count[pidx] = rpo_pso_stagnation_count_after_learning(stagnation_count[pidx], accepted)
        end
        return nothing
    end

    reset_swarm!(
        current_n_waypoints,
        current_search_margin;
        seed_points=warmstart_path,
        seed_curve_type=:polyline,
        use_warmstart_bounds=warmstart_path !== nothing,
    )
    curr_cost, curr_obs, _ = evaluate_swarm!()

    @inbounds for iter in 1:cfg.n_iters
        iter_start_ns = time_ns()
        if cfg.reexplore_enable &&
            iter >= max(1, cfg.reexplore_trigger_iter) &&
            gbest_components !== nothing &&
            gbest_components.violation_count > 0

            seed_path = rpo_position_to_path(gbest, start, goal, current_n_waypoints)
            new_n_waypoints = reexplore_waypoint_count()
            margin_scale = max(1.0, cfg.reexplore_search_margin_scale)
            new_search_margin = max(
                current_search_margin * margin_scale,
                current_search_margin + max(cfg.sample_ds_m, effective_safe_distance, 1.0e-6),
            )
            if new_n_waypoints > current_n_waypoints || new_search_margin > current_search_margin
                reset_swarm!(new_n_waypoints, new_search_margin; seed_points=seed_path, seed_curve_type=cfg.curve_type)
                curr_cost, curr_obs, eval_timed_out = evaluate_swarm!(iter_start_ns=iter_start_ns)
                early_stop_best_cost = Inf
                early_stop_stale_iters = 0
                if eval_timed_out
                    record_iteration!(iter)
                    record_iteration_timeout!(iter, :reexplore_evaluate, iter_start_ns)
                    break
                end
            end
        end
        if iteration_timed_out(iter_start_ns)
            record_iteration!(iter)
            record_iteration_timeout!(iter, cfg.reexplore_enable ? :reexplore : :iteration_start, iter_start_ns)
            break
        end

        apply_stagnation_learning!(curr_cost, curr_obs)
        if iteration_timed_out(iter_start_ns)
            record_iteration!(iter)
            record_iteration_timeout!(iter, :stagnation_learning, iter_start_ns)
            break
        end

        rpo_pso_cull_swarm!(
            positions,
            velocities,
            pbest,
            pbest_cost,
            lo_rep,
            hi_rep,
            gbest,
            start,
            goal,
            current_n_waypoints,
            cfg,
            iter,
            rng,
        )
        if iteration_timed_out(iter_start_ns)
            record_iteration!(iter)
            record_iteration_timeout!(iter, :cull, iter_start_ns)
            break
        end

        record_iteration!(iter)

        if rpo_pso_material_improvement(gbest_cost, early_stop_best_cost, cfg)
            early_stop_best_cost = gbest_cost
            early_stop_stale_iters = 0
        elseif cfg.early_stopping_enable &&
            cfg.early_stopping_patience > 0 &&
            iter >= cfg.early_stopping_min_iters &&
            rpo_pso_early_stopping_feasible(gbest_components, cfg)

            early_stop_stale_iters += 1
            if early_stop_stale_iters >= cfg.early_stopping_patience
                early_stop_iter = iter
                break
            end
        end
        if iteration_timed_out(iter_start_ns)
            record_iteration_timeout!(iter, :early_stopping, iter_start_ns)
            break
        end

        weights = rpo_pso_iteration_weights(cfg, iter)
        @threads for pidx in 1:cfg.n_particles
            local_rng = thread_rngs[threadid()]
            for d in 1:dim
                velocities[d, pidx] = weights.w_inertia * velocities[d, pidx] +
                    weights.c1 * rand(local_rng) * (pbest[d, pidx] - positions[d, pidx]) +
                    weights.c2 * rand(local_rng) * (gbest[d] - positions[d, pidx])
                positions[d, pidx] = clamp(positions[d, pidx] + velocities[d, pidx], lo_rep[d], hi_rep[d])
            end
        end
        if iteration_timed_out(iter_start_ns)
            record_iteration_timeout!(iter, :velocity_update, iter_start_ns)
            break
        end
        if iter < cfg.n_iters
            curr_cost, curr_obs, eval_timed_out = evaluate_swarm!(iter_start_ns=iter_start_ns)
            if eval_timed_out
                record_iteration_timeout!(iter, :swarm_evaluate, iter_start_ns)
                break
            end
        end
    end

    final_cfg = cfg_for_current()
    path = rpo_position_to_path(gbest, start, goal, current_n_waypoints)
    refined, refined_cost, improved = rpo_post_refine_path(
        path,
        geometry,
        final_cfg;
        safe_distance_m=effective_safe_distance,
        objective_evaluator=objective_evaluator,
    )
    comps = rpo_path_objective_components(
        refined,
        geometry,
        final_cfg;
        safe_distance_m=effective_safe_distance,
        objective_evaluator=objective_evaluator,
    )
    return (
        path=refined,
        cost=refined_cost,
        components=comps,
        config=final_cfg,
        adaptive=adaptive,
        refinement_improved=improved,
        early_stopped=early_stop_iter > 0,
        early_stop_iter=early_stop_iter,
        iteration_timed_out=optimization_timed_out,
        iteration_timeout_iter=iteration_timeout_iter,
        iteration_timeout_phase=iteration_timeout_phase,
        iteration_timeout_events=iteration_timeout_events,
        warmstart=warmstart,
        cost_history=cost_history,
    )
end
