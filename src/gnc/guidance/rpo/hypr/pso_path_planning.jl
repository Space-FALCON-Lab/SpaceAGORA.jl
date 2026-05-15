function rpo_pso_iteration_weights(cfg::RPOPSOConfig, iter::Int)
    if !cfg.schedule_enable || cfg.n_iters <= 1
        return (w_inertia=cfg.w_inertia, c1=cfg.c1, c2=cfg.c2)
    end
    transition_fraction = clamp(cfg.schedule_transition_fraction, 1.0e-6, 1.0)
    progress = clamp((iter - 1) / ((cfg.n_iters - 1) * transition_fraction), 0.0, 1.0)
    smooth = progress^2 * (3.0 - 2.0 * progress)
    w_end = max(cfg.schedule_w_min, cfg.schedule_w_end_fraction * cfg.w_inertia)
    c1_end = clamp(cfg.schedule_c1_end_fraction * cfg.c1, cfg.schedule_c_min, cfg.schedule_c_max)
    c2_end = clamp(cfg.schedule_c2_end_fraction * cfg.c2, cfg.schedule_c_min, cfg.schedule_c_max)
    return (
        w_inertia=cfg.w_inertia + smooth * (w_end - cfg.w_inertia),
        c1=cfg.c1 + smooth * (c1_end - cfg.c1),
        c2=cfg.c2 + smooth * (c2_end - cfg.c2),
    )
end

function rpo_pso_cull_swarm!(
    positions,
    velocities,
    pbest,
    pbest_cost,
    lo_rep,
    hi_rep,
    gbest,
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

function rpo_pso_effective_safe_distance(cfg::RPOPSOConfig, safe_distance_m)
    safe_distance_m === nothing && return cfg.safe_distance_m
    safe = Float64(safe_distance_m)
    return safe > 0.0 || cfg.safe_distance_m <= 0.0 ? safe : cfg.safe_distance_m
end

function rpo_pso_protected_particle_mask(costs, elite_fraction)
    n = length(costs)
    mask = falses(n)
    n == 0 && return mask
    elite_count = clamp(Int(ceil(clamp(elite_fraction, 0.0, 1.0) * n)), 1, n)
    ranked = sortperm(costs)
    for idx in ranked[1:elite_count]
        mask[idx] = true
    end
    return mask
end

function rpo_pso_plan_path(start_rtn, goal_rtn, geometry, base_cfg::RPOPSOConfig; safe_distance_m=nothing, rng=Random.default_rng())
    adaptive_safe_distance = safe_distance_m === nothing ||
        (Float64(safe_distance_m) == 0.0 && base_cfg.safe_distance_m > 0.0) ?
        base_cfg.safe_distance_m :
        Float64(safe_distance_m)
    cfg, adaptive = rpo_adaptive_pso_config(base_cfg, start_rtn, goal_rtn, geometry; safe_distance_m=adaptive_safe_distance)
    validate_rpo_pso_config(cfg)
    effective_safe_distance = rpo_pso_effective_safe_distance(cfg, safe_distance_m)

    current_n_waypoints = max(0, cfg.n_waypoints)
    if current_n_waypoints == 0
        path = hcat(SVector{3, Float64}(start_rtn), SVector{3, Float64}(goal_rtn))
        comps = rpo_normalized_path_cost_components(path, geometry, cfg; safe_distance_m=effective_safe_distance)
        return (path=path, cost=comps.total, components=comps, config=cfg, adaptive=adaptive, cost_history=[comps.total])
    end

    start = SVector{3, Float64}(start_rtn)
    goal = SVector{3, Float64}(goal_rtn)
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
    thread_rngs = [MersenneTwister(rand(rng, UInt)) for _ in 1:max(1, maxthreadid())]

    function cfg_for_current()
        return rpo_pso_config(cfg; n_waypoints=current_n_waypoints, search_margin_m=current_search_margin)
    end

    function seed_control_points(seed_points, n_waypoints)
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
        seeded = rpo_resample_polyline_points(seed_points, n_waypoints + 2)
        seeded[:, 1] .= start
        seeded[:, end] .= goal
        return seeded
    end

    function reset_swarm!(new_n_waypoints, new_search_margin; seed_points=nothing)
        current_n_waypoints = max(1, Int(new_n_waypoints))
        current_search_margin = max(0.0, Float64(new_search_margin))
        local_cfg = cfg_for_current()
        lo, hi = rpo_pso_bounds(start, goal, local_cfg)
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

        seeded = seed_control_points(seed_points, current_n_waypoints)
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

    function evaluate_position(pos; cost_cutoff=Inf)
        path = rpo_position_to_path(pos, start, goal, current_n_waypoints)
        return rpo_normalized_path_cost_components(
            path,
            geometry,
            cfg;
            safe_distance_m=effective_safe_distance,
            cost_cutoff=cost_cutoff,
        )
    end

    function evaluate_swarm!()
        curr_cost = fill(Inf, cfg.n_particles)
        curr_obs = fill(Inf, cfg.n_particles)
        curr_components = Vector{Any}(undef, cfg.n_particles)
        @threads for pidx in 1:cfg.n_particles
            comps = evaluate_position(view(positions, :, pidx); cost_cutoff=pbest_cost[pidx])
            curr_cost[pidx] = comps.total
            curr_obs[pidx] = comps.J_obs
            curr_components[pidx] = comps
        end
        @inbounds for pidx in 1:cfg.n_particles
            comps = curr_components[pidx]
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
        end
        return curr_cost, curr_obs
    end

    function reexplore_waypoint_count()
        scaled = Int(ceil(max(1.0, cfg.reexplore_waypoint_scale) * current_n_waypoints))
        incremented = current_n_waypoints + max(1, Int(cfg.reexplore_waypoint_increment))
        target = max(current_n_waypoints + 1, scaled, incremented)
        cfg.reexplore_max_waypoints > 0 && (target = min(target, cfg.reexplore_max_waypoints))
        return target
    end

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
                    stagnation_count[pidx] = 0
                    accepted = true
                    break
                end
                n_trials += 1
                n_trials >= max_blocks && break
            end
            accepted || (stagnation_count[pidx] = 0)
        end
        return nothing
    end

    reset_swarm!(current_n_waypoints, current_search_margin)
    curr_cost, curr_obs = evaluate_swarm!()

    @inbounds for iter in 1:cfg.n_iters
        if cfg.reexplore_enable &&
            iter >= max(1, cfg.reexplore_trigger_iter) &&
            gbest_components !== nothing &&
            gbest_components.J_obs > 0

            seed_path = rpo_position_to_path(gbest, start, goal, current_n_waypoints)
            new_n_waypoints = reexplore_waypoint_count()
            margin_scale = max(1.0, cfg.reexplore_search_margin_scale)
            new_search_margin = max(
                current_search_margin * margin_scale,
                current_search_margin + max(cfg.sample_ds_m, effective_safe_distance, 1.0e-6),
            )
            if new_n_waypoints > current_n_waypoints || new_search_margin > current_search_margin
                reset_swarm!(new_n_waypoints, new_search_margin; seed_points=seed_path)
                curr_cost, curr_obs = evaluate_swarm!()
            end
        end

        apply_stagnation_learning!(curr_cost, curr_obs)

        rpo_pso_cull_swarm!(
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

        push!(cost_history, gbest_cost)

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
        if iter < cfg.n_iters
            curr_cost, curr_obs = evaluate_swarm!()
        end
    end

    final_cfg = cfg_for_current()
    path = rpo_position_to_path(gbest, start, goal, current_n_waypoints)
    refined, refined_cost, improved = rpo_post_refine_path(path, geometry, final_cfg; safe_distance_m=effective_safe_distance)
    comps = rpo_normalized_path_cost_components(refined, geometry, final_cfg; safe_distance_m=effective_safe_distance)
    return (
        path=refined,
        cost=refined_cost,
        components=comps,
        config=final_cfg,
        adaptive=adaptive,
        refinement_improved=improved,
        cost_history=cost_history,
    )
end
