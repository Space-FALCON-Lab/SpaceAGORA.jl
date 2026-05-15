function rpo_refinement_config(cfg::RPOPSOConfig)
    ds = min(cfg.sample_ds_m, cfg.refinement_sample_ds_m)
    return rpo_pso_config(cfg; sample_ds_m=ds)
end

function rpo_refinement_better(candidate, current, cfg::RPOPSOConfig)
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

function rpo_most_unsafe_sample(path, geometry, cfg::RPOPSOConfig; safe_distance_m::Real=0.0)
    samples = rpo_sample_path(path, cfg.sample_ds_m; curve_type=cfg.curve_type)
    min_clearance = Inf
    min_sample = samples[:, 1]
    min_nearest = samples[:, 1]
    @inbounds for j in 1:size(samples, 2)
        stats = rpo_clearance_to_station(SVector{3, Float64}(samples[:, j]), geometry)
        if stats.clearance < min_clearance
            min_clearance = stats.clearance
            min_sample = samples[:, j]
            min_nearest = stats.nearest_point
        end
    end
    return (sample=min_sample, nearest=min_nearest, clearance=min_clearance)
end

function rpo_repair_unsafe_waypoint(path, geometry, cfg::RPOPSOConfig, current_components; safe_distance_m::Real=0.0)
    current_components.J_obs <= 0 && return nothing
    pts = Matrix{Float64}(path)
    size(pts, 2) <= 2 && return nothing
    unsafe = rpo_most_unsafe_sample(pts, geometry, cfg; safe_distance_m=safe_distance_m)
    deficit = Float64(safe_distance_m) - unsafe.clearance
    deficit <= 0.0 && return nothing
    direction = unsafe.sample .- unsafe.nearest
    nrm = norm(direction)
    nrm <= 1.0e-9 && return nothing
    direction ./= nrm
    interior = 2:(size(pts, 2) - 1)
    nearest_wp = interior[argmin([norm(pts[:, j] - unsafe.sample) for j in interior])]
    candidate = copy(pts)
    step = min(max(deficit + cfg.refinement_sample_ds_m, cfg.refinement_sample_ds_m), cfg.refinement_merge_distance_m)
    candidate[:, nearest_wp] .+= step .* direction
    return candidate
end

function rpo_try_remove_waypoints(path, geometry, cfg::RPOPSOConfig, current_components; safe_distance_m::Real=0.0)
    current = Matrix{Float64}(path)
    size(current, 2) <= 2 && return current, current_components, false
    for j in 2:(size(current, 2) - 1)
        keep = [k for k in 1:size(current, 2) if k != j]
        candidate = current[:, keep]
        comps = rpo_normalized_path_cost_components(candidate, geometry, cfg; safe_distance_m=safe_distance_m)
        if rpo_refinement_better(comps, current_components, cfg)
            return candidate, comps, true
        end
    end
    return current, current_components, false
end

function rpo_post_refine_path(path, geometry, cfg::RPOPSOConfig; safe_distance_m::Real=0.0)
    current = Matrix{Float64}(path)
    decision_cfg = rpo_refinement_config(cfg)
    current_components = rpo_normalized_path_cost_components(current, geometry, decision_cfg; safe_distance_m=safe_distance_m)
    if !cfg.refinement_enable || cfg.refinement_rounds == 0
        return current, rpo_path_cost(current, geometry, cfg; safe_distance_m=safe_distance_m), false
    end

    improved = false
    for _ in 1:cfg.refinement_rounds
        changed = false

        repair = rpo_repair_unsafe_waypoint(
            current,
            geometry,
            decision_cfg,
            current_components;
            safe_distance_m=safe_distance_m,
        )
        if repair !== nothing
            comps = rpo_normalized_path_cost_components(repair, geometry, decision_cfg; safe_distance_m=safe_distance_m)
            if rpo_refinement_better(comps, current_components, decision_cfg)
                current = repair
                current_components = comps
                improved = true
                changed = true
            end
        end

        removed_path, removed_components, removed = rpo_try_remove_waypoints(
            current,
            geometry,
            decision_cfg,
            current_components;
            safe_distance_m=safe_distance_m,
        )
        if removed
            current = removed_path
            current_components = removed_components
            improved = true
            changed = true
        end

        changed || break
    end

    final_cost = rpo_path_cost(current, geometry, cfg; safe_distance_m=safe_distance_m)
    return current, final_cost, improved
end
