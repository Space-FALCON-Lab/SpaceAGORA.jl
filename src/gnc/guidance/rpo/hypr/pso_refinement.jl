"""Return the refinement subconfiguration from an RPO PSO config."""
function rpo_refinement_config(cfg::RPOPSOConfig)
    ds = min(cfg.sample_ds_m, cfg.refinement_sample_ds_m)
    return rpo_pso_config(cfg; sample_ds_m=ds)
end

"""Compare candidate and current objective components under refinement tolerances."""
function rpo_refinement_better(candidate, current, cfg::RPOPSOConfig)
    candidate.violation_count == 0 || return false
    isfinite(candidate.total) || return false
    hasproperty(candidate, :retimed_feasible) &&
        !candidate.retimed_feasible && return false
    candidate.J_obs <= current.J_obs + 1.0e-9 || return false
    abs_improvement = current.total - candidate.total
    rel_improvement = abs_improvement / max(abs(current.total), 1.0e-12)
    return abs_improvement > cfg.refinement_min_abs_cost_improvement &&
        rel_improvement > cfg.refinement_min_rel_cost_improvement
end

"""Clamp candidate refinement waypoints into the configured search envelope."""
function rpo_refinement_clamp_path(path, cfg::RPOPSOConfig)
    pts = Matrix{Float64}(path)
    lo, hi = rpo_pso_bounds(pts[:, 1], pts[:, end], cfg)
    @inbounds for j in 2:(size(pts, 2) - 1)
        for axis in 1:3
            pts[axis, j] = clamp(pts[axis, j], lo[axis], hi[axis])
        end
    end
    return pts
end

"""Sample a straight segment for refinement collision checks."""
function rpo_refinement_segment_samples(a, b, ds::Real)
    av = SVector{3, Float64}(a)
    bv = SVector{3, Float64}(b)
    n = max(1, Int(ceil(norm(bv - av) / max(Float64(ds), 1.0e-9))))
    samples = zeros(3, n + 1)
    @inbounds for k in 0:n
        α = k / n
        samples[:, k + 1] .= (1.0 - α) * av + α * bv
    end
    return samples
end

"""Check whether a shortcut segment remains outside the keepout region."""
function rpo_refinement_segment_is_safe(a, b, geometry, cfg::RPOPSOConfig; safe_distance_m::Real=0.0)
    ds = rpo_hypr_sampling_density_m(cfg, safe_distance_m)
    samples = if cfg.adaptive_sampling_enable
        min_ds = rpo_adaptive_sampling_min_ds_m(ds, geometry, cfg; safe_distance_m=safe_distance_m)
        rpo_adaptive_segment_samples(
            a,
            b,
            geometry;
            safe_distance_m=safe_distance_m,
            min_ds_m=min_ds,
            max_ds_m=max(cfg.adaptive_sampling_max_ds_m, min_ds),
            far_clearance_m=cfg.adaptive_sampling_far_clearance_m,
            power=cfg.adaptive_sampling_power,
        )
    else
        rpo_refinement_segment_samples(a, b, ds)
    end
    required_clearance = Float64(safe_distance_m) + cfg.refinement_straight_clearance_margin_m
    stats = rpo_clearance_stats_from_samples(
        samples,
        geometry,
        required_clearance,
    )
    return stats.min_clearance + cfg.clearance_feasibility_tol_m >= required_clearance
end

"""Try shortcutting sampled path points while preserving collision safety."""
function rpo_refinement_shortcut_samples(samples, geometry, cfg::RPOPSOConfig; safe_distance_m::Real=0.0)
    current = Matrix{Float64}(samples)
    size(current, 2) <= 2 && return current
    max_passes = max(cfg.refinement_waypoint_passes, 1)
    for _ in 1:max_passes
        changed = false
        i = 1
        while i <= size(current, 2) - 2
            accepted = false
            for j in size(current, 2):-1:(i + 2)
                old_len = rpo_path_length(current[:, i:j])
                chord_len = norm(current[:, j] - current[:, i])
                if chord_len <= old_len + 1.0e-9 &&
                        rpo_refinement_segment_is_safe(current[:, i], current[:, j], geometry, cfg; safe_distance_m=safe_distance_m)
                    keep = vcat(1:i, j:size(current, 2))
                    current = current[:, keep]
                    changed = true
                    accepted = true
                    break
                end
            end
            accepted || (i += 1)
        end
        changed || break
    end
    return current
end

"""Evaluate a Bernstein basis polynomial for Bezier fitting."""
function rpo_refinement_bernstein(n::Int, j::Int, u::Float64)
    return binomial(n, j) * (1.0 - u)^(n - j) * u^j
end

"""Assign normalized fitting parameters to sampled path points."""
function rpo_refinement_sample_params(samples)
    pts = Matrix{Float64}(samples)
    n = size(pts, 2)
    n <= 1 && return zeros(n)
    s = rpo_arc_length_params(pts)
    total = s[end]
    if total <= eps(Float64)
        return collect(range(0.0, 1.0; length=n))
    end
    return s ./ total
end

"""Fit Bezier control points to samples while keeping start and goal fixed."""
function rpo_fit_bezier_fixed_endpoints(samples, n_control::Int, cfg::RPOPSOConfig)
    pts = Matrix{Float64}(samples)
    n_control = max(Int(n_control), 2)
    start = SVector{3, Float64}(pts[:, 1])
    goal = SVector{3, Float64}(pts[:, end])
    n_internal = n_control - 2
    controls = zeros(3, n_control)
    controls[:, 1] .= start
    controls[:, end] .= goal
    n_internal == 0 && return controls

    degree = n_control - 1
    u = rpo_refinement_sample_params(pts)
    A = zeros(length(u), n_internal)
    rhs = zeros(length(u), 3)
    @inbounds for i in eachindex(u)
        b0 = rpo_refinement_bernstein(degree, 0, u[i])
        bn = rpo_refinement_bernstein(degree, degree, u[i])
        rhs[i, :] .= pts[:, i] .- b0 .* start .- bn .* goal
        for j in 1:n_internal
            A[i, j] = rpo_refinement_bernstein(degree, j, u[i])
        end
    end
    normal = A' * A
    ridge = max(eps(Float64), 1.0e-10 * maximum(diag(normal); init=1.0))
    @inbounds for j in 1:n_internal
        normal[j, j] += ridge
    end
    fit = normal \ (A' * rhs)
    @inbounds for j in 1:n_internal
        controls[:, j + 1] .= fit[j, :]
    end
    return rpo_refinement_clamp_path(controls, cfg)
end

"""Fit a global Bezier seed and project its handles onto the RRT guide's safe side."""
function rpo_fit_bezier_clearance_preserving(
    samples,
    n_control::Int,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=cfg.safe_distance_m,
)
    guide_points = Matrix{Float64}(samples)
    controls = rpo_fit_bezier_fixed_endpoints(
        guide_points, n_control, cfg,
    )
    safe = max(0.0, Float64(safe_distance_m))
    safe <= 0.0 && return (
        path=controls, min_clearance=Inf, feasible=true, passes=0,
    )

    repair_ds = min(0.05, max(0.01, 0.1 * safe))
    validation_ds = min(0.02, repair_ds)
    target_clearance = safe + min(0.02, 0.1 * safe)
    n_grid = max(
        101,
        Int(ceil(rpo_path_length(guide_points) / repair_ds)) + 1,
    )
    guide = rpo_resample_polyline_points(guide_points, n_grid)
    degree = size(controls, 2) - 1
    max_passes = 32
    completed_passes = 0

    @inbounds for pass in 1:max_passes
        corrections = 0
        for sample_index in 2:(n_grid - 1)
            u = (sample_index - 1) / (n_grid - 1)
            point = rpo_bezier_point(controls, u)
            clearance = rpo_clearance_distance_to_station(point, geometry)
            clearance + 1.0e-9 >= target_clearance && continue

            guide_point = SVector{3, Float64}(view(guide, :, sample_index))
            guide_info = rpo_clearance_to_station(guide_point, geometry)
            guide_info.distance > eps(Float64) || continue
            normal = (guide_point - guide_info.nearest_point) /
                guide_info.distance
            obstacle_radius = guide_info.distance - guide_info.clearance
            deficit = obstacle_radius + target_clearance -
                dot(normal, point - guide_info.nearest_point)

            if deficit <= 0.0
                point_info = rpo_clearance_to_station(point, geometry)
                point_info.distance > eps(Float64) || continue
                normal = (point - point_info.nearest_point) /
                    point_info.distance
                obstacle_radius = point_info.distance - point_info.clearance
                deficit = obstacle_radius + target_clearance -
                    dot(normal, point - point_info.nearest_point)
            end
            deficit > 0.0 || continue

            basis_norm_sq = 0.0
            for control_index in 1:(degree - 1)
                basis = rpo_refinement_bernstein(
                    degree, control_index, u,
                )
                basis_norm_sq += basis * basis
            end
            basis_norm_sq > 1.0e-14 || continue
            for control_index in 1:(degree - 1)
                basis = rpo_refinement_bernstein(
                    degree, control_index, u,
                )
                controls[:, control_index + 1] .+=
                    (deficit * basis / basis_norm_sq) .* normal
            end
            corrections += 1
        end
        completed_passes = pass
        corrections == 0 && break
    end

    validation_samples = rpo_sample_path_bezier(
        controls, validation_ds,
    )
    stats = rpo_clearance_stats_from_samples(
        validation_samples, geometry, safe,
    )
    return (
        path=controls,
        min_clearance=stats.min_clearance,
        feasible=stats.min_clearance + cfg.clearance_feasibility_tol_m >= safe,
        passes=completed_passes,
    )
end

"""Project a point onto a segment during handle-tightening refinement."""
function rpo_refinement_project_to_segment(q, a, b)
    qv = SVector{3, Float64}(q)
    av = SVector{3, Float64}(a)
    bv = SVector{3, Float64}(b)
    ab = bv - av
    denom = dot(ab, ab)
    denom <= eps(Float64) && return av
    α = clamp(dot(qv - av, ab) / denom, 0.0, 1.0)
    return av + α * ab
end

"""Evaluate and accept a refinement candidate only when it improves the current path."""
function rpo_try_accept_refinement(
    candidate,
    geometry,
    cfg::RPOPSOConfig,
    current_components;
    safe_distance_m::Real=0.0,
    objective_evaluator=nothing,
)
    clamped = rpo_refinement_clamp_path(candidate, cfg)
    comps = rpo_path_objective_components(
        clamped,
        geometry,
        cfg;
        safe_distance_m=safe_distance_m,
        objective_evaluator=objective_evaluator,
    )
    if rpo_refinement_better(comps, current_components, cfg)
        return clamped, comps, true
    end
    return nothing, current_components, false
end

"""Shortcut a sampled path and refit it to the active Bezier degree."""
function rpo_refine_shortcut_refit(
    path,
    geometry,
    cfg::RPOPSOConfig,
    current_components;
    safe_distance_m::Real=0.0,
    objective_evaluator=nothing,
)
    samples = rpo_sample_path(
        path,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=rpo_hypr_refinement_sampling_density_m(cfg, safe_distance_m),
        curve_type=cfg.curve_type,
    )
    shortcut = rpo_refinement_shortcut_samples(samples, geometry, cfg; safe_distance_m=safe_distance_m)
    size(shortcut, 2) == size(samples, 2) && return Matrix{Float64}(path), current_components, false
    candidate = cfg.curve_type == :cubic_bezier ? shortcut :
        rpo_fit_bezier_fixed_endpoints(shortcut, size(path, 2), cfg)
    return rpo_try_accept_refinement(
        candidate,
        geometry,
        cfg,
        current_components;
        safe_distance_m=safe_distance_m,
        objective_evaluator=objective_evaluator,
    )
end

"""Pull Bezier handles toward local path chords when doing so improves the objective."""
function rpo_refine_tighten_handles(
    path,
    geometry,
    cfg::RPOPSOConfig,
    current_components;
    safe_distance_m::Real=0.0,
    objective_evaluator=nothing,
)
    current = Matrix{Float64}(path)
    size(current, 2) <= 2 && return current, current_components, false
    start = SVector{3, Float64}(current[:, 1])
    goal = SVector{3, Float64}(current[:, end])
    improved = false
    λ_values = (0.25, 0.5, 0.75)
    @inbounds for j in 2:(size(current, 2) - 1)
        α = (j - 1) / (size(current, 2) - 1)
        q = SVector{3, Float64}(current[:, j])
        accepted_handle = false
        targets = (
            rpo_refinement_project_to_segment(q, start, goal),
            rpo_refinement_project_to_segment(q, current[:, j - 1], current[:, j + 1]),
            (1.0 - α) * start + α * goal,
        )
        for target in targets
            for λ in λ_values
                candidate = copy(current)
                candidate[:, j] .= (1.0 - λ) .* q .+ λ .* target
                accepted, comps, did_accept = rpo_try_accept_refinement(
                    candidate,
                    geometry,
                    cfg,
                    current_components;
                    safe_distance_m=safe_distance_m,
                    objective_evaluator=objective_evaluator,
                )
                if did_accept
                    current = accepted
                    current_components = comps
                    improved = true
                    accepted_handle = true
                    break
                end
            end
            accepted_handle && break
        end
    end
    return current, current_components, improved
end

"""Attempt a lower-degree Bezier representation without worsening the path objective."""
function rpo_refine_lower_degree(
    path,
    geometry,
    cfg::RPOPSOConfig,
    current_components;
    safe_distance_m::Real=0.0,
    objective_evaluator=nothing,
)
    current = Matrix{Float64}(path)
    cfg.curve_type == :cubic_bezier &&
        return current, current_components, false
    size(current, 2) <= 3 && return current, current_components, false
    improved = false
    dense = rpo_sample_path(
        current,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=rpo_hypr_refinement_sampling_density_m(cfg, safe_distance_m),
        curve_type=cfg.curve_type,
    )
    min_control = 2
    for n_control in (size(current, 2) - 1):-1:min_control
        candidate = rpo_fit_bezier_fixed_endpoints(dense, n_control, cfg)
        accepted, comps, did_accept = rpo_try_accept_refinement(
            candidate,
            geometry,
            cfg,
            current_components;
            safe_distance_m=safe_distance_m,
            objective_evaluator=objective_evaluator,
        )
        if did_accept
            current = accepted
            current_components = comps
            dense = rpo_sample_path(
                current,
                cfg,
                geometry;
                safe_distance_m=safe_distance_m,
                base_ds_m=rpo_hypr_refinement_sampling_density_m(cfg, safe_distance_m),
                curve_type=cfg.curve_type,
            )
            improved = true
        end
    end
    return current, current_components, improved
end

"""Run the configured sequence of RPO post-refinement passes."""
function rpo_post_refine_path(
    path,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    objective_evaluator=nothing,
)
    current = Matrix{Float64}(path)
    decision_cfg = rpo_refinement_config(cfg)
    current_components = rpo_path_objective_components(
        current,
        geometry,
        decision_cfg;
        safe_distance_m=safe_distance_m,
        objective_evaluator=objective_evaluator,
    )
    if !cfg.refinement_enable || cfg.refinement_rounds == 0
        final_components = rpo_path_objective_components(
            current,
            geometry,
            cfg;
            safe_distance_m=safe_distance_m,
            objective_evaluator=objective_evaluator,
        )
        return current, final_components.total, false
    end

    improved = false
    for _ in 1:cfg.refinement_rounds
        changed = false

        candidate, comps, did_accept = rpo_refine_shortcut_refit(
            current,
            geometry,
            decision_cfg,
            current_components;
            safe_distance_m=safe_distance_m,
            objective_evaluator=objective_evaluator,
        )
        if did_accept
            current = candidate
            current_components = comps
            improved = true
            changed = true
        end

        candidate, comps, did_accept = rpo_refine_tighten_handles(
            current,
            geometry,
            decision_cfg,
            current_components;
            safe_distance_m=safe_distance_m,
            objective_evaluator=objective_evaluator,
        )
        if did_accept
            current = candidate
            current_components = comps
            improved = true
            changed = true
        end

        candidate, comps, did_accept = rpo_refine_lower_degree(
            current,
            geometry,
            decision_cfg,
            current_components;
            safe_distance_m=safe_distance_m,
            objective_evaluator=objective_evaluator,
        )
        if did_accept
            current = candidate
            current_components = comps
            improved = true
            changed = true
        end

        changed || break
    end

    final_components = rpo_path_objective_components(
        current,
        geometry,
        cfg;
        safe_distance_m=safe_distance_m,
        objective_evaluator=objective_evaluator,
    )
    return current, final_components.total, improved
end
