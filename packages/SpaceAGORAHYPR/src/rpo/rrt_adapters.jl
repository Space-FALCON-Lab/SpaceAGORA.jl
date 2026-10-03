# HYPR owns configuration, objective selection, refinement and Bezier fitting.
# RRT owns its search settings, trees and search loop.
"""Restore the existing configured result fields and field order."""
function _rpo_rrt_configured_result(result, cfg::RPOPSOConfig)
    return (path=result.path, raw_path=result.raw_path, cost=result.cost,
        raw_cost=result.raw_cost, components=result.components,
        raw_components=result.raw_components, config=cfg, adaptive=result.adaptive,
        refinement_improved=result.refinement_improved, cost_history=result.cost_history,
        history=result.history, iterations=result.iterations, objective=result.objective,
        path_found=result.path_found)
end

"""Plan RRT-Connect with the existing HYPR-configured objective and refinement."""
function rpo_rrt_connect_plan_path(start_rtn, goal_rtn, geometry, cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    settings::RPORRTConnectSettings=RPORRTConnectSettings(),
    max_runtime_s::Real=Inf, rng=Random.default_rng(), post_refine::Bool=true)
    local_cfg = rpo_pso_config(cfg; curve_type=:polyline)
    start = SVector{3, Float64}(start_rtn)
    goal = SVector{3, Float64}(goal_rtn)
    bounds = rpo_pso_bounds(start, goal, local_cfg)
    result = rpo_rrt_connect_plan_path(start, goal, geometry;
        bounds=bounds, settings=settings, safe_distance_m=safe_distance_m,
        max_runtime_s=max_runtime_s, rng=rng,
        evaluate_components=p -> rpo_normalized_path_cost_components(p, geometry, local_cfg; safe_distance_m=safe_distance_m),
        evaluate_cost=p -> rpo_path_cost(p, geometry, local_cfg; safe_distance_m=safe_distance_m),
        refine_path=post_refine ? p -> rpo_post_refine_path(p, geometry, local_cfg; safe_distance_m=safe_distance_m) : nothing)
    return _rpo_rrt_configured_result(result, local_cfg)
end

"""Plan RRT* with the existing HYPR-configured objective and refinement."""
function rpo_rrt_star_plan_path(start_rtn, goal_rtn, geometry, cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    settings::RPORRTStarSettings=RPORRTStarSettings(),
    max_runtime_s::Real=Inf, rng=Random.default_rng())
    local_cfg = rpo_pso_config(cfg; curve_type=:polyline)
    start = SVector{3, Float64}(start_rtn)
    goal = SVector{3, Float64}(goal_rtn)
    bounds = rpo_pso_bounds(start, goal, local_cfg)
    result = rpo_rrt_star_plan_path(start, goal, geometry;
        bounds=bounds, settings=settings, safe_distance_m=safe_distance_m,
        max_runtime_s=max_runtime_s, rng=rng,
        evaluate_components=p -> rpo_normalized_path_cost_components(p, geometry, local_cfg; safe_distance_m=safe_distance_m),
        evaluate_cost=p -> rpo_path_cost(p, geometry, local_cfg; safe_distance_m=safe_distance_m),
        refine_path=p -> rpo_post_refine_path(p, geometry, local_cfg; safe_distance_m=safe_distance_m))
    return _rpo_rrt_configured_result(result, local_cfg)
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

