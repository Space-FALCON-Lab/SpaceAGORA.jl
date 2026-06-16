"""Dynamic spherical obstacle definition used by RPO replanning tests and demos."""
Base.@kwdef struct RPOReplanningSphere
    center_rtn::SVector{3, Float64}
    radius_m::Float64
    appear_time_s::Float64 = 0.0
    disappear_time_s::Float64 = Inf
    velocity_rtn_mps::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
    label::String = "dynamic_sphere"
end

"""Dynamic spherical obstacle definition used by RPO replanning tests and demos."""
function RPOReplanningSphere(center_rtn, radius_m::Real; appear_time_s::Real=0.0, disappear_time_s::Real=Inf, velocity_rtn_mps=zeros(3), label::AbstractString="dynamic_sphere")
    radius = Float64(radius_m)
    radius > 0.0 || throw(ArgumentError("RPOReplanningSphere radius_m must be positive."))
    appear = Float64(appear_time_s)
    disappear = Float64(disappear_time_s)
    disappear >= appear || throw(ArgumentError("RPOReplanningSphere disappear_time_s must be >= appear_time_s."))
    return RPOReplanningSphere(
        SVector{3, Float64}(center_rtn),
        radius,
        appear,
        disappear,
        SVector{3, Float64}(velocity_rtn_mps),
        String(label),
    )
end

"""Runtime replanning policy for dynamic RPO obstacles, goal changes, and tracking errors."""
struct RPOReplanningConfig
    enabled::Bool
    spheres::Vector{RPOReplanningSphere}
    desired_goal_rtn::Union{Nothing, SVector{3, Float64}}
    goal_change_tolerance_m::Float64
    safe_distance_m::Float64
    retime_clearance_m::Float64
    min_replan_interval_s::Float64
    hysteresis_samples::Int
    sphere_surface_samples::Int
    remaining_sample_ds_m::Float64
    tracking_error_retime_m::Float64
    tracking_error_replan_m::Float64
    rng_seed::Int
end

"""Read a replanning field from objects or named tuples with a default fallback."""
function _rpo_replanning_property(value, name::Symbol, default)
    return hasproperty(value, name) ? getproperty(value, name) : default
end

"""Normalize supported sphere-like inputs to an RPOReplanningSphere."""
function _rpo_replanning_sphere(value)
    value isa RPOReplanningSphere && return value
    center = _rpo_replanning_property(value, :center_rtn, _rpo_replanning_property(value, :center, nothing))
    center === nothing && throw(ArgumentError("Replanning sphere data must provide center_rtn or center."))
    radius = _rpo_replanning_property(value, :radius_m, _rpo_replanning_property(value, :radius, nothing))
    radius === nothing && throw(ArgumentError("Replanning sphere data must provide radius_m or radius."))
    return RPOReplanningSphere(
        center,
        radius;
        appear_time_s=_rpo_replanning_property(value, :appear_time_s, _rpo_replanning_property(value, :appear_time, 0.0)),
        disappear_time_s=_rpo_replanning_property(value, :disappear_time_s, _rpo_replanning_property(value, :disappear_time, Inf)),
        velocity_rtn_mps=_rpo_replanning_property(value, :velocity_rtn_mps, _rpo_replanning_property(value, :velocity, zeros(3))),
        label=String(_rpo_replanning_property(value, :label, "dynamic_sphere")),
    )
end

"""Runtime replanning policy for dynamic RPO obstacles, goal changes, and tracking errors."""
function RPOReplanningConfig(; enabled::Bool=true, spheres=RPOReplanningSphere[], desired_goal_rtn=nothing, goal_change_tolerance_m::Real=1.0e-6, safe_distance_m::Real=0.0, retime_clearance_m::Real=0.25, min_replan_interval_s::Real=0.0, hysteresis_samples::Integer=1, sphere_surface_samples::Integer=96, remaining_sample_ds_m::Real=0.10, tracking_error_retime_m::Real=Inf, tracking_error_replan_m::Real=Inf, rng_seed::Integer=740)
    sphere_vec = RPOReplanningSphere[_rpo_replanning_sphere(s) for s in spheres]
    desired_goal = desired_goal_rtn === nothing ? nothing : SVector{3, Float64}(desired_goal_rtn)
    goal_tol = Float64(goal_change_tolerance_m)
    safe = Float64(safe_distance_m)
    retime = Float64(retime_clearance_m)
    min_interval = Float64(min_replan_interval_s)
    sample_ds = Float64(remaining_sample_ds_m)
    tracking_error_retime = max(0.0, Float64(tracking_error_retime_m))
    tracking_error_replan = max(0.0, Float64(tracking_error_replan_m))
    goal_tol >= 0.0 || throw(ArgumentError("goal_change_tolerance_m must be nonnegative."))
    safe >= 0.0 || throw(ArgumentError("safe_distance_m must be nonnegative."))
    retime >= safe || throw(ArgumentError("retime_clearance_m must be >= safe_distance_m."))
    min_interval >= 0.0 || throw(ArgumentError("min_replan_interval_s must be nonnegative."))
    hysteresis_samples >= 1 || throw(ArgumentError("hysteresis_samples must be >= 1."))
    sphere_surface_samples >= 12 || throw(ArgumentError("sphere_surface_samples must be >= 12."))
    sample_ds > 0.0 || throw(ArgumentError("remaining_sample_ds_m must be positive."))
    !isnan(tracking_error_retime) || throw(ArgumentError("tracking_error_retime_m must not be NaN."))
    !isnan(tracking_error_replan) || throw(ArgumentError("tracking_error_replan_m must not be NaN."))
    return RPOReplanningConfig(
        enabled,
        sphere_vec,
        desired_goal,
        goal_tol,
        safe,
        retime,
        min_interval,
        Int(hysteresis_samples),
        Int(sphere_surface_samples),
        sample_ds,
        tracking_error_retime,
        tracking_error_replan,
        Int(rng_seed),
    )
end

"""Return a dynamic replanning sphere center at a simulation time."""
function rpo_replanning_sphere_center(sphere::RPOReplanningSphere, t::Real)
    return sphere.center_rtn + max(0.0, Float64(t) - sphere.appear_time_s) * sphere.velocity_rtn_mps
end

"""Filter replanning spheres that are active at a simulation time."""
function rpo_active_replanning_spheres(config::RPOReplanningConfig, t::Real)
    tt = Float64(t)
    return RPOReplanningSphere[
        RPOReplanningSphere(
            rpo_replanning_sphere_center(sphere, tt),
            sphere.radius_m;
            appear_time_s=sphere.appear_time_s,
            disappear_time_s=sphere.disappear_time_s,
            velocity_rtn_mps=sphere.velocity_rtn_mps,
            label=sphere.label,
        )
        for sphere in config.spheres
        if tt >= sphere.appear_time_s && tt <= sphere.disappear_time_s
    ]
end

"""Build a stable signature for the active dynamic obstacle set."""
function rpo_replanning_signature(spheres::AbstractVector{RPOReplanningSphere})
    isempty(spheres) && return UInt(0)
    acc = UInt(0x9e3779b97f4a7c15)
    for sphere in spheres
        acc ⊻= UInt(hash((sphere.label, round.(Tuple(sphere.center_rtn); digits=4), round(sphere.radius_m; digits=4))))
    end
    return acc
end

"""Sample surface points for a spherical replanning obstacle."""
function rpo_sphere_surface_points(sphere::RPOReplanningSphere; n_points::Integer=96)
    n = Int(n_points)
    pts = zeros(3, n)
    golden = pi * (3.0 - sqrt(5.0))
    for k in 0:(n - 1)
        z = 1.0 - 2.0 * (k + 0.5) / n
        r = sqrt(max(0.0, 1.0 - z * z))
        θ = golden * k
        pts[:, k + 1] .= sphere.center_rtn .+ sphere.radius_m .* SVector{3, Float64}(r * cos(θ), r * sin(θ), z)
    end
    return pts
end

"""Return RPO station geometry augmented with active dynamic obstacle samples."""
function rpo_geometry_with_replanning_spheres(geometry::RPOReferenceGeometry, spheres::AbstractVector{RPOReplanningSphere}; sphere_surface_samples::Integer=96)
    isempty(spheres) && return geometry
    n_station = size(geometry.station.points_body, 2)
    n_sphere = Int(sphere_surface_samples) * length(spheres)
    points = zeros(3, n_station + n_sphere)
    points[:, 1:n_station] .= geometry.station.points_body
    offset = n_station
    for sphere in spheres
        pts = rpo_sphere_surface_points(sphere; n_points=sphere_surface_samples)
        points[:, (offset + 1):(offset + size(pts, 2))] .= pts
        offset += size(pts, 2)
    end
    station = RPOStationGeometry(
        points;
        keepout_radius_m=geometry.station.keepout_radius_m,
        name=geometry.station.name * "+replanning",
    )
    return RPOReferenceGeometry(station; chaser=geometry.chaser)
end

"""Return the remaining reference path from the current chaser position."""
function rpo_remaining_reference_path(plan::RPOPlan, current_rtn; sample_ds_m::Real=0.10)
    current = SVector{3, Float64}(current_rtn)
    if size(plan.r_ref_rtn, 2) == 0
        return reshape(collect(current), 3, 1)
    end
    best_idx = 1
    best_dist = Inf
    @inbounds for j in 1:size(plan.r_ref_rtn, 2)
        dx = plan.r_ref_rtn[1, j] - current[1]
        dy = plan.r_ref_rtn[2, j] - current[2]
        dz = plan.r_ref_rtn[3, j] - current[3]
        d2 = dx * dx + dy * dy + dz * dz
        if d2 < best_dist
            best_dist = d2
            best_idx = j
        end
    end
    tail = plan.r_ref_rtn[:, best_idx:end]
    path = zeros(3, size(tail, 2) + 1)
    path[:, 1] .= current
    path[:, 2:end] .= tail
    return rpo_sample_path_polyline(path, sample_ds_m)
end

"""Compute the distance from the current chaser position to the active reference path."""
function rpo_reference_tracking_error(plan::RPOPlan, current_rtn)
    current = SVector{3, Float64}(current_rtn)
    size(plan.r_ref_rtn, 2) == 0 && return Inf
    best_dist = Inf
    @inbounds for j in 1:size(plan.r_ref_rtn, 2)
        dx = plan.r_ref_rtn[1, j] - current[1]
        dy = plan.r_ref_rtn[2, j] - current[2]
        dz = plan.r_ref_rtn[3, j] - current[3]
        best_dist = min(best_dist, sqrt(dx * dx + dy * dy + dz * dz))
    end
    return best_dist
end

"""Decide whether to keep, retime, or replan the active RPO path."""
function rpo_replanning_decision(plan::RPOPlan, current_rtn, geometry::RPOReferenceGeometry, config::RPOReplanningConfig, t::Real)
    if !config.enabled || !plan.valid
        return (action=:none, reason=:disabled, min_clearance=Inf, spheres=RPOReplanningSphere[], geometry=geometry)
    end
    spheres = rpo_active_replanning_spheres(config, t)
    augmented = isempty(spheres) ? geometry : rpo_geometry_with_replanning_spheres(
        geometry,
        spheres;
        sphere_surface_samples=config.sphere_surface_samples,
    )
    if config.desired_goal_rtn !== nothing
        n_ref = size(plan.r_ref_rtn, 2)
        current_goal = n_ref == 0 ? SVector{3, Float64}(current_rtn) : SVector{3, Float64}(plan.r_ref_rtn[:, n_ref])
        if norm(current_goal - config.desired_goal_rtn) > config.goal_change_tolerance_m
            return (action=:replan, reason=:goal_changed, min_clearance=Inf, spheres=spheres, geometry=augmented)
        end
    end
    if isempty(spheres)
        tracking_error = rpo_reference_tracking_error(plan, current_rtn)
        if tracking_error > config.tracking_error_replan_m
            return (action=:replan, reason=:tracking_error, min_clearance=Inf, spheres=spheres, geometry=geometry)
        end
        if tracking_error > config.tracking_error_retime_m
            return (action=:retime, reason=:tracking_error, min_clearance=Inf, spheres=spheres, geometry=geometry)
        end
        return (action=:none, reason=:no_active_spheres, min_clearance=Inf, spheres=spheres, geometry=geometry)
    end
    remaining = rpo_remaining_reference_path(plan, current_rtn; sample_ds_m=config.remaining_sample_ds_m)
    stats = rpo_path_clearance_stats(remaining, augmented; safe_distance_m=config.safe_distance_m)
    if stats.min_clearance < config.safe_distance_m
        return (action=:replan, reason=:unsafe_remaining_path, min_clearance=stats.min_clearance, spheres=spheres, geometry=augmented)
    elseif stats.min_clearance < config.retime_clearance_m
        return (action=:retime, reason=:low_clearance_remaining_path, min_clearance=stats.min_clearance, spheres=spheres, geometry=augmented)
    end
    return (action=:none, reason=:clearance_ok, min_clearance=stats.min_clearance, spheres=spheres, geometry=augmented)
end

"""Create a retimed RPO plan from a path after replanning."""
function rpo_plan_from_path(path_rtn, geometry, cfg::RPOPSOConfig, safe_distance_m::Real, t::Real; cost=NaN, diagnostics=NamedTuple())
    t_ref, r_ref, v_ref = rpo_reference_from_path(path_rtn, geometry, cfg; safe_distance_m=safe_distance_m)
    return RPOPlan(
        valid=true,
        t_ref_s=t_ref,
        r_ref_rtn=r_ref,
        v_ref_rtn=v_ref,
        path_rtn=Matrix{Float64}(path_rtn),
        cost=Float64(cost),
        diagnostics=merge((planned_at_s=Float64(t),), diagnostics),
    )
end

"""Retain the current RPO route but rebuild timing from the current state."""
function rpo_retime_existing_plan(plan::RPOPlan, current_rtn, geometry, cfg::RPOPSOConfig, safe_distance_m::Real, t::Real)
    current = SVector{3, Float64}(current_rtn)
    tail = rpo_remaining_reference_path(plan, current; sample_ds_m=cfg.sample_ds_m)
    if size(tail, 2) >= 2
        tail[:, 1] .= current
    end
    return rpo_plan_from_path(
        tail,
        geometry,
        cfg,
        safe_distance_m,
        t;
        cost=plan.cost,
        diagnostics=merge(plan.diagnostics, (replanning_action=:retime,)),
    )
end
