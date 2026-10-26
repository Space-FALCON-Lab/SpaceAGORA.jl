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
    observation_window_s::Float64
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
function RPOReplanningConfig(; enabled::Bool=true, spheres=RPOReplanningSphere[], desired_goal_rtn=nothing, goal_change_tolerance_m::Real=1.0e-6, safe_distance_m::Real=0.0, retime_clearance_m::Real=0.25, min_replan_interval_s::Real=0.0, hysteresis_samples::Integer=1, sphere_surface_samples::Integer=96, remaining_sample_ds_m::Real=0.10, tracking_error_retime_m::Real=Inf, tracking_error_replan_m::Real=Inf, rng_seed::Integer=740, observation_window_s::Real=0.0)
    window = Float64(observation_window_s)
    isfinite(window) && window >= 0 || throw(ArgumentError("observation_window_s must be finite and nonnegative."))
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
        window,
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
            appear_time_s=tt,
            disappear_time_s=sphere.disappear_time_s + config.observation_window_s,
            velocity_rtn_mps=sphere.velocity_rtn_mps,
            label=sphere.label,
        )
        for sphere in config.spheres
        if tt >= sphere.appear_time_s && tt <= sphere.disappear_time_s + config.observation_window_s
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

"""Unexpired reference, including the interpolated position at the current time."""
function rpo_future_reference_path(plan::RPOPlan, t::Real)
    n = size(plan.r_ref_rtn, 2)
    n == 0 && return zeros(3, 0)
    elapsed = max(0.0, Float64(t) - get(plan.diagnostics, :planned_at_s, 0.0))
    times = plan.t_ref_s
    length(times) == n || throw(ArgumentError("RPO reference times and positions must match."))
    # Compare absolute deadlines to avoid subtraction roundoff at a callback tick.
    (n == 1 || Float64(t) >= get(plan.diagnostics, :planned_at_s, 0.0) + last(times)) &&
        return plan.r_ref_rtn[:, end:end]
    j = clamp(searchsortedlast(times, elapsed), 1, n - 1)
    alpha = clamp((elapsed-times[j]) / max(times[j+1]-times[j], eps(Float64)), 0.0, 1.0)
    now = (1-alpha)*plan.r_ref_rtn[:,j] + alpha*plan.r_ref_rtn[:,j+1]
    return hcat(now, plan.r_ref_rtn[:,j+1:end])
end

"""Return the unexpired path connected to the current position."""
function rpo_remaining_reference_path(plan::RPOPlan, current_rtn; sample_ds_m::Real=0.10, t::Real=get(plan.diagnostics, :planned_at_s, 0.0))
    tail = rpo_future_reference_path(plan, t)
    current = SVector{3,Float64}(current_rtn)
    size(tail,2)==0 && return reshape(collect(current),3,1)
    isfinite(sample_ds_m) && sample_ds_m > 0 || throw(ArgumentError("sample_ds_m must be positive and finite."))
    # Resampling the complete chain at uniform arc lengths can omit a corner
    # and create a shortcut through an obstacle. Subdivide each segment instead.
    chain = hcat(current, tail)
    pieces = Matrix{Float64}[]
    for j in 1:size(chain,2)-1
        segment = rpo_sample_path_polyline(chain[:,j:j+1], sample_ds_m)
        push!(pieces, j==1 ? segment : segment[:,2:end])
    end
    return hcat(pieces...)
end

"""Minimum distance to the future piecewise-linear reference (paper Algorithm 1)."""
function rpo_reference_tracking_error(plan::RPOPlan, current_rtn; t::Real=get(plan.diagnostics, :planned_at_s, 0.0))
    tail = rpo_future_reference_path(plan, t)
    n = size(tail, 2)
    n == 0 && return Inf
    p = SVector{3,Float64}(current_rtn)
    best = norm(p - tail[:,1])
    for j in 1:n-1
        a = SVector{3,Float64}(tail[:,j]); d = SVector{3,Float64}(tail[:,j+1]) - a
        alpha = clamp(dot(p-a,d) / max(dot(d,d), eps(Float64)), 0.0, 1.0)
        best = min(best, norm(p-a-alpha*d))
    end
    return best
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
    tracking_error = rpo_reference_tracking_error(plan, current_rtn; t=t)
    remaining = rpo_remaining_reference_path(plan, current_rtn; sample_ds_m=config.remaining_sample_ds_m, t=t)
    stats = rpo_path_clearance_stats(remaining, augmented; safe_distance_m=config.safe_distance_m)
    if stats.min_clearance < config.safe_distance_m || tracking_error > config.tracking_error_replan_m
        reason = stats.min_clearance < config.safe_distance_m ? :unsafe_remaining_path : :tracking_error
        return (action=:replan, reason=reason, min_clearance=stats.min_clearance, spheres=spheres, geometry=augmented)
    elseif stats.min_clearance < config.retime_clearance_m || tracking_error > config.tracking_error_retime_m
        # An expired reference has no remaining traversal to retime.
        # Safety, goal changes and excessive tracking error are checked first.
        # Reinstalling it every tick would indefinitely reset its completion clock.
        deadline = get(plan.diagnostics, :planned_at_s, 0.0) + last(plan.t_ref_s)
        if Float64(t) >= deadline &&
            tracking_error <= config.tracking_error_retime_m
            return (action=:none, reason=:reference_complete, min_clearance=stats.min_clearance, spheres=spheres, geometry=augmented)
        end
        reason = stats.min_clearance < config.retime_clearance_m ? :low_clearance_remaining_path : :tracking_error
        return (action=:retime, reason=reason, min_clearance=stats.min_clearance, spheres=spheres, geometry=augmented)
    end
    return (action=:none, reason=:clearance_ok, min_clearance=stats.min_clearance, spheres=spheres, geometry=augmented)
end

"""Retain geometric samples before time integration changes their spacing.

The planner's actual curve representation and geometric sampling policy define
this polyline. Never infer control points from a time-sampled reference. Keeping
it separate prevents repeated retiming from turning tiny time steps into new
geometric corners at the measured-state connection.
"""
function _rpo_canonical_retiming_geometry(path_rtn, geometry, cfg::RPOPSOConfig, safe_distance_m)
    size(path_rtn, 2) <= 1 && return Matrix{Float64}(path_rtn)
    # Explicit polyline vertices are geometric corners, even if separated by
    # less than one sampling step. Do not resample them out of the route.
    cfg.curve_type == :polyline && return _rpo_retime_route_vertices(path_rtn)
    samples = rpo_sample_path(path_rtn, cfg, geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=rpo_retime_sampling_ds_m(cfg, safe_distance_m),
        curve_type=cfg.curve_type)
    return _rpo_retime_route_vertices(samples)
end

"""Create a retimed RPO plan from a path after replanning."""
function rpo_plan_from_path(path_rtn, geometry, cfg::RPOPSOConfig, safe_distance_m::Real, t::Real; cost=NaN, diagnostics=NamedTuple())
    t_ref, r_ref, v_ref = rpo_reference_from_path(path_rtn, geometry, cfg; safe_distance_m=safe_distance_m)
    nominal = haskey(diagnostics, :retiming_geometry_rtn) ? diagnostics.retiming_geometry_rtn :
        _rpo_canonical_retiming_geometry(path_rtn, geometry, cfg, safe_distance_m)
    return RPOPlan(
        valid=true,
        t_ref_s=t_ref,
        r_ref_rtn=r_ref,
        v_ref_rtn=v_ref,
        path_rtn=Matrix{Float64}(path_rtn),
        cost=Float64(cost),
        diagnostics=merge((retiming_geometry_rtn=nominal,), diagnostics, (planned_at_s=Float64(t),)),
    )
end

"""Remove only monotone collinear sampling points, within the duplicate tolerance.

Time samples are not geometric waypoints. Retiming must not create a new turn
at the next time sample when that sample lies inside a straight route segment.
Every retained chord is checked against all removed points; reversals and bends
larger than the existing 1e-10 m duplicate tolerance remain explicit vertices.
"""
function _rpo_retime_route_vertices(route; tolerance_m=1.0e-10)
    n = size(route, 2)
    n <= 2 && return Matrix{Float64}(route)
    keep = falses(n)
    keep[1] = keep[end] = true
    intervals = [(1, n)]
    while !isempty(intervals)
        lo, hi = pop!(intervals)
        hi <= lo + 1 && continue
        a = SVector{3,Float64}(route[:,lo])
        d = SVector{3,Float64}(route[:,hi]) - a
        length2 = dot(d,d)
        split = 0
        if length2 <= tolerance_m^2
            # Coincident endpoints can enclose an excursion; never erase it.
            split = lo + argmax([norm(route[:,j] - a) for j in lo+1:hi-1])
            norm(route[:,split] - a) <= tolerance_m && continue
        else
            worst = tolerance_m
            last_alpha = 0.0
            for j in lo+1:hi-1
                q = SVector{3,Float64}(route[:,j]) - a
                alpha = dot(q,d) / length2
                deviation = norm(q - clamp(alpha,0.0,1.0)*d)
                if deviation > worst
                    worst = deviation
                    split = j
                end
                if alpha < last_alpha || alpha > 1.0
                    split = j
                    break
                end
                last_alpha = alpha
            end
            split == 0 && continue
        end
        keep[split] = true
        push!(intervals, (lo,split), (split,hi))
    end
    return Matrix{Float64}(route[:,keep])
end

"""Check the retained route contract before retiming or releasing a plan."""
function _rpo_retiming_geometry_valid(plan::RPOPlan)
    route = get(plan.diagnostics, :retiming_geometry_rtn, plan.r_ref_rtn)
    route isa AbstractMatrix || return false
    size(route, 1) == 3 && size(route, 2) > 0 || return false
    size(plan.r_ref_rtn, 2) > 0 || return false
    all(x -> x isa Real && isfinite(x), route) || return false
    return norm(route[:, end] - plan.r_ref_rtn[:, end]) <= 1.0e-10
end

"""Retain the current RPO route but rebuild timing from the current state."""
function rpo_retime_existing_plan(plan::RPOPlan, current_rtn, geometry, cfg::RPOPSOConfig, safe_distance_m::Real, t::Real; current_velocity_rtn=nothing)
    _rpo_retiming_geometry_valid(plan) || throw(ArgumentError(
        "Retiming geometry must be finite 3-by-N points ending at the reference goal."))
    current = SVector{3, Float64}(current_rtn)
    # Algorithm 1 retimes from geometric progress, then stitches the current
    # state. Keep this separate from its future-in-time tracking-error metric.
    route = _rpo_retime_route_vertices(get(plan.diagnostics, :retiming_geometry_rtn, plan.r_ref_rtn))
    n=size(route,2)
    n>0 || throw(ArgumentError("Cannot retime an empty reference."))
    best=Inf;segment=1;projection=copy(route[:,1])
    for j in 1:n-1
        a=SVector{3,Float64}(route[:,j]);d=SVector{3,Float64}(route[:,j+1])-a
        alpha=clamp(dot(current-a,d)/max(dot(d,d),eps(Float64)),0.,1.)
        q=a+alpha*d;distance=norm(current-q)
        if distance<best
            best=distance;segment=j;projection=collect(q)
        end
    end
    tail=n==1 ? Matrix{Float64}(route) : hcat(projection,route[:,segment+1:end])
    # The connection from the measured position consumes reference time too.
    # Overwriting only the first sample after retiming would compress an entire
    # displacement into one update. Keep the projected route separately so a
    # later retime never incorporates this temporary connector recursively.
    # The projection locates progress; it is not a waypoint to visit first.
    # Replace that synthetic first point before retiming, retaining every
    # subsequent geometric route vertex. Collinear time samples are not turns.
    # Visiting the projection explicitly can
    # create an artificial sharp corner at every update and stall progress.
    connected = norm(current - tail[:, 1]) <= 1.0e-10 ? tail :
        (size(tail, 2) == 1 ? hcat(current, tail) : hcat(current, tail[:, 2:end]))
    retimed=rpo_plan_from_path(connected,geometry,rpo_pso_config(cfg;curve_type=:polyline),safe_distance_m,t;
        cost=plan.cost,diagnostics=merge(plan.diagnostics,(replanning_action=:retime,retiming_geometry_rtn=tail)))
    if size(retimed.r_ref_rtn,2)>1
        current_velocity_rtn !== nothing && (retimed.v_ref_rtn[:,1] .= current_velocity_rtn)
    end
    return retimed
end
