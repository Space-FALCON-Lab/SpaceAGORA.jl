"""Return the Euclidean length of an RPO waypoint path."""
function rpo_path_length(points)
    return hypr_path_length(points)
end

"""Evaluate an RPO Bezier control polygon at a normalized parameter."""
function rpo_bezier_point(points, t::Real)
    return hypr_bezier_point(points, t)
end

"""Evaluate an RPO Bezier control polygon in caller-provided buffers."""
function rpo_bezier_point!(out, work, points, t::Float64)
    return hypr_bezier_point!(out, work, points, t)
end

"""Sample a Bezier RPO path using approximately uniform arc-length spacing."""
function rpo_sample_path_bezier(points, ds::Real)
    pts = Matrix{Float64}(points)
    n_pts = size(pts, 2)
    n_pts <= 1 && return pts
    dx = pts[1, end] - pts[1, 1]
    dy = pts[2, end] - pts[2, 1]
    dz = pts[3, end] - pts[3, 1]
    chord_len = sqrt(dx * dx + dy * dy + dz * dz)
    n = max(1, Int(ceil(max(rpo_path_length(pts), chord_len) / Float64(ds))))
    out = zeros(3, n + 1)
    work = similar(pts)
    point = zeros(3)
    for j in 0:n
        rpo_bezier_point!(point, work, pts, Float64(j) / n)
        out[1, j + 1] = point[1]
        out[2, j + 1] = point[2]
        out[3, j + 1] = point[3]
    end
    out[1, 1] = pts[1, 1]
    out[2, 1] = pts[2, 1]
    out[3, 1] = pts[3, 1]
    out[1, end] = pts[1, end]
    out[2, end] = pts[2, end]
    out[3, end] = pts[3, end]
    return out
end

"""Resample a polyline to a fixed number of points."""
function rpo_resample_polyline_points(points, n_samples::Int)
    pts = Matrix{Float64}(points)
    n_samples = max(n_samples, 2)
    s = rpo_arc_length_params(pts)
    out = zeros(3, n_samples)
    total = s[end]
    for j in 1:n_samples
        out[:, j] .= rpo_interpolate_along_path(pts, s, total * (j - 1) / (n_samples - 1))
    end
    return out
end

"""Sample a polyline RPO path at the requested spacing."""
function rpo_sample_path_polyline(points, ds::Real)
    pts = Matrix{Float64}(points)
    total = rpo_path_length(pts)
    n = max(2, Int(ceil(total / Float64(ds))) + 1)
    return rpo_resample_polyline_points(pts, n)
end

"""Return the station keepout radius inflated by the requested safety margin."""
function rpo_inflated_obstacle_radius_m(geometry, safe_distance_m::Real)
    return geometry.station.keepout_radius_m +
        maximum(geometry.chaser.half_extents_body) +
        max(Float64(safe_distance_m), 0.0)
end

"""Compute the minimum adaptive sampling spacing allowed near obstacles."""
function rpo_adaptive_sampling_min_ds_m(
    base_ds::Real,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
)
    min_ds = max(Float64(base_ds), 1.0e-9)
    cfg.adaptive_sampling_enable || return min_ds
    inflated_radius = rpo_inflated_obstacle_radius_m(geometry, safe_distance_m)
    if Float64(safe_distance_m) > 0.0
        min_ds = min(min_ds, cfg.adaptive_sampling_safe_distance_fraction * Float64(safe_distance_m))
    end
    if inflated_radius > 0.0
        min_ds = min(min_ds, cfg.adaptive_sampling_obstacle_guard_fraction * inflated_radius)
    end
    return max(min_ds, 1.0e-9)
end

"""Choose a local adaptive sampling step from clearance, geometry, and speed estimates."""
function rpo_adaptive_sampling_step_m(
    clearance::Real,
    min_ds::Real,
    max_ds::Real,
    far_clearance_m::Real,
    power::Real;
    safe_distance_m::Real=0.0,
)
    min_step = max(Float64(min_ds), 1.0e-9)
    max_step = max(Float64(max_ds), min_step)
    far = max(Float64(far_clearance_m), min_step)
    clearance_excess = max(Float64(clearance) - Float64(safe_distance_m), 0.0)
    u = clamp(clearance_excess / far, 0.0, 1.0)
    step = min_step + (u^Float64(power)) * (max_step - min_step)
    return min(step, max_step, clearance_excess + min_step)
end

"""Sample one segment with adaptive spacing and preserve endpoint coverage."""
function rpo_adaptive_segment_samples(
    a,
    b,
    geometry;
    safe_distance_m::Real=0.0,
    min_ds_m::Real,
    max_ds_m::Real,
    far_clearance_m::Real,
    power::Real=1.0,
)
    q0 = SVector{3, Float64}(a)
    q1 = SVector{3, Float64}(b)
    delta = q1 - q0
    dist = norm(delta)
    dist <= eps(Float64) && return reshape(collect(q0), 3, 1)

    dir = delta / dist
    samples = Vector{SVector{3, Float64}}()
    push!(samples, q0)
    s = 0.0
    max_steps = max(2, Int(ceil(dist / max(Float64(min_ds_m), 1.0e-9))) + 2)
    steps = 0
    while s < dist - 1.0e-12 && steps < max_steps
        q = q0 + s * dir
        clearance = rpo_clearance_distance_to_station(q, geometry)
        ds = rpo_adaptive_sampling_step_m(
            clearance,
            min_ds_m,
            max_ds_m,
            far_clearance_m,
            power;
            safe_distance_m=safe_distance_m,
        )
        s = min(dist, s + ds)
        push!(samples, q0 + s * dir)
        steps += 1
    end
    samples[end] = q1

    out = zeros(3, length(samples))
    @inbounds for (j, q) in enumerate(samples)
        out[:, j] .= q
    end
    return out
end

"""Sample a polyline path with clearance-adaptive spacing."""
function rpo_sample_path_polyline_adaptive(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    base_ds_m::Real=cfg.sample_ds_m,
)
    pts = Matrix{Float64}(points)
    size(pts, 2) <= 1 && return pts
    min_ds = rpo_adaptive_sampling_min_ds_m(base_ds_m, geometry, cfg; safe_distance_m=safe_distance_m)
    max_ds = max(cfg.adaptive_sampling_max_ds_m, min_ds)
    samples = Matrix{Float64}[]
    @inbounds for j in 1:(size(pts, 2) - 1)
        seg = rpo_adaptive_segment_samples(
            pts[:, j],
            pts[:, j + 1],
            geometry;
            safe_distance_m=safe_distance_m,
            min_ds_m=min_ds,
            max_ds_m=max_ds,
            far_clearance_m=cfg.adaptive_sampling_far_clearance_m,
            power=cfg.adaptive_sampling_power,
        )
        j > 1 && (seg = seg[:, 2:end])
        push!(samples, seg)
    end
    return hcat(samples...)
end

"""Estimate Bezier speed at a normalized parameter using a finite local chord."""
function rpo_bezier_speed_estimate(points, work, point, t::Float64)
    t2 = min(1.0, t + 1.0e-4)
    t2 == t && (t2 = max(0.0, t - 1.0e-4))
    t2 == t && return 0.0
    p2 = zeros(3)
    rpo_bezier_point!(p2, work, points, t2)
    return norm(p2 - point) / abs(t2 - t)
end

"""Sample a Bezier path with adaptive spacing tied to local clearance and speed."""
function rpo_sample_path_bezier_adaptive(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    base_ds_m::Real=cfg.sample_ds_m,
)
    pts = Matrix{Float64}(points)
    size(pts, 2) <= 1 && return pts
    min_ds = rpo_adaptive_sampling_min_ds_m(base_ds_m, geometry, cfg; safe_distance_m=safe_distance_m)
    max_ds = max(cfg.adaptive_sampling_max_ds_m, min_ds)
    length_ref = max(rpo_path_length(pts), norm(pts[:, end] - pts[:, 1]), min_ds)
    samples = Vector{Vector{Float64}}()
    work = similar(pts)
    point = zeros(3)
    t = 0.0
    rpo_bezier_point!(point, work, pts, t)
    push!(samples, copy(point))
    max_steps = max(2, Int(ceil(length_ref / min_ds)) + 2)
    steps = 0
    while t < 1.0 - 1.0e-12 && steps < max_steps
        clearance = rpo_clearance_distance_to_station(point, geometry)
        ds = rpo_adaptive_sampling_step_m(
            clearance,
            min_ds,
            max_ds,
            cfg.adaptive_sampling_far_clearance_m,
            cfg.adaptive_sampling_power;
            safe_distance_m=safe_distance_m,
        )
        speed = max(rpo_bezier_speed_estimate(pts, work, point, t), length_ref, 1.0e-9)
        dt = clamp(ds / speed, 1.0e-6, 1.0 - t)
        candidate = zeros(3)
        while true
            rpo_bezier_point!(candidate, work, pts, t + dt)
            if norm(candidate - point) <= 1.25 * ds || dt <= 1.0e-6 || t + dt >= 1.0
                break
            end
            dt *= 0.5
        end
        t += dt
        point .= candidate
        push!(samples, copy(point))
        steps += 1
    end
    if t < 1.0
        rpo_bezier_point!(point, work, pts, 1.0)
        push!(samples, copy(point))
    end

    out = zeros(3, length(samples))
    @inbounds for (j, q) in enumerate(samples)
        out[:, j] .= q
    end
    out[:, 1] .= pts[:, 1]
    out[:, end] .= pts[:, end]
    return out
end

"""Sample an RPO candidate path using the configured curve representation and spacing policy."""
function rpo_sample_path(points, ds::Real; curve_type::Symbol=:bezier)
    curve_type == :bezier && return rpo_sample_path_bezier(points, ds)
    curve_type == :polyline && return rpo_sample_path_polyline(points, ds)
    throw(ArgumentError("Unsupported RPO path curve_type=$(curve_type). Use :bezier or :polyline."))
end

"""Sample an RPO candidate path using the configured curve representation and spacing policy."""
function rpo_sample_path(
    points,
    cfg::RPOPSOConfig,
    geometry;
    safe_distance_m::Real=cfg.safe_distance_m,
    base_ds_m::Real=cfg.sample_ds_m,
    curve_type::Symbol=cfg.curve_type,
)
    if !cfg.adaptive_sampling_enable
        return rpo_sample_path(points, base_ds_m; curve_type=curve_type)
    end
    curve_type == :bezier && return rpo_sample_path_bezier_adaptive(
        points,
        geometry,
        cfg;
        safe_distance_m=safe_distance_m,
        base_ds_m=base_ds_m,
    )
    curve_type == :polyline && return rpo_sample_path_polyline_adaptive(
        points,
        geometry,
        cfg;
        safe_distance_m=safe_distance_m,
        base_ds_m=base_ds_m,
    )
    throw(ArgumentError("Unsupported RPO path curve_type=$(curve_type). Use :bezier or :polyline."))
end
