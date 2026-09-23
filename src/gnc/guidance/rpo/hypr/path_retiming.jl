"""Compute cumulative arc-length parameters for sampled path points."""
function rpo_arc_length_params(points)
    pts = Matrix{Float64}(points)
    s = zeros(size(pts, 2))

    @inbounds for j in 2:size(pts, 2)
        ds = norm(pts[:, j] - pts[:, j - 1])
        s[j] = s[j - 1] + ds
    end

    return s
end


"""Drop consecutive path samples that are closer than the tolerance."""
function rpo_remove_near_duplicate_samples(points; tol::Real=1.0e-10)
    pts = Matrix{Float64}(points)
    n = size(pts, 2)

    n <= 1 && return pts

    keep = Int[1]

    @inbounds for j in 2:n
        if norm(pts[:, j] - pts[:, keep[end]]) > Float64(tol)
            push!(keep, j)
        end
    end

    if length(keep) < n
        @warn "RPO retiming removed near-duplicate path samples." removed=(n - length(keep)) kept=length(keep)
    end

    return pts[:, keep]
end


"""Interpolate a point at a requested arc-length coordinate."""
function rpo_interpolate_along_path(points, s_vals, s_query::Real)
    pts = Matrix{Float64}(points)
    s = Vector{Float64}(s_vals)
    sq = Float64(s_query)

    n = length(s)

    n == 1 && return copy(pts[:, 1])
    sq <= s[1] && return copy(pts[:, 1])
    sq >= s[end] && return copy(pts[:, end])

    idx = clamp(searchsortedlast(s, sq), 1, n - 1)

    # Move forward if this segment has zero length.
    while idx < n && s[idx + 1] - s[idx] <= eps(Float64)
        idx += 1
    end

    idx >= n && return copy(pts[:, end])

    denom = s[idx + 1] - s[idx]

    if denom <= eps(Float64)
        return copy(pts[:, idx])
    end

    α = clamp((sq - s[idx]) / denom, 0.0, 1.0)

    return (1.0 - α) .* pts[:, idx] .+ α .* pts[:, idx + 1]
end


"""Estimate curvature at each sampled path point from neighboring samples."""
function rpo_curvature_from_samples(samples, s_vals)
    pts = Matrix{Float64}(samples)
    s = Vector{Float64}(s_vals)

    n = size(pts, 2)
    κ = zeros(n)

    n < 3 && return κ

    @inbounds for j in 2:(n - 1)
        ds1 = s[j] - s[j - 1]
        ds2 = s[j + 1] - s[j]
        ds = s[j + 1] - s[j - 1]

        if ds1 <= eps(Float64) || ds2 <= eps(Float64) || ds <= eps(Float64)
            κ[j] = 0.0
            continue
        end

        r′ = (pts[:, j + 1] - pts[:, j - 1]) / ds

        r″ = 2.0 .* (
            (pts[:, j + 1] - pts[:, j]) / ds2 -
            (pts[:, j] - pts[:, j - 1]) / ds1
        ) / ds

        rpn = norm(r′)

        if rpn > eps(Float64)
            κ[j] = norm(cross(r′, r″)) / rpn^3
        else
            κ[j] = 0.0
        end
    end

    κ[1] = κ[2]
    κ[end] = κ[end - 1]

    return κ
end


"""
Distance available for the stopping-distance limit of the retimer.

In `:manuscript` mode this is Sec. III.E's d_avail = max(0, c - d_safe), where c
is the clearance to the keep-out surface. The legacy retimer uses the raw
distance to the nearest station point.
"""
@inline function rpo_retime_available_distance(cfg::RPOPSOConfig, clearance::Real, distance::Real, safe_distance_m::Real)
    cfg.hypr_mode === :manuscript && return max(0.0, Float64(clearance) - Float64(safe_distance_m))
    return max(0.0, Float64(distance))
end

"""
Pointwise retiming speed limit of Sec. III.E at one sample: the clearance-limited
speed from the stopping-distance inequality and the curvature-limited speed,
scaled by the safety factor k_v. In `:manuscript` mode the global cap enters
the minimum before scaling, v = k_v min(v_clear, v_curv, v_max); the legacy
retimer applies the cap after scaling.
"""
@inline function rpo_retime_pointwise_speed(cfg::RPOPSOConfig, d_avail::Real, curvature::Real)
    amax = Float64(cfg.retime_a_max_mps2)
    reaction_time = Float64(cfg.retime_reaction_time_s)
    d = Float64(d_avail)
    κ = Float64(curvature)
    v_clear = if d <= 0.0 || amax <= 0.0
        0.0
    else
        -amax * reaction_time + sqrt((amax * reaction_time)^2 + 2.0 * amax * d)
    end
    v_clear = max(0.0, v_clear)
    v_curve = κ <= 0.0 || amax <= 0.0 ? Inf : sqrt(amax / κ)
    if cfg.hypr_mode === :manuscript
        return Float64(cfg.retime_speed_scale) * min(v_clear, v_curve, Float64(cfg.retime_max_speed_mps))
    end
    v = Float64(cfg.retime_speed_scale) * min(v_clear, v_curve)
    max_speed = Float64(cfg.retime_max_speed_mps)
    isfinite(max_speed) && (v = min(v, max_speed))
    return v
end

"""Collision-sampling spacing the retimer uses: the planner's density in `:manuscript` mode, `sample_ds_m` otherwise."""
function rpo_retime_sampling_ds_m(cfg::RPOPSOConfig, safe_distance_m::Real)
    cfg.hypr_mode === :manuscript && return rpo_hypr_sampling_density_m(cfg, safe_distance_m)
    return cfg.sample_ds_m
end

"""Retiming path samples into position, velocity, acceleration, and time references."""
function rpo_retime_path(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    fallback_speed_mps::Real=1.0e-3,
    duplicate_tol_m::Real=1.0e-10,
)
    if cfg.retime_accel_limit_enable
        ref = rpo_retimed_reference(
            points,
            geometry,
            cfg;
            safe_distance_m=safe_distance_m,
            fallback_speed_mps=fallback_speed_mps,
            duplicate_tol_m=duplicate_tol_m,
        )
        return ref.r_rtn, ref.s_m, ref.speed_mps
    end
    raw_samples = rpo_sample_path(
        points,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=rpo_retime_sampling_ds_m(cfg, safe_distance_m),
        curve_type=cfg.curve_type,
    )
    samples = rpo_remove_near_duplicate_samples(raw_samples; tol=duplicate_tol_m)

    s_samples = rpo_arc_length_params(samples)
    total = s_samples[end]

    if total <= eps(Float64)
        @warn "RPO retiming received a zero-length path. Returning a single-point reference."
        return samples[:, 1:1], [0.0], [0.0]
    end

    κ = rpo_curvature_from_samples(samples, s_samples)

    n = length(s_samples)
    v_max = zeros(n)
    geometry_distance = zeros(n)

    max_speed = Float64(cfg.retime_max_speed_mps)
    cfg_min_speed = Float64(cfg.retime_min_speed_mps)

    # This is the minimum speed used when the path has locally infeasible geometry distance.
    # It prevents the batch run from crashing, but it does not make the path collision-free.
    fallback_speed = max(Float64(fallback_speed_mps), eps(Float64))

    if isfinite(max_speed)
        fallback_speed = min(fallback_speed, max_speed)
    end

    fallback_speed = max(fallback_speed, eps(Float64))

    infeasible_idxs = Int[]

    @inbounds for j in eachindex(s_samples)
        station = rpo_clearance_to_station(samples[:, j], geometry)
        geometry_distance[j] = station.distance

        d_avail = rpo_retime_available_distance(cfg, station.clearance, station.distance, safe_distance_m)

        v = rpo_retime_pointwise_speed(cfg, d_avail, κ[j])

        # If the geometry distance model says v = 0, warn and use a tiny fallback speed
        # instead of allowing the retimer to stall.
        if v <= 0.0
            push!(infeasible_idxs, j)
            v = fallback_speed
        end

        # Preserve user-configured minimum speed if it is positive.
        v_max[j] = max(v, cfg_min_speed)
    end

    if !isempty(infeasible_idxs)
        k = first(infeasible_idxs)

        @warn "RPO retiming encountered infeasible zero-speed samples. Continuing with fallback speed. The path likely intersects the RPO geometry." (
            count = length(infeasible_idxs),
            first_idx = k,
            first_s_m = s_samples[k],
            total_s_m = total,
            first_position_m = samples[:, k],
            geometry_distance_m = geometry_distance[k],
            curvature_1pm = κ[k],
            fallback_speed_mps = fallback_speed,
        )
    end

    if maximum(v_max) <= 0.0
        @warn "RPO retiming found no feasible positive speeds. Using constant fallback speed for entire path." fallback_speed_mps=fallback_speed
        fill!(v_max, fallback_speed)
    end

    r_hist = Vector{Vector{Float64}}()
    s_hist = Float64[]
    v_hist = Float64[]

    s = 0.0
    steps = 0

    while true
        idx = clamp(searchsortedlast(s_samples, s), 1, length(s_samples))

        v = v_max[idx]

        # Defensive fallback. This should rarely trigger because v_max was already repaired.
        if v <= 0.0 || !isfinite(v)
            @warn "RPO retiming hit invalid local speed during propagation. Replacing with fallback speed." (
                idx = idx,
                s_m = s,
                total_s_m = total,
                v_mps = v,
                fallback_speed_mps = fallback_speed,
                geometry_distance_m = geometry_distance[idx],
                curvature_1pm = κ[idx],
            )

            v = fallback_speed
        end

        push!(s_hist, s)
        push!(v_hist, v)
        push!(r_hist, rpo_interpolate_along_path(samples, s_samples, s))

        s >= total - 1.0e-9 && break

        steps += 1

        if steps > cfg.retime_max_steps
            @warn "RPO retiming exceeded maximum step count. Forcing final endpoint into returned trajectory." (
                steps = steps,
                max_steps = cfg.retime_max_steps,
                current_s_m = s,
                total_s_m = total,
            )

            if s < total
                push!(s_hist, total)
                push!(v_hist, fallback_speed)
                push!(r_hist, copy(samples[:, end]))
            end

            break
        end

        ds = v * Float64(cfg.retime_dt_s)

        if ds <= eps(Float64) || !isfinite(ds)
            @warn "RPO retiming computed invalid arc-length step. Using fallback step." (
                s_m = s,
                total_s_m = total,
                v_mps = v,
                dt_s = cfg.retime_dt_s,
                fallback_speed_mps = fallback_speed,
            )

            ds = fallback_speed * Float64(cfg.retime_dt_s)
        end

        s = min(total, s + ds)
    end

    mat = zeros(3, length(r_hist))

    for (j, r) in enumerate(r_hist)
        mat[:, j] .= r
    end

    return mat, s_hist, v_hist
end


# Acceleration-limited retiming (`retime_accel_limit_enable`), a documented
# addition to Sec. III.E: the pointwise limits are kept, and forward and
# backward passes then bound the tangential acceleration and fix the boundary
# speeds, so the reference starts at the chaser's speed and ends at rest.

const _RPO_GAUSS5_NODES = (-0.906179845938664, -0.5384693101056831, 0.0, 0.5384693101056831, 0.906179845938664)
const _RPO_GAUSS5_WEIGHTS = (0.23692688505618908, 0.47862867049936647, 0.5688888888888889, 0.47862867049936647, 0.23692688505618908)

"""Curve being retimed: Bezier control points with their first and second hodographs, or polyline vertices."""
struct RPORetimeCurve
    curve_type::Symbol
    ctrl::Matrix{Float64}
    d1::Matrix{Float64}
    d2::Matrix{Float64}
    work::Matrix{Float64}
end

"""Build the retiming curve data for a control polygon; each evaluation reuses the curve's own work buffer."""
function RPORetimeCurve(points, curve_type::Symbol)
    ctrl = Matrix{Float64}(points)
    m = size(ctrl, 2)
    d1 = m >= 2 ? (m - 1) .* (ctrl[:, 2:end] .- ctrl[:, 1:(end - 1)]) : zeros(3, 0)
    d2 = m >= 3 ? (m - 2) .* (d1[:, 2:end] .- d1[:, 1:(end - 1)]) : zeros(3, 0)
    return RPORetimeCurve(curve_type, ctrl, d1, d2, zeros(3, max(m, 1)))
end

"""Evaluate a 3-D Bezier polygon at `u` by de Casteljau's recursion without allocating."""
@inline function _rpo_bezier_eval3!(work::Matrix{Float64}, ctrl::Matrix{Float64}, u::Float64)
    m = size(ctrl, 2)
    m == 0 && return SVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for j in 1:m
        work[1, j] = ctrl[1, j]
        work[2, j] = ctrl[2, j]
        work[3, j] = ctrl[3, j]
    end
    a = 1.0 - u
    @inbounds for r in 1:(m - 1)
        for j in 1:(m - r)
            work[1, j] = a * work[1, j] + u * work[1, j + 1]
            work[2, j] = a * work[2, j] + u * work[2, j + 1]
            work[3, j] = a * work[3, j] + u * work[3, j + 1]
        end
    end
    return @inbounds SVector{3, Float64}(work[1, 1], work[2, 1], work[3, 1])
end

"""Point on the retimed Bezier curve at parameter `u`."""
@inline _rpo_curve_point(curve::RPORetimeCurve, u::Float64) = _rpo_bezier_eval3!(curve.work, curve.ctrl, u)
"""First derivative dr/du of the retimed Bezier curve."""
@inline _rpo_curve_d1(curve::RPORetimeCurve, u::Float64) = _rpo_bezier_eval3!(curve.work, curve.d1, u)
"""Second derivative d²r/du² of the retimed Bezier curve."""
@inline _rpo_curve_d2(curve::RPORetimeCurve, u::Float64) = _rpo_bezier_eval3!(curve.work, curve.d2, u)

"""Column indices kept after dropping consecutive samples closer than `tol`; the last sample is always kept."""
function rpo_near_duplicate_keep_indices(samples, tol::Real)
    n = size(samples, 2)
    keep = Int[]
    n == 0 && return keep
    push!(keep, 1)
    @inbounds for j in 2:n
        i = keep[end]
        dx = samples[1, j] - samples[1, i]
        dy = samples[2, j] - samples[2, i]
        dz = samples[3, j] - samples[3, i]
        sqrt(dx * dx + dy * dy + dz * dz) > Float64(tol) && push!(keep, j)
    end
    if keep[end] != n
        length(keep) == 1 ? push!(keep, n) : (keep[end] = n)
    end
    return keep
end

"""
    rpo_accel_limited_speeds(s, v_point, a_max; v_start=0.0, v_end=0.0)

Limit pointwise speeds `v_point` at arc lengths `s` so the tangential
acceleration never exceeds `a_max`: the profile starts at `v_start`, ends at
`v_end`, and a forward and a backward pass bound the change of v² between
neighbouring samples by 2 a_max Δs.
"""
function rpo_accel_limited_speeds(s, v_point, a_max::Real; v_start::Real=0.0, v_end::Real=0.0)
    n = length(s)
    v = Vector{Float64}(v_point)
    n == 0 && return v
    a = Float64(a_max)
    v[1] = Float64(v_start)
    n == 1 && return v
    v[n] = Float64(v_end)
    @inbounds for j in 1:(n - 1)
        ds = max(s[j + 1] - s[j], 0.0)
        v[j + 1] = min(v[j + 1], sqrt(v[j]^2 + 2.0 * a * ds))
    end
    @inbounds for j in (n - 1):-1:1
        ds = max(s[j + 1] - s[j], 0.0)
        v[j] = min(v[j], sqrt(v[j + 1]^2 + 2.0 * a * ds))
    end
    return v
end

"""
    rpo_retime_profile(curve, samples, params, clearances, geometry, cfg; safe_distance_m=0.0)

Acceleration-limited speed profile over the samples of one path.

Pointwise limits follow Sec. III.E (`rpo_retime_pointwise_speed`); forward and
backward passes then bound the tangential acceleration by `retime_a_max_mps2`,
starting at `retime_initial_speed_mps` and ending at rest. Within a sample
interval v² is linear in arc length (constant acceleration), so crossing it
takes 2Δs/(v_i + v_{i+1}). For a Bezier curve, arc length comes from
Gauss-Legendre quadrature of |r'(u)| and curvature from the curve's
derivatives at the samples; a polyline uses chord lengths and the three-point
estimate. Clearances that are `NaN` are computed here. Interior samples whose
pointwise limit is zero (inside the safety distance) get `fallback_speed_mps`,
as in the legacy retimer, so the reference cannot stall there.
"""
function rpo_retime_profile(
    curve::RPORetimeCurve,
    samples,
    params,
    clearances,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    fallback_speed_mps::Real=1.0e-3,
    duplicate_tol_m::Real=1.0e-10,
    warn::Bool=true,
)
    bezier = curve.curve_type == :bezier && !isempty(params)
    keep = rpo_near_duplicate_keep_indices(samples, duplicate_tol_m)
    pts = Matrix{Float64}(samples[:, keep])
    u = bezier ? Vector{Float64}(params[keep]) : Float64[]
    clear = Vector{Float64}(clearances[keep])
    if size(pts, 2) == 2
        # One interval cannot go from rest to rest at constant acceleration; split it.
        if bezier
            u_mid = 0.5 * (u[1] + u[2])
            mid = _rpo_curve_point(curve, u_mid)
            insert!(u, 2, u_mid)
        else
            mid = 0.5 .* (SVector{3, Float64}(pts[:, 1]) .+ SVector{3, Float64}(pts[:, 2]))
        end
        pts = hcat(pts[:, 1], collect(mid), pts[:, 2])
        insert!(clear, 2, NaN)
    end
    n = size(pts, 2)
    s = zeros(n)
    rprime = zeros(n)
    κ = zeros(n)
    if bezier
        @inbounds for j in 1:n
            d1 = _rpo_curve_d1(curve, u[j])
            rprime[j] = norm(d1)
            if size(curve.d2, 2) > 0 && rprime[j] > 1.0e-12
                d2 = _rpo_curve_d2(curve, u[j])
                κ[j] = norm(cross(d1, d2)) / rprime[j]^3
            end
        end
        @inbounds for j in 1:(n - 1)
            half = 0.5 * (u[j + 1] - u[j])
            mid = 0.5 * (u[j + 1] + u[j])
            acc = 0.0
            for q in 1:5
                acc += _RPO_GAUSS5_WEIGHTS[q] * norm(_rpo_curve_d1(curve, mid + half * _RPO_GAUSS5_NODES[q]))
            end
            s[j + 1] = s[j] + half * acc
        end
    else
        @inbounds for j in 2:n
            dx = pts[1, j] - pts[1, j - 1]
            dy = pts[2, j] - pts[2, j - 1]
            dz = pts[3, j] - pts[3, j - 1]
            s[j] = s[j - 1] + sqrt(dx * dx + dy * dy + dz * dz)
        end
        @inbounds for j in 2:(n - 1)
            ds1 = s[j] - s[j - 1]
            ds2 = s[j + 1] - s[j]
            ds = s[j + 1] - s[j - 1]
            (ds1 <= eps(Float64) || ds2 <= eps(Float64)) && continue
            r1 = SVector{3, Float64}(pts[1, j + 1] - pts[1, j - 1], pts[2, j + 1] - pts[2, j - 1], pts[3, j + 1] - pts[3, j - 1]) / ds
            fwd = SVector{3, Float64}(pts[1, j + 1] - pts[1, j], pts[2, j + 1] - pts[2, j], pts[3, j + 1] - pts[3, j]) / ds2
            bwd = SVector{3, Float64}(pts[1, j] - pts[1, j - 1], pts[2, j] - pts[2, j - 1], pts[3, j] - pts[3, j - 1]) / ds1
            r2 = 2.0 .* (fwd - bwd) ./ ds
            rn = norm(r1)
            rn > eps(Float64) && (κ[j] = norm(cross(r1, r2)) / rn^3)
        end
        n >= 3 && (κ[1] = κ[2]; κ[n] = κ[n - 1])
    end

    body_margin = geometry.station.keepout_radius_m + maximum(geometry.chaser.half_extents_body)
    fallback = max(Float64(fallback_speed_mps), eps(Float64))
    isfinite(cfg.retime_max_speed_mps) && (fallback = max(min(fallback, cfg.retime_max_speed_mps), eps(Float64)))
    v_point = zeros(n)
    fallback_count = 0
    @inbounds for j in 1:n
        if isnan(clear[j])
            clear[j] = rpo_clearance_distance_to_station(SVector{3, Float64}(pts[1, j], pts[2, j], pts[3, j]), geometry)
        end
        d_avail = rpo_retime_available_distance(cfg, clear[j], clear[j] + body_margin, safe_distance_m)
        v = rpo_retime_pointwise_speed(cfg, d_avail, κ[j])
        if 1 < j < n
            if v <= 0.0
                fallback_count += 1
                v = fallback
            end
            v = max(v, cfg.retime_min_speed_mps)
        end
        v_point[j] = v
    end
    if warn && fallback_count > 0
        @warn "RPO retiming encountered infeasible zero-speed samples. Continuing with fallback speed. The path likely intersects the RPO geometry." count=fallback_count fallback_speed_mps=fallback
    end

    v_start = Float64(cfg.retime_initial_speed_mps)
    v = rpo_accel_limited_speeds(s, v_point, cfg.retime_a_max_mps2; v_start=v_start, v_end=0.0)
    if warn && v[1] < v_start - 1.0e-12
        @warn "RPO retiming cannot keep the initial speed within the acceleration limit; the reference starts slower." initial_speed_mps=v_start reference_start_mps=v[1]
    end

    t = zeros(n)
    a_seg = zeros(max(n - 1, 0))
    @inbounds for j in 1:(n - 1)
        ds = s[j + 1] - s[j]
        vs = v[j] + v[j + 1]
        dt = ds <= 0.0 ? 0.0 : (vs > 0.0 ? 2.0 * ds / vs : Inf)
        t[j + 1] = t[j] + dt
        a_seg[j] = ds > 0.0 ? (v[j + 1]^2 - v[j]^2) / (2.0 * ds) : 0.0
    end
    return (
        curve=curve,
        bezier=bezier,
        samples=pts,
        params=u,
        s=s,
        t=t,
        v=v,
        v_point=v_point,
        curvature=κ,
        clearance=clear,
        rprime=rprime,
        a_seg=a_seg,
        fallback_count=fallback_count,
        length_m=s[end],
        duration_s=t[end],
    )
end

"""Curve parameter at arc length `sq` inside profile interval `j` (cubic Hermite with du/ds = 1/|r'(u)| at the samples)."""
@inline function _rpo_profile_param(profile, j::Int, sq::Float64)
    s = profile.s
    u = profile.params
    h = s[j + 1] - s[j]
    h <= 0.0 && return u[j]
    τ = clamp((sq - s[j]) / h, 0.0, 1.0)
    du = u[j + 1] - u[j]
    rp0 = profile.rprime[j]
    rp1 = profile.rprime[j + 1]
    m0 = rp0 > 1.0e-12 ? h / rp0 : du
    m1 = rp1 > 1.0e-12 ? h / rp1 : du
    τ2 = τ * τ
    τ3 = τ2 * τ
    uq = (2.0 * τ3 - 3.0 * τ2 + 1.0) * u[j] + (τ3 - 2.0 * τ2 + τ) * m0 + (-2.0 * τ3 + 3.0 * τ2) * u[j + 1] + (τ3 - τ2) * m1
    return clamp(uq, u[j], u[j + 1])
end

"""Position and unit tangent at arc length `sq` inside profile interval `j`."""
@inline function _rpo_profile_point_tangent(profile, j::Int, sq::Float64)
    if profile.bezier
        uq = _rpo_profile_param(profile, j, sq)
        p = _rpo_curve_point(profile.curve, uq)
        d = _rpo_curve_d1(profile.curve, uq)
        dn = norm(d)
        return p, dn > 1.0e-12 ? d / dn : SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    pts = profile.samples
    s = profile.s
    a = @inbounds SVector{3, Float64}(pts[1, j], pts[2, j], pts[3, j])
    b = @inbounds SVector{3, Float64}(pts[1, j + 1], pts[2, j + 1], pts[3, j + 1])
    h = s[j + 1] - s[j]
    h <= 0.0 && return a, SVector{3, Float64}(0.0, 0.0, 0.0)
    α = clamp((sq - s[j]) / h, 0.0, 1.0)
    return (1.0 - α) * a + α * b, (b - a) / h
end

"""Profile interval, arc length, speed and tangential acceleration at time `tq`, scanning forward from `j`."""
@inline function _rpo_profile_state_at_time(profile, tq::Float64, j::Int)
    t = profile.t
    n = length(t)
    if tq >= t[n]
        return n - 1, profile.s[n], profile.v[n], 0.0
    end
    @inbounds while j < n - 1 && t[j + 1] <= tq
        j += 1
    end
    τ = tq - t[j]
    a = profile.a_seg[j]
    sq = min(profile.s[j] + profile.v[j] * τ + 0.5 * a * τ * τ, profile.s[j + 1])
    return j, sq, max(0.0, profile.v[j] + a * τ), a
end

"""Number of fixed steps that brings the reference to the goal: t_K = K Δt ≥ duration."""
@inline function _rpo_profile_step_count(duration_s::Float64, dt::Float64)
    duration_s <= 0.0 && return 0
    return max(1, ceil(Int, duration_s / dt - 1.0e-9))
end

"""
    rpo_retimed_reference_from_profile(profile, dt)

Sample an acceleration-limited profile at the fixed step `dt`: t_k = k dt for
k = 0..K with t_K the first step at or after arrival. Positions lie on the
curve, velocities are speed times the unit tangent, and the last sample is
the goal at rest.
"""
function rpo_retimed_reference_from_profile(profile, dt::Real)
    step = Float64(dt)
    n = length(profile.s)
    K = _rpo_profile_step_count(profile.duration_s, step)
    t_s = collect(0:K) .* step
    r = zeros(3, K + 1)
    v = zeros(3, K + 1)
    s_ref = zeros(K + 1)
    speed = zeros(K + 1)
    accel = zeros(K + 1)
    j = 1
    @inbounds for k in 0:K
        if n == 1
            r[:, k + 1] .= profile.samples[:, 1]
            continue
        end
        j, sq, vq, a = _rpo_profile_state_at_time(profile, t_s[k + 1], j)
        p, tangent = _rpo_profile_point_tangent(profile, j, sq)
        r[1, k + 1] = p[1]
        r[2, k + 1] = p[2]
        r[3, k + 1] = p[3]
        v[1, k + 1] = vq * tangent[1]
        v[2, k + 1] = vq * tangent[2]
        v[3, k + 1] = vq * tangent[3]
        s_ref[k + 1] = sq
        speed[k + 1] = vq
        accel[k + 1] = a
    end
    r[:, end] .= profile.samples[:, end]
    v[:, end] .= 0.0
    s_ref[end] = profile.s[end]
    speed[end] = 0.0
    accel[end] = 0.0
    return (t_s=t_s, r_rtn=r, v_rtn=v, s_m=s_ref, speed_mps=speed, accel_tangential_mps2=accel, profile=profile)
end

"""
    rpo_retimed_reference(points, geometry, cfg; safe_distance_m=0.0)

Acceleration-limited reference for a control polygon: sample the path as the
planner does, build `rpo_retime_profile`, and sample it at `retime_dt_s`.
"""
function rpo_retimed_reference(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    fallback_speed_mps::Real=1.0e-3,
    duplicate_tol_m::Real=1.0e-10,
    warn::Bool=true,
)
    samples, params, clearances = rpo_sample_path_with_params(
        points,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=rpo_retime_sampling_ds_m(cfg, safe_distance_m),
        curve_type=cfg.curve_type,
    )
    profile = rpo_retime_profile(
        RPORetimeCurve(points, cfg.curve_type),
        samples,
        params,
        clearances,
        geometry,
        cfg;
        safe_distance_m=safe_distance_m,
        fallback_speed_mps=fallback_speed_mps,
        duplicate_tol_m=duplicate_tol_m,
        warn=warn,
    )
    return rpo_retimed_reference_from_profile(profile, cfg.retime_dt_s)
end

"""
    rpo_reference_accel_demand(ref, mean_motion)

Acceleration an acceleration-limited reference (`rpo_retimed_reference`) asks
of the actuators at each sample: the tangential acceleration along the curve,
the centripetal term v²κ toward the centre of curvature, and the HCW terms,
u = r̈ - [3n² x + 2n ẏ, -2n ẋ, -n² z]. Returns per-sample magnitudes
`tangential_mps2`, `centripetal_mps2` and `hcw_mps2`, and `u_rtn` (3 x K).
A polyline has no curvature between its vertices, so its centripetal term is
zero here.
"""
function rpo_reference_accel_demand(ref, mean_motion::Real)
    profile = ref.profile
    K = length(ref.t_s)
    n = Float64(mean_motion)
    u_rtn = zeros(3, K)
    tangential = zeros(K)
    centripetal = zeros(K)
    hcw = zeros(K)
    ns = length(profile.s)
    j = 1
    @inbounds for k in 1:K
        speed = ref.speed_mps[k]
        a_t = ref.accel_tangential_mps2[k]
        x, y, z = ref.r_rtn[1, k], ref.r_rtn[2, k], ref.r_rtn[3, k]
        vx, vy, vz = ref.v_rtn[1, k], ref.v_rtn[2, k], ref.v_rtn[3, k]
        tangent = SVector{3, Float64}(0.0, 0.0, 0.0)
        κvec = SVector{3, Float64}(0.0, 0.0, 0.0)
        if ns >= 2
            sq = ref.s_m[k]
            while j < ns - 1 && profile.s[j + 1] < sq
                j += 1
            end
            if profile.bezier
                uq = _rpo_profile_param(profile, j, sq)
                d1 = _rpo_curve_d1(profile.curve, uq)
                d1n = norm(d1)
                if d1n > 1.0e-12
                    tangent = d1 / d1n
                    if size(profile.curve.d2, 2) > 0
                        d2 = _rpo_curve_d2(profile.curve, uq)
                        κvec = (d2 - dot(d2, tangent) * tangent) / d1n^2
                    end
                end
            else
                _, tangent = _rpo_profile_point_tangent(profile, j, sq)
            end
        end
        acc = a_t * tangent + speed^2 * κvec
        g = SVector{3, Float64}(3.0 * n * n * x + 2.0 * n * vy, -2.0 * n * vx, -n * n * z)
        u = acc - g
        u_rtn[1, k] = u[1]
        u_rtn[2, k] = u[2]
        u_rtn[3, k] = u[3]
        tangential[k] = abs(a_t)
        centripetal[k] = speed^2 * norm(κvec)
        hcw[k] = norm(g)
    end
    return (tangential_mps2=tangential, centripetal_mps2=centripetal, hcw_mps2=hcw, u_rtn=u_rtn)
end
