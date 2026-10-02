# Retiming calculations from explicit inputs; configured policy belongs to HYPR.
# Definitions retain their GuidanceHooks identity and existing geometry contract.

"""
    rpo_retime_samples(raw_samples, geometry; max_speed_mps, min_speed_mps,
        dt_s, max_steps, available_distance, pointwise_speed, ...)

Advance along already sampled positions using supplied speed policies. Distances
are metres, speeds m/s and `dt_s` seconds. `available_distance(clearance, distance,
safe_distance)` returns available metres; `pointwise_speed(distance, curvature)`
returns m/s with curvature in 1/m. Policies must preserve inputs and avoid hidden
random draws. Existing fallback-speed and maximum-step behavior are retained;
fallback output is not a collision-free or dynamically feasible certificate.
"""
function rpo_retime_samples(
    raw_samples,
    geometry;
    max_speed_mps,
    min_speed_mps,
    dt_s,
    max_steps,
    available_distance,
    pointwise_speed,
    safe_distance_m::Real=0.0,
    fallback_speed_mps::Real=1.0e-3,
    duplicate_tol_m::Real=1.0e-10,
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

    max_speed = Float64(max_speed_mps)
    cfg_min_speed = Float64(min_speed_mps)

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

        d_avail = available_distance(station.clearance, station.distance, safe_distance_m)

        v = pointwise_speed(d_avail, κ[j])

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

        if steps > max_steps
            @warn "RPO retiming exceeded maximum step count. Forcing final endpoint into returned trajectory." (
                steps = steps,
                max_steps = max_steps,
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

        ds = v * Float64(dt_s)

        if ds <= eps(Float64) || !isfinite(ds)
            @warn "RPO retiming computed invalid arc-length step. Using fallback step." (
                s_m = s,
                total_s_m = total,
                v_mps = v,
                dt_s = dt_s,
                fallback_speed_mps = fallback_speed,
            )

            ds = fallback_speed * Float64(dt_s)
        end

        s = min(total, s + ds)
    end

    mat = zeros(3, length(r_hist))

    for (j, r) in enumerate(r_hist)
        mat[:, j] .= r
    end

    return mat, s_hist, v_hist
end

"""
    rpo_retime_profile(curve, samples, params, clearances, geometry;
        max_speed_mps, min_speed_mps, initial_speed_mps, a_max_mps2,
        available_distance, pointwise_speed, ...)

Construct an acceleration-limited profile from supplied samples and policies.
Units and policy callback arguments match `rpo_retime_samples`; acceleration is
m/s². Preserve the existing duplicate removal, two-point split, curve quadrature,
clearance queries, forward/backward limits, terminal rest and fallback behavior.
Missing (`NaN`) clearances are evaluated through the existing geometry owner.
Sampling policy and configured defaults belong to the caller. This internal
boundary makes no broader numeric-type or physical-feasibility guarantee.
"""
function rpo_retime_profile(
    curve::RPORetimeCurve,
    samples,
    params,
    clearances,
    geometry;
    max_speed_mps,
    min_speed_mps,
    initial_speed_mps,
    a_max_mps2,
    available_distance,
    pointwise_speed,
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
    isfinite(max_speed_mps) && (fallback = max(min(fallback, max_speed_mps), eps(Float64)))
    v_point = zeros(n)
    fallback_count = 0
    @inbounds for j in 1:n
        if isnan(clear[j])
            clear[j] = rpo_clearance_distance_to_station(SVector{3, Float64}(pts[1, j], pts[2, j], pts[3, j]), geometry)
        end
        d_avail = available_distance(clear[j], clear[j] + body_margin, safe_distance_m)
        v = pointwise_speed(d_avail, κ[j])
        if 1 < j < n
            if v <= 0.0
                fallback_count += 1
                v = fallback
            end
            v = max(v, min_speed_mps)
        end
        v_point[j] = v
    end
    if warn && fallback_count > 0
        @warn "RPO retiming encountered infeasible zero-speed samples. Continuing with fallback speed. The path likely intersects the RPO geometry." count=fallback_count fallback_speed_mps=fallback
    end

    v_start = Float64(initial_speed_mps)
    v = rpo_accel_limited_speeds(s, v_point, a_max_mps2; v_start=v_start, v_end=0.0)
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
