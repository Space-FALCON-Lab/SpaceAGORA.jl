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


"""Retiming path samples into position, velocity, acceleration, and time references."""
function rpo_retime_path(
    points,
    geometry,
    cfg::RPOPSOConfig;
    safe_distance_m::Real=0.0,
    fallback_speed_mps::Real=1.0e-3,
    duplicate_tol_m::Real=1.0e-10,
)
    raw_samples = rpo_sample_path(
        points,
        cfg,
        geometry;
        safe_distance_m=safe_distance_m,
        base_ds_m=cfg.sample_ds_m,
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

    amax = Float64(cfg.retime_a_max_mps2)
    reaction_time = Float64(cfg.retime_reaction_time_s)
    speed_scale = Float64(cfg.retime_speed_scale)
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
        geometry_distance[j] = rpo_clearance_to_station(samples[:, j], geometry).distance

        d_avail = max(0.0, geometry_distance[j])

        v_clear = if d_avail <= 0.0 || amax <= 0.0
            0.0
        else
            -amax * reaction_time + sqrt((amax * reaction_time)^2 + 2.0 * amax * d_avail)
        end

        v_clear = max(0.0, v_clear)

        v_curve = if κ[j] <= 0.0 || amax <= 0.0
            Inf
        else
            sqrt(amax / κ[j])
        end

        v = speed_scale * min(v_clear, v_curve)

        if isfinite(max_speed)
            v = min(v, max_speed)
        end

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
