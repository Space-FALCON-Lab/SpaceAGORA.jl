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
function rpo_remove_near_duplicate_samples(
    points;
    tol::Real=1.0e-10,
    warn_removed::Bool=true,
)
    pts = Matrix{Float64}(points)
    n = size(pts, 2)

    n <= 1 && return pts

    keep = Int[1]

    @inbounds for j in 2:n
        if norm(pts[:, j] - pts[:, keep[end]]) > Float64(tol)
            push!(keep, j)
        end
    end

    if warn_removed && length(keep) < n
        @warn "RPO retiming removed near-duplicate path samples." removed=(n - length(keep)) kept=length(keep)
    end

    return pts[:, keep]
end


"""Interpolate a point at a requested arc-length coordinate into a caller-provided buffer."""
function rpo_interpolate_along_path!(out, pts, s, s_query::Real)
    sq = Float64(s_query)

    n = length(s)

    if n == 1 || sq <= s[1]
        @inbounds for axis in 1:3
            out[axis] = pts[axis, 1]
        end
        return out
    elseif sq >= s[end]
        @inbounds for axis in 1:3
            out[axis] = pts[axis, end]
        end
        return out
    end

    idx = clamp(searchsortedlast(s, sq), 1, n - 1)

    # Move forward if this segment has zero length.
    while idx < n && s[idx + 1] - s[idx] <= eps(Float64)
        idx += 1
    end

    if idx >= n
        @inbounds for axis in 1:3
            out[axis] = pts[axis, end]
        end
        return out
    end

    denom = s[idx + 1] - s[idx]

    if denom <= eps(Float64)
        @inbounds for axis in 1:3
            out[axis] = pts[axis, idx]
        end
        return out
    end

    α = clamp((sq - s[idx]) / denom, 0.0, 1.0)
    @inbounds for axis in 1:3
        out[axis] = (1.0 - α) * pts[axis, idx] + α * pts[axis, idx + 1]
    end
    return out
end


"""Interpolate a point at a requested arc-length coordinate."""
function rpo_interpolate_along_path(points, s_vals, s_query::Real)
    pts = Matrix{Float64}(points)
    s = Vector{Float64}(s_vals)
    out = zeros(3)
    return rpo_interpolate_along_path!(out, pts, s, s_query)
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

        rx = (pts[1, j + 1] - pts[1, j - 1]) / ds
        ry = (pts[2, j + 1] - pts[2, j - 1]) / ds
        rz = (pts[3, j + 1] - pts[3, j - 1]) / ds
        rxx = 2.0 * (
            (pts[1, j + 1] - pts[1, j]) / ds2 -
            (pts[1, j] - pts[1, j - 1]) / ds1
        ) / ds
        ryy = 2.0 * (
            (pts[2, j + 1] - pts[2, j]) / ds2 -
            (pts[2, j] - pts[2, j - 1]) / ds1
        ) / ds
        rzz = 2.0 * (
            (pts[3, j + 1] - pts[3, j]) / ds2 -
            (pts[3, j] - pts[3, j - 1]) / ds1
        ) / ds
        tangent_norm = sqrt(rx * rx + ry * ry + rz * rz)

        if tangent_norm > eps(Float64)
            cross_x = ry * rzz - rz * ryy
            cross_y = rz * rxx - rx * rzz
            cross_z = rx * ryy - ry * rxx
            cross_norm = sqrt(
                cross_x * cross_x + cross_y * cross_y + cross_z * cross_z,
            )
            κ[j] = cross_norm / tangent_norm^3
        else
            κ[j] = 0.0
        end
    end

    κ[1] = κ[2]
    κ[end] = κ[end - 1]

    return κ
end


"""Compact retiming data used to materialize or stream an RPO reference."""
struct RPORetimingProfile
    samples::Matrix{Float64}
    s_samples::Vector{Float64}
    s_ref::Vector{Float64}
    speed_ref::Vector{Float64}
    dt_s::Float64
end


"""Largest finite-difference HCW feedforward command in a retiming profile."""
function rpo_hcw_feedforward_max_acceleration(
    profile::RPORetimingProfile,
    mean_motion_radps::Real,
)
    n_reference = length(profile.s_ref)
    n_steps = max(n_reference - 1, 0)
    n_steps == 0 && return 0.0

    dt = profile.dt_s
    mean_motion = Float64(mean_motion_radps)
    current_position = zeros(3)
    next_position = zeros(3)
    following_position = zeros(3)
    current_velocity = zeros(3)
    next_velocity = zeros(3)
    rpo_interpolate_along_path!(
        current_position, profile.samples, profile.s_samples, profile.s_ref[1],
    )
    rpo_interpolate_along_path!(
        next_position, profile.samples, profile.s_samples, profile.s_ref[2],
    )
    # The reference departs from rest; interior velocities use forward differences.

    maximum_acceleration = 0.0
    @inbounds for step_index in 1:n_steps
        if step_index < n_steps
            rpo_interpolate_along_path!(
                following_position,
                profile.samples,
                profile.s_samples,
                profile.s_ref[step_index + 2],
            )
            for axis in 1:3
                next_velocity[axis] =
                    (following_position[axis] - next_position[axis]) / dt
            end
        else
            next_velocity .= 0.0
        end

        ax = (next_velocity[1] - current_velocity[1]) / dt -
            (3.0 * mean_motion^2 * current_position[1] +
             2.0 * mean_motion * current_velocity[2])
        ay = (next_velocity[2] - current_velocity[2]) / dt +
            2.0 * mean_motion * current_velocity[1]
        az = (next_velocity[3] - current_velocity[3]) / dt +
            mean_motion^2 * current_position[3]
        maximum_acceleration = max(
            maximum_acceleration, sqrt(ax * ax + ay * ay + az * az),
        )

        if step_index < n_steps
            current_position .= next_position
            next_position .= following_position
            current_velocity .= next_velocity
        end
    end
    return maximum_acceleration
end


"""Build a uniformly sampled profile from a continuous path-speed envelope."""
function rpo_profile_from_speed_envelope(
    samples::Matrix{Float64},
    s_samples::Vector{Float64},
    speed_envelope::Vector{Float64},
    dt_s::Float64,
    maximum_steps::Int,
)
    n_samples = length(s_samples)
    n_samples == length(speed_envelope) || throw(DimensionMismatch(
        "speed envelope must match path samples",
    ))
    n_samples == 1 && return RPORetimingProfile(
        samples[:, 1:1], [0.0], [0.0], [0.0], dt_s,
    )

    segment_duration = zeros(n_samples - 1)
    cumulative_time = zeros(n_samples)
    @inbounds for index in 1:(n_samples - 1)
        ds = s_samples[index + 1] - s_samples[index]
        average_speed = 0.5 * (
            speed_envelope[index] + speed_envelope[index + 1]
        )
        average_speed > 0.0 || return nothing
        segment_duration[index] = ds / average_speed
        cumulative_time[index + 1] =
            cumulative_time[index] + segment_duration[index]
    end
    nominal_duration = cumulative_time[end]
    isfinite(nominal_duration) && nominal_duration > 0.0 || return nothing
    n_steps = ceil(Int, nominal_duration / dt_s)
    n_steps = max(n_steps, 2)
    n_steps <= maximum_steps || return nothing

    # Rounding duration up to a whole number of control updates only slows the
    # envelope, preserving its acceleration bounds and landing on the endpoint.
    time_dilation = n_steps * dt_s / nominal_duration
    s_ref = Vector{Float64}(undef, n_steps + 1)
    speed_ref = Vector{Float64}(undef, n_steps + 1)
    @inbounds for step in 0:n_steps
        nominal_time = min(step * dt_s / time_dilation, nominal_duration)
        segment = clamp(
            searchsortedlast(cumulative_time, nominal_time), 1, n_samples - 1,
        )
        local_time = nominal_time - cumulative_time[segment]
        duration = segment_duration[segment]
        v0 = speed_envelope[segment]
        v1 = speed_envelope[segment + 1]
        acceleration = (v1 - v0) / duration
        s_ref[step + 1] = s_samples[segment] +
            v0 * local_time + 0.5 * acceleration * local_time^2
        speed_ref[step + 1] =
            (v0 + acceleration * local_time) / time_dilation
    end
    s_ref[1] = s_samples[1]
    s_ref[end] = s_samples[end]
    speed_ref[end] = speed_ref[end - 1]
    return RPORetimingProfile(
        samples, s_samples, s_ref, speed_ref, dt_s,
    )
end


"""Construct and validate an HCW-feasible speed envelope along a sampled path."""
function rpo_hcw_acceleration_constrained_profile(
    samples::Matrix{Float64},
    s_samples::Vector{Float64},
    local_speed_limit::Vector{Float64},
    dt_s::Float64,
    mean_motion_radps::Real,
    acceleration_limit_mps2::Real,
    maximum_steps::Int;
    maximum_passes::Int=8,
    corner_stop_fallback::Bool=true,
)
    limit = Float64(acceleration_limit_mps2)
    limit > 0.0 || throw(ArgumentError("acceleration limit must be positive"))
    # A two-point path needs an interior speed knot between its resting endpoints.
    if length(s_samples) == 2
        samples = hcat(samples[:, 1], 0.5 .* (samples[:, 1] + samples[:, 2]), samples[:, 2])
        s_samples = [s_samples[1], 0.5 * sum(s_samples), s_samples[2]]
        local_speed_limit = [local_speed_limit[1], minimum(local_speed_limit), local_speed_limit[2]]
    end
    speed_envelope = copy(local_speed_limit)
    speed_envelope[1] = 0.0
    speed_envelope[end] = 0.0
    tangential_limit = 0.5 * limit

    # Bound speed changes in both directions before materializing time samples.
    @inbounds for index in 2:length(speed_envelope)
        ds = s_samples[index] - s_samples[index - 1]
        reachable = sqrt(max(
            0.0, speed_envelope[index - 1]^2 + 2.0 * tangential_limit * ds,
        ))
        speed_envelope[index] = min(speed_envelope[index], reachable)
    end
    @inbounds for index in (length(speed_envelope) - 1):-1:1
        ds = s_samples[index + 1] - s_samples[index]
        reachable = sqrt(max(
            0.0, speed_envelope[index + 1]^2 + 2.0 * tangential_limit * ds,
        ))
        speed_envelope[index] = min(speed_envelope[index], reachable)
    end

    speed_scale = 1.0
    tolerance = max(1.0e-12, 1.0e-8 * limit)
    for _ in 1:maximum_passes
        candidate = rpo_profile_from_speed_envelope(
            samples,
            s_samples,
            speed_scale .* speed_envelope,
            dt_s,
            maximum_steps,
        )
        candidate === nothing && break
        maximum_acceleration = rpo_hcw_feedforward_max_acceleration(
            candidate, mean_motion_radps,
        )
        maximum_acceleration <= limit + tolerance && return candidate
        speed_scale *= min(0.98, 0.98 * sqrt(limit / maximum_acceleration))
    end
    if corner_stop_fallback && length(s_samples) > 2
        # A polyline can be traversed by stopping at its corners. Add a midpoint
        # on each segment so even adjacent stops have room to accelerate.
        n = length(s_samples)
        fallback_samples = Matrix{Float64}(undef, 3, 2 * n - 1)
        fallback_s = Vector{Float64}(undef, 2 * n - 1)
        fallback_speed = Vector{Float64}(undef, 2 * n - 1)
        for index in 1:n
            k = 2 * index - 1
            fallback_samples[:, k] .= samples[:, index]
            fallback_s[k] = s_samples[index]
            fallback_speed[k] = local_speed_limit[index]
            if index == 1 || index == n
                fallback_speed[k] = 0.0
            else
                incoming = (samples[:, index] - samples[:, index - 1]) /
                    (s_samples[index] - s_samples[index - 1])
                outgoing = (samples[:, index + 1] - samples[:, index]) /
                    (s_samples[index + 1] - s_samples[index])
                norm(outgoing - incoming) > 1.0e-8 && (fallback_speed[k] = 0.0)
            end
            if index < n
                fallback_samples[:, k + 1] .= 0.5 .* (samples[:, index] + samples[:, index + 1])
                fallback_s[k + 1] = 0.5 * (s_samples[index] + s_samples[index + 1])
                fallback_speed[k + 1] = min(local_speed_limit[index], local_speed_limit[index + 1])
            end
        end
        return rpo_hcw_acceleration_constrained_profile(
            fallback_samples, fallback_s, fallback_speed, dt_s,
            mean_motion_radps, acceleration_limit_mps2, maximum_steps;
            maximum_passes=maximum_passes, corner_stop_fallback=false,
        )
    end
    return nothing
end


"""Build a retiming profile from already-sampled path points."""
function rpo_retiming_profile_from_samples(
    sampled_points,
    geometry,
    cfg::RPOPSOConfig;
    fallback_speed_mps::Real=1.0e-3,
    duplicate_tol_m::Real=1.0e-10,
    geometry_distance=nothing,
    samples_are_clean::Bool=false,
    retime_dt_s::Real=cfg.retime_dt_s,
    retime_mean_motion_radps=nothing,
    retime_command_limit_mps2=nothing,
)
    raw_samples = Matrix{Float64}(sampled_points)
    samples = samples_are_clean ? raw_samples :
        rpo_remove_near_duplicate_samples(raw_samples; tol=duplicate_tol_m)
    supplied_geometry_distance = geometry_distance === nothing ? nothing :
        Vector{Float64}(geometry_distance)
    if supplied_geometry_distance !== nothing && size(samples, 2) != length(supplied_geometry_distance)
        throw(DimensionMismatch("geometry_distance must match the cleaned sample count."))
    end

    s_samples = rpo_arc_length_params(samples)
    total = s_samples[end]
    dt = Float64(retime_dt_s)
    dt > 0.0 || throw(ArgumentError("retime_dt_s must be positive."))

    if total <= eps(Float64)
        @warn "RPO retiming received a zero-length path. Returning a single-point reference."
        return RPORetimingProfile(samples[:, 1:1], [0.0], [0.0], [0.0], dt)
    end

    κ = rpo_curvature_from_samples(samples, s_samples)

    n = length(s_samples)
    v_max = zeros(n)
    local_geometry_distance = supplied_geometry_distance === nothing ? zeros(n) :
        supplied_geometry_distance

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
        if supplied_geometry_distance === nothing
            local_geometry_distance[j] = rpo_clearance_to_station(samples[:, j], geometry).distance
        end

        d_avail = max(0.0, local_geometry_distance[j])

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
            geometry_distance_m = local_geometry_distance[k],
            curvature_1pm = κ[k],
            fallback_speed_mps = fallback_speed,
        )
    end

    if maximum(v_max) <= 0.0
        @warn "RPO retiming found no feasible positive speeds. Using constant fallback speed for entire path." fallback_speed_mps=fallback_speed
        fill!(v_max, fallback_speed)
    end

    if retime_mean_motion_radps !== nothing &&
       retime_command_limit_mps2 !== nothing
        command_limit = min(
            Float64(retime_command_limit_mps2), Float64(cfg.retime_a_max_mps2),
        )
        return rpo_hcw_acceleration_constrained_profile(
            samples,
            s_samples,
            v_max,
            dt,
            Float64(retime_mean_motion_radps),
            command_limit,
            cfg.retime_max_steps,
        )
    end

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
                geometry_distance_m = local_geometry_distance[idx],
                curvature_1pm = κ[idx],
            )

            v = fallback_speed
        end

        push!(s_hist, s)
        push!(v_hist, v)
        s >= total - 1.0e-9 && break

        steps += 1

        if steps > cfg.retime_max_steps
            if s < total
                push!(s_hist, total)
                push!(v_hist, fallback_speed)
            end

            break
        end

        ds = v * dt

        if ds <= eps(Float64) || !isfinite(ds)
            @warn "RPO retiming computed invalid arc-length step. Using fallback step." (
                s_m = s,
                total_s_m = total,
                v_mps = v,
                dt_s = dt,
                fallback_speed_mps = fallback_speed,
            )

            ds = fallback_speed * dt
        end

        s = min(total, s + ds)
    end

    # A moving rest-to-rest reference needs separate departure and braking intervals.
    if length(s_hist) == 2
        insert!(s_hist, 2, 0.5 * total)
        insert!(v_hist, 2, 0.5 * total / dt)
    end
    return RPORetimingProfile(samples, s_samples, s_hist, v_hist, dt)
end


"""Materialize position samples from a compact RPO retiming profile."""
function rpo_retime_path_from_profile(profile::RPORetimingProfile)
    mat = zeros(3, length(profile.s_ref))
    @inbounds for j in eachindex(profile.s_ref)
        rpo_interpolate_along_path!(
            view(mat, :, j), profile.samples, profile.s_samples, profile.s_ref[j],
        )
    end
    return mat, copy(profile.s_ref), copy(profile.speed_ref)
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
    profile = rpo_retiming_profile_from_samples(
        raw_samples,
        geometry,
        cfg;
        fallback_speed_mps=fallback_speed_mps,
        duplicate_tol_m=duplicate_tol_m,
    )
    return rpo_retime_path_from_profile(profile)
end
