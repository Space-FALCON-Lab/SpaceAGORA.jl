# RPO curve geometry used by sampling, retiming and reference evaluation.
# Definitions retain their GuidanceHooks identity.
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
