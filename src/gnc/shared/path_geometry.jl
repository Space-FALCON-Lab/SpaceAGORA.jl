# Common path geometry for RPO and robot-arm planning.
# Loaded into HYPRUtils to preserve existing qualified access.

"""Return the Euclidean arc length of a waypoint matrix."""
function hypr_path_length(points)
    pts = Matrix{Float64}(points)
    size(pts, 2) < 2 && return 0.0
    total = 0.0
    @inbounds for j in 1:(size(pts, 2) - 1)
        total += norm(pts[:, j + 1] - pts[:, j])
    end
    return total
end

"""Evaluate a Bezier curve at a normalized parameter without mutating caller-owned work buffers."""
function hypr_bezier_point(points, t::Real)
    pts = Matrix{Float64}(points)
    out = zeros(size(pts, 1))
    work = similar(pts)
    return hypr_bezier_point!(out, work, pts, Float64(t))
end

"""Evaluate a Bezier curve in-place using caller-provided output and work buffers."""
function hypr_bezier_point!(out, work, points, t::Float64)
    n = size(points, 2)
    work[:, 1:n] .= points
    @inbounds for r in 1:(n - 1)
        for j in 1:(n - r)
            work[:, j] .= (1 - t) .* work[:, j] .+ t .* work[:, j + 1]
        end
    end
    out .= work[:, 1]
    return out
end

"""Sample either a Bezier or polyline path at a fixed number of points."""
function hypr_sample_count_path(points, n_samples::Int; curve_type::Symbol=:bezier)
    n_samples >= 2 || throw(ArgumentError("n_samples must be at least 2."))
    pts = Matrix{Float64}(points)
    n_dims = size(pts, 1)
    n_points = size(pts, 2)
    samples = zeros(n_dims, n_samples)
    if curve_type == :bezier
        degree = n_points - 1
        @inbounds for k in 1:n_samples
            s = (k - 1) / max(n_samples - 1, 1)
            for i in 0:degree
                coeff = binomial(degree, i) * (1.0 - s)^(degree - i) * s^i
                samples[:, k] .+= coeff .* pts[:, i + 1]
            end
        end
    elseif curve_type == :polyline
        @inbounds for k in 1:n_samples
            u = (n_points - 1) * (k - 1) / max(n_samples - 1, 1)
            seg = clamp(floor(Int, u) + 1, 1, n_points - 1)
            alpha = clamp(u - (seg - 1), 0.0, 1.0)
            samples[:, k] .= (1.0 - alpha) .* pts[:, seg] .+ alpha .* pts[:, seg + 1]
        end
    else
        throw(ArgumentError("curve_type must be :bezier or :polyline."))
    end
    samples[:, 1] .= pts[:, 1]
    samples[:, end] .= pts[:, end]
    return samples
end
