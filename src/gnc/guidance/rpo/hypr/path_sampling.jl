function rpo_path_length(points)
    pts = Matrix{Float64}(points)
    size(pts, 2) <= 1 && return 0.0
    len = 0.0
    @inbounds for j in 1:(size(pts, 2) - 1)
        dx = pts[1, j + 1] - pts[1, j]
        dy = pts[2, j + 1] - pts[2, j]
        dz = pts[3, j + 1] - pts[3, j]
        len += sqrt(dx * dx + dy * dy + dz * dz)
    end
    return len
end

function rpo_bezier_point(points, t::Real)
    pts = Matrix{Float64}(points)
    work = similar(pts)
    out = zeros(3)
    rpo_bezier_point!(out, work, pts, Float64(t))
    return out
end

function rpo_bezier_point!(out, work, points, t::Float64)
    n_pts = size(points, 2)
    @inbounds for j in 1:n_pts
        work[1, j] = points[1, j]
        work[2, j] = points[2, j]
        work[3, j] = points[3, j]
    end
    omt = 1.0 - t
    @inbounds for level in 1:(n_pts - 1)
        for j in 1:(n_pts - level)
            work[1, j] = omt * work[1, j] + t * work[1, j + 1]
            work[2, j] = omt * work[2, j] + t * work[2, j + 1]
            work[3, j] = omt * work[3, j] + t * work[3, j + 1]
        end
    end
    out[1] = work[1, 1]
    out[2] = work[2, 1]
    out[3] = work[3, 1]
    return out
end

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

function rpo_sample_path_polyline(points, ds::Real)
    pts = Matrix{Float64}(points)
    total = rpo_path_length(pts)
    n = max(2, Int(ceil(total / Float64(ds))) + 1)
    return rpo_resample_polyline_points(pts, n)
end

function rpo_sample_path(points, ds::Real; curve_type::Symbol=:bezier)
    curve_type == :bezier && return rpo_sample_path_bezier(points, ds)
    curve_type == :polyline && return rpo_sample_path_polyline(points, ds)
    throw(ArgumentError("Unsupported RPO path curve_type=$(curve_type). Use :bezier or :polyline."))
end
