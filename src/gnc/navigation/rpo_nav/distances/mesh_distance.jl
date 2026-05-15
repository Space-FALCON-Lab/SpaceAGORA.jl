@inline function _rpo_station_point(points::Matrix{Float64}, idx::Int)
    return SVector{3, Float64}(points[1, idx], points[2, idx], points[3, idx])
end

function _nearest_station_point_kdtree(
    node::Nothing,
    p::SVector{3, Float64},
    points::Matrix{Float64},
    best_idx::Int,
    best_d2::Float64,
)
    return best_idx, best_d2
end

function _nearest_station_point_kdtree(
    node::RPOStationKDNode,
    p::SVector{3, Float64},
    points::Matrix{Float64},
    best_idx::Int,
    best_d2::Float64,
)
    q = _rpo_station_point(points, node.idx)
    d2 = sum(abs2, p - q)
    if d2 < best_d2
        best_idx = node.idx
        best_d2 = d2
    end

    axis = node.axis
    diff = p[axis] - points[axis, node.idx]
    near_branch = diff <= 0.0 ? node.left : node.right
    far_branch = diff <= 0.0 ? node.right : node.left
    best_idx, best_d2 = _nearest_station_point_kdtree(near_branch, p, points, best_idx, best_d2)
    if diff * diff < best_d2
        best_idx, best_d2 = _nearest_station_point_kdtree(far_branch, p, points, best_idx, best_d2)
    end
    return best_idx, best_d2
end

function _nearest_station_distance_sq_kdtree(
    node::Nothing,
    p::SVector{3, Float64},
    points::Matrix{Float64},
    best_d2::Float64,
)
    return best_d2
end

function _nearest_station_distance_sq_kdtree(
    node::RPOStationKDNode,
    p::SVector{3, Float64},
    points::Matrix{Float64},
    best_d2::Float64,
)
    dx = p[1] - points[1, node.idx]
    dy = p[2] - points[2, node.idx]
    dz = p[3] - points[3, node.idx]
    d2 = dx * dx + dy * dy + dz * dz
    d2 < best_d2 && (best_d2 = d2)

    axis = node.axis
    diff = p[axis] - points[axis, node.idx]
    near_branch = diff <= 0.0 ? node.left : node.right
    far_branch = diff <= 0.0 ? node.right : node.left
    best_d2 = _nearest_station_distance_sq_kdtree(near_branch, p, points, best_d2)
    if diff * diff < best_d2
        best_d2 = _nearest_station_distance_sq_kdtree(far_branch, p, points, best_d2)
    end
    return best_d2
end

@inline function nearest_station_distance_sq(p_body, station::RPOStationGeometry)
    p = SVector{3, Float64}(p_body)
    station.kd_root === nothing && throw(ArgumentError("RPO station KD-tree is empty."))
    return _nearest_station_distance_sq_kdtree(station.kd_root, p, station.points_body, Inf)
end

@inline function nearest_station_point(p_body, station::RPOStationGeometry)
    p = SVector{3, Float64}(p_body)
    station.kd_root === nothing && throw(ArgumentError("RPO station KD-tree is empty."))
    best_idx, best_d2 = _nearest_station_point_kdtree(station.kd_root, p, station.points_body, 1, Inf)
    q = _rpo_station_point(station.points_body, best_idx)
    return q, sqrt(best_d2), best_idx
end
