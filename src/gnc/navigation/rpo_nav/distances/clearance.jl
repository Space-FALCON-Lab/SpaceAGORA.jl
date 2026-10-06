"""Return signed clearance from a body-frame point to the station keepout surface."""
function rpo_clearance_to_station(p_body, geometry::RPOReferenceGeometry)
    nearest, distance, idx = nearest_station_point(p_body, geometry.station)
    clearance = distance - geometry.station.keepout_radius_m - norm(geometry.chaser.half_extents_body)
    return (clearance=clearance, distance=distance, nearest_point=nearest, nearest_index=idx)
end

"""Return only the signed clearance distance to the station keepout surface."""
@inline function rpo_clearance_distance_to_station(p_body, geometry::RPOReferenceGeometry)
    distance = sqrt(nearest_station_distance_sq(p_body, geometry.station))
    return distance - geometry.station.keepout_radius_m - norm(geometry.chaser.half_extents_body)
end

"""Compute minimum clearance and violation counts along a body-frame path."""
function rpo_path_clearance_stats(path_body, geometry::RPOReferenceGeometry; safe_distance_m::Real=0.0)
    path = Matrix{Float64}(path_body)
    size(path, 1) == 3 || throw(ArgumentError("RPO path must be a 3 x N matrix."))
    margin = Float64(safe_distance_m)
    min_clearance = Inf
    violation_count = 0
    @inbounds for j in 1:size(path, 2)
        p = SVector{3, Float64}(path[1, j], path[2, j], path[3, j])
        clearance = rpo_clearance_distance_to_station(p, geometry)
        min_clearance = min(min_clearance, clearance)
        violation_count += clearance < margin ? 1 : 0
    end
    return (
        min_clearance=min_clearance,
        violation_count=violation_count,
        violation_fraction=violation_count / max(size(path, 2), 1),
    )
end

"""Nearest station-point distance to an entire segment, pruning by KD split planes."""
function _rpo_segment_distance_sq(node, a, b, points, best)
    node === nothing && return best
    q = _rpo_station_point(points, node.idx)
    ab = b - a
    length_sq = sum(abs2, ab)
    t = length_sq == 0.0 ? 0.0 : clamp(dot(q - a, ab) / length_sq, 0.0, 1.0)
    best = min(best, sum(abs2, q - (a + t * ab)))
    axis = node.axis
    split = q[axis]
    # A whole branch lies on one side of the split plane. Its distance to
    # the segment cannot be smaller than this coordinate separation.
    left_gap = max(min(a[axis], b[axis]) - split, 0.0)
    right_gap = max(split - max(a[axis], b[axis]), 0.0)
    if left_gap^2 <= best
        best = _rpo_segment_distance_sq(node.left, a, b, points, best)
    end
    if right_gap^2 <= best
        best = _rpo_segment_distance_sq(node.right, a, b, points, best)
    end
    return best
end

"""Signed capsule clearance, excluding the caller's additional safety margin.

The capsule encloses the chaser box at every attitude using its half-diagonal.
The station point cloud is inflated by the station keepout radius. This measures
an entire straight segment, including its endpoints; it does not bound CAD
surface features absent from the point cloud.
"""
function rpo_capsule_clearance_to_station(start_body, end_body, geometry::RPOReferenceGeometry)
    a = SVector{3, Float64}(start_body)
    b = SVector{3, Float64}(end_body)
    station = geometry.station
    station.kd_root === nothing && throw(ArgumentError("RPO station KD-tree is empty."))
    initial = nearest_station_distance_sq(a, station)
    distance_sq = _rpo_segment_distance_sq(station.kd_root, a, b, station.points_body, initial)
    return sqrt(distance_sq) - station.keepout_radius_m - norm(geometry.chaser.half_extents_body)
end
