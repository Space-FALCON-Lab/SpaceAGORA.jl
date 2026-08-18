"""Return signed clearance from a body-frame point to the station keepout surface."""
function rpo_clearance_to_station(p_body, geometry::RPOReferenceGeometry)
    nearest, distance, idx = nearest_station_point(p_body, geometry.station)
    clearance = distance - geometry.station.keepout_radius_m - maximum(geometry.chaser.half_extents_body)
    return (clearance=clearance, distance=distance, nearest_point=nearest, nearest_index=idx)
end

"""Return only the signed clearance distance to the station keepout surface."""
@inline function rpo_clearance_distance_to_station(p_body, geometry::RPOReferenceGeometry)
    distance = sqrt(nearest_station_distance_sq(p_body, geometry.station))
    return distance - geometry.station.keepout_radius_m - maximum(geometry.chaser.half_extents_body)
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
