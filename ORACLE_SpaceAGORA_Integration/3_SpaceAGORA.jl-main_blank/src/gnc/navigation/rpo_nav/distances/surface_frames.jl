"""Estimate a station surface normal from the nearest point-cloud sample."""
function rpo_surface_normal_from_pointcloud(p_body, geometry::RPOReferenceGeometry)
    nearest, distance, _ = nearest_station_point(p_body, geometry.station)
    v = SVector{3, Float64}(p_body) - nearest
    if norm(v) <= eps(Float64)
        return SVector{3, Float64}(1.0, 0.0, 0.0)
    end
    return normalize(v)
end

