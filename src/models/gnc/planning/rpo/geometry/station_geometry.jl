"""KD-tree node for nearest-neighbor queries on station point clouds."""
mutable struct RPOStationKDNode
    idx::Int
    axis::Int
    left::Union{Nothing, RPOStationKDNode}
    right::Union{Nothing, RPOStationKDNode}
end

"""
    RPOStationGeometry

Target-station geometry for RPO clearance queries: body-frame point samples in
metres (`points_body`, 3 x N), a KD-tree over them for nearest-point queries, the
keep-out radius around every point (`keepout_radius_m`) and a `name`. Build one
with `RPOStationGeometry(points_body; keepout_radius_m, name)` and pair it with an
`RPOCubeSatGeometry` chaser in an `RPOReferenceGeometry`.
"""
struct RPOStationGeometry
    points_body::Matrix{Float64}
    kd_root::Union{Nothing, RPOStationKDNode}
    keepout_radius_m::Float64
    name::String
end

"""Build a KD-tree over station point-cloud indices."""
function _rpo_build_station_kdtree(points::Matrix{Float64}, indices::Vector{Int}, depth::Int=0)
    isempty(indices) && return nothing
    axis = mod(depth, 3) + 1
    sort!(indices; by=idx -> points[axis, idx])
    mid = cld(length(indices), 2)
    left_indices = mid > 1 ? indices[1:(mid - 1)] : Int[]
    right_indices = mid < length(indices) ? indices[(mid + 1):end] : Int[]
    return RPOStationKDNode(
        indices[mid],
        axis,
        _rpo_build_station_kdtree(points, left_indices, depth + 1),
        _rpo_build_station_kdtree(points, right_indices, depth + 1),
    )
end

"""
    RPOStationGeometry(points_body; keepout_radius_m=0.0, name="station")

Build station geometry and its KD-tree from body-frame point samples, a 3 x N
matrix in metres with at least one point. `keepout_radius_m` must be nonnegative.
"""
function RPOStationGeometry(points_body; keepout_radius_m::Real=0.0, name::AbstractString="station")
    points = Matrix{Float64}(points_body)
    size(points, 1) == 3 || throw(ArgumentError("RPOStationGeometry point cloud must be a 3 x N matrix."))
    size(points, 2) > 0 || throw(ArgumentError("RPOStationGeometry point cloud must contain at least one point."))
    keepout = Float64(keepout_radius_m)
    keepout >= 0.0 || throw(ArgumentError("RPOStationGeometry keepout_radius_m must be nonnegative."))
    kd_root = _rpo_build_station_kdtree(points, collect(1:size(points, 2)))
    return RPOStationGeometry(points, kd_root, keepout, String(name))
end
