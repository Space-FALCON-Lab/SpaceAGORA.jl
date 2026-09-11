module RPOStationAssets

using LinearAlgebra
using Random

export station_geometry_path, station_cad_path, load_rpo_station_pointcloud, load_rpo_station_cad_triangles, load_rpo_station_cad_pointcloud

const _ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const _STATION_GEOMETRY_ROOT = joinpath(_ROOT, "data", "rpo", "station_geometry")

function station_geometry_path(kind::Symbol=:demo)
    if kind == :demo
        return joinpath(_STATION_GEOMETRY_ROOT, "demo", "station_pointcloud.csv")
    elseif kind in (:gateway, :gateway_core)
        return station_cad_path(kind)
    end
    throw(ArgumentError("Unknown RPO station geometry kind $(kind). Use :demo, :gateway, or add an artifact-backed loader."))
end

function station_cad_path(kind::Symbol=:gateway)
    if kind in (:gateway, :gateway_core)
        return joinpath(_STATION_GEOMETRY_ROOT, "gateway", "Gateway_Core.stl")
    end
    throw(ArgumentError("Unknown RPO station CAD kind $(kind). Use :gateway or add an artifact-backed loader."))
end

function load_rpo_station_pointcloud(kind::Symbol=:demo)
    kind == :demo || throw(ArgumentError("RPO station pointcloud kind $(kind) is not available. Use station_cad_path for CAD assets."))
    path = station_geometry_path(kind)
    isfile(path) || throw(ArgumentError("RPO station pointcloud not found at $(path)."))
    rows = readlines(path)
    points = Float64[]
    for row in rows
        isempty(strip(row)) && continue
        startswith(strip(row), "#") && continue
        vals = parse.(Float64, split(row, ","))
        length(vals) == 3 || throw(ArgumentError("Expected x,y,z rows in $(path)."))
        append!(points, vals)
    end
    n = length(points) ÷ 3
    return reshape(points, 3, n)
end

function _stl_is_binary(path::AbstractString)
    stat(path).size < 84 && return false
    header = Vector{UInt8}(undef, 84)
    open(path, "r") do io
        read!(io, header)
    end
    ntri = reinterpret(UInt32, header[81:84])[1]
    return stat(path).size == 84 + 50 * Int(ntri)
end

function _load_stl_triangles(path::AbstractString)
    if _stl_is_binary(path)
        return open(path, "r") do io
            header = Vector{UInt8}(undef, 80)
            read!(io, header)
            count_bytes = Vector{UInt8}(undef, 4)
            read!(io, count_bytes)
            ntri = Int(reinterpret(UInt32, count_bytes)[1])
            triangles = zeros(Float64, 3, 3 * ntri)
            for tri_idx in 1:ntri
                normal = Vector{Float32}(undef, 3)
                read!(io, normal)
                base = 3 * (tri_idx - 1)
                v1 = Vector{Float32}(undef, 3)
                v2 = Vector{Float32}(undef, 3)
                v3 = Vector{Float32}(undef, 3)
                read!(io, v1)
                read!(io, v2)
                read!(io, v3)
                triangles[:, base + 1] .= v1
                triangles[:, base + 2] .= v2
                triangles[:, base + 3] .= v3
                read(io, UInt16)
            end
            triangles
        end
    end

    verts = Vector{Vector{Float64}}()
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            if startswith(s, "vertex")
                parts = split(s)
                push!(verts, [parse(Float64, parts[2]), parse(Float64, parts[3]), parse(Float64, parts[4])])
            end
        end
    end
    triangles = zeros(Float64, 3, length(verts))
    for (idx, vertex) in enumerate(verts)
        triangles[:, idx] .= vertex
    end
    return triangles
end

function _center_triangles!(triangles::Matrix{Float64})
    bounds_min = vec(minimum(triangles; dims=2))
    bounds_max = vec(maximum(triangles; dims=2))
    triangles .-= (bounds_min .+ bounds_max) ./ 2.0
    return triangles
end

function _sample_triangle_point(v1, v2, v3, rng)
    u = rand(rng)
    v = rand(rng)
    if u + v > 1.0
        u = 1.0 - u
        v = 1.0 - v
    end
    return v1 .+ u .* (v2 .- v1) .+ v .* (v3 .- v1)
end

function _stl_pointcloud(path::AbstractString; n_points::Integer=10000, rng=Random.default_rng(), center::Bool=true, scale::Real=1.0)
    n_points > 0 || throw(ArgumentError("n_points must be positive."))
    triangles = _load_stl_triangles(path)
    scale != 1.0 && (triangles .*= Float64(scale))
    center && _center_triangles!(triangles)

    ntri = size(triangles, 2) ÷ 3
    ntri > 0 || throw(ArgumentError("STL at $(path) contains no triangles."))
    areas = zeros(Float64, ntri)
    for tri_idx in 1:ntri
        base = 3 * (tri_idx - 1)
        v1 = triangles[:, base + 1]
        v2 = triangles[:, base + 2]
        v3 = triangles[:, base + 3]
        areas[tri_idx] = 0.5 * norm(cross(v2 - v1, v3 - v1))
    end
    total_area = sum(areas)
    total_area > 0.0 || throw(ArgumentError("STL at $(path) has zero total triangle area."))
    cdf = cumsum(areas) ./ total_area

    points = zeros(Float64, 3, n_points)
    for point_idx in 1:n_points
        tri_idx = searchsortedfirst(cdf, rand(rng))
        tri_idx = clamp(tri_idx, 1, ntri)
        base = 3 * (tri_idx - 1)
        points[:, point_idx] .= _sample_triangle_point(
            triangles[:, base + 1],
            triangles[:, base + 2],
            triangles[:, base + 3],
            rng,
        )
    end
    return points
end

function load_rpo_station_cad_triangles(kind::Symbol=:gateway; center::Bool=true, scale::Real=1.0)
    path = station_cad_path(kind)
    isfile(path) || throw(ArgumentError("RPO station CAD not found at $(path)."))
    triangles = _load_stl_triangles(path)
    scale != 1.0 && (triangles .*= Float64(scale))
    center && _center_triangles!(triangles)
    return triangles
end

function load_rpo_station_cad_pointcloud(kind::Symbol=:gateway; n_points::Integer=10000, rng=MersenneTwister(42), center::Bool=true, scale::Real=1.0)
    path = station_cad_path(kind)
    isfile(path) || throw(ArgumentError("RPO station CAD not found at $(path)."))
    return _stl_pointcloud(path; n_points=n_points, rng=rng, center=center, scale=scale)
end

end # module RPOStationAssets
