# Triangle-soup readers for CAD/model files (STL, OBJ, glTF/GLB with embedded
# buffers) and the geometry helpers built on them: bounding boxes, area-weighted
# surface point clouds, and the model-axes transform (`scale`, XYZ Euler
# rotation in degrees) shared with the interactive viewer. Owned by the
# structure layer so both the visualization scene (centeing, RPO station
# point clouds) and the aerodynamic panel method (`aerodynamic_mesh_surrogate.jl`)
# read the same geometry from the same file.

# 3D model overrides the viewer can parse in the browser, by file extension.
const MODEL_FORMATS = Dict{String, Tuple{String, String}}(
    ".stl" => ("stl", "model/stl"),
    ".obj" => ("obj", "model/obj"),
    ".glb" => ("glb", "model/gltf-binary"),
    ".gltf" => ("gltf", "model/gltf+json"),
)
"""
    model_format(path) -> (format, mime)

Viewer model format from the file extension: `"stl"`, `"obj"`, `"glb"` or
`"gltf"`. Throws for anything else.
"""
function model_format(path::AbstractString)::Tuple{String, String}
    ext = lowercase(splitext(String(path))[2])
    haskey(MODEL_FORMATS, ext) || throw(ArgumentError("Unsupported 3D model format $(repr(ext)) for $(path); use .stl, .obj, .glb or .gltf."))
    return MODEL_FORMATS[ext]
end

# Extensions three's GLTFLoader can only handle with an extra decoder the
# page does not carry (Draco, meshopt, KTX2/Basis textures).
const GLTF_UNSUPPORTED_REQUIRED = ("KHR_draco_mesh_compression", "EXT_meshopt_compression", "KHR_texture_basisu")

"""
    gltf_required_extensions(path) -> Vector{String}

`extensionsRequired` of a `.glb` or `.gltf` file (empty when none).
"""
function gltf_required_extensions(path::AbstractString)::Vector{String}
    bytes = read(path)
    json_text = if length(bytes) >= 20 && bytes[1:4] == Vector{UInt8}("glTF")
        chunk_len = Int(reinterpret(UInt32, bytes[13:16])[1])
        String(bytes[21:min(20 + chunk_len, length(bytes))])
    else
        String(bytes)
    end
    parsed = try
        JSON.parse(json_text)
    catch
        return String[]
    end
    parsed isa AbstractDict || return String[]
    return String[String(x) for x in get(parsed, "extensionsRequired", Any[])]
end

@inline function _rotation_xyz_deg(rot)::SMatrix{3, 3, Float64}
    rx, ry, rz = deg2rad(Float64(rot[1])), deg2rad(Float64(rot[2])), deg2rad(Float64(rot[3]))
    Rx = @SMatrix [1.0 0.0 0.0; 0.0 cos(rx) -sin(rx); 0.0 sin(rx) cos(rx)]
    Ry = @SMatrix [cos(ry) 0.0 sin(ry); 0.0 1.0 0.0; -sin(ry) 0.0 cos(ry)]
    Rz = @SMatrix [cos(rz) -sin(rz) 0.0; sin(rz) cos(rz) 0.0; 0.0 0.0 1.0]
    # three.js Euler order XYZ: v' = Rx * Ry * Rz * v
    return Rx * Ry * Rz
end

# ---------------------------------------------------------------------------
# STL
# ---------------------------------------------------------------------------

function _stl_binary_triangles(bytes::Vector{UInt8})::Matrix{Float64}
    ntri = Int(reinterpret(UInt32, bytes[81:84])[1])
    out = Matrix{Float64}(undef, 3, 3 * ntri)
    off = 84
    @inbounds for t in 1:ntri
        base = off + 50 * (t - 1) + 12  # skip the normal
        for v in 1:3
            for c in 1:3
                out[c, 3 * (t - 1) + v] = Float64(reinterpret(Float32, bytes[base + 12 * (v - 1) + 4 * (c - 1) + 1:base + 12 * (v - 1) + 4 * c])[1])
            end
        end
    end
    return out
end

function _stl_ascii_triangles(text::String)::Matrix{Float64}
    verts = Float64[]
    for line in eachline(IOBuffer(text))
        s = strip(line)
        startswith(s, "vertex") || continue
        parts = split(s)
        length(parts) >= 4 || continue
        append!(verts, (parse(Float64, parts[2]), parse(Float64, parts[3]), parse(Float64, parts[4])))
    end
    n = length(verts) ÷ 9
    return reshape(verts[1:9 * n], 3, 3 * n)
end

function _load_stl_triangles(path::AbstractString)::Matrix{Float64}
    bytes = read(path)
    if length(bytes) >= 84
        ntri = Int(reinterpret(UInt32, bytes[81:84])[1])
        length(bytes) == 84 + 50 * ntri && return _stl_binary_triangles(bytes)
    end
    return _stl_ascii_triangles(String(bytes))
end

# ---------------------------------------------------------------------------
# OBJ (geometry only; polygons are fan-triangulated)
# ---------------------------------------------------------------------------

function _load_obj_triangles(path::AbstractString)::Matrix{Float64}
    vertices = SVector{3, Float64}[]
    tris = Float64[]
    for line in eachline(path)
        s = strip(line)
        if startswith(s, "v ")
            parts = split(s)
            push!(vertices, SVector{3, Float64}(parse(Float64, parts[2]), parse(Float64, parts[3]), parse(Float64, parts[4])))
        elseif startswith(s, "f ")
            idx = Int[]
            for token in split(s)[2:end]
                i = parse(Int, first(split(token, "/")))
                push!(idx, i < 0 ? length(vertices) + 1 + i : i)
            end
            for k in 2:(length(idx) - 1)
                for v in (vertices[idx[1]], vertices[idx[k]], vertices[idx[k + 1]])
                    append!(tris, v)
                end
            end
        end
    end
    return reshape(tris, 3, length(tris) ÷ 3)
end

# ---------------------------------------------------------------------------
# glTF 2.0 / GLB (embedded buffers only)
# ---------------------------------------------------------------------------

const _GLTF_COMPONENT_TYPES = Dict{Int, DataType}(5120 => Int8, 5121 => UInt8, 5122 => Int16, 5123 => UInt16, 5125 => UInt32, 5126 => Float32)
const _GLTF_TYPE_COUNTS = Dict{String, Int}("SCALAR" => 1, "VEC2" => 2, "VEC3" => 3, "VEC4" => 4, "MAT4" => 16)

function _gltf_document(path::AbstractString)
    bytes = read(path)
    if length(bytes) >= 20 && bytes[1:4] == Vector{UInt8}("glTF")
        json_len = Int(reinterpret(UInt32, bytes[13:16])[1])
        doc = JSON.parse(String(bytes[21:20 + json_len]))
        bin = UInt8[]
        pos = 20 + json_len + 1
        if pos + 8 <= length(bytes)
            bin_len = Int(reinterpret(UInt32, bytes[pos:pos + 3])[1])
            bin = bytes[pos + 8:min(pos + 7 + bin_len, length(bytes))]
        end
        return doc, [bin]
    end
    doc = JSON.parse(String(bytes))
    buffers = Vector{UInt8}[]
    for b in get(doc, "buffers", Any[])
        uri = get(b, "uri", nothing)
        uri === nothing && throw(ArgumentError("$(basename(path)): glTF buffer without a data URI; embed the buffers or use GLB."))
        startswith(uri, "data:") || throw(ArgumentError("$(basename(path)): glTF buffer references an external file ($(uri)); embed the buffers or use GLB."))
        push!(buffers, base64decode(split(uri, ",", limit=2)[2]))
    end
    return doc, buffers
end

function _gltf_accessor(doc, buffers, index::Int)::Matrix{Float64}
    acc = doc["accessors"][index + 1]
    T = _GLTF_COMPONENT_TYPES[Int(acc["componentType"])]
    ncomp = _GLTF_TYPE_COUNTS[String(acc["type"])]
    count = Int(acc["count"])
    out = Matrix{Float64}(undef, ncomp, count)
    haskey(acc, "bufferView") || (fill!(out, 0.0); return out)
    view = doc["bufferViews"][Int(acc["bufferView"]) + 1]
    buf = buffers[Int(get(view, "buffer", 0)) + 1]
    base = Int(get(view, "byteOffset", 0)) + Int(get(acc, "byteOffset", 0))
    elem = sizeof(T) * ncomp
    stride = Int(get(view, "byteStride", elem))
    @inbounds for i in 1:count
        start = base + (i - 1) * stride
        for c in 1:ncomp
            s = start + (c - 1) * sizeof(T) + 1
            out[c, i] = Float64(reinterpret(T, buf[s:s + sizeof(T) - 1])[1])
        end
    end
    return out
end

function _gltf_node_matrix(node)::SMatrix{4, 4, Float64}
    if haskey(node, "matrix")
        m = Float64.(node["matrix"])
        return SMatrix{4, 4, Float64}(m...)  # column-major, as glTF stores it
    end
    t = Float64.(get(node, "translation", [0.0, 0.0, 0.0]))
    q = Float64.(get(node, "rotation", [0.0, 0.0, 0.0, 1.0]))
    s = Float64.(get(node, "scale", [1.0, 1.0, 1.0]))
    R = rot(SVector{4, Float64}(q))'  # active rotation from the scalar-last quaternion
    M = SMatrix{4, 4, Float64}(
        R[1, 1] * s[1], R[2, 1] * s[1], R[3, 1] * s[1], 0.0,
        R[1, 2] * s[2], R[2, 2] * s[2], R[3, 2] * s[2], 0.0,
        R[1, 3] * s[3], R[2, 3] * s[3], R[3, 3] * s[3], 0.0,
        t[1], t[2], t[3], 1.0,
    )
    return M
end

function _gltf_collect_triangles!(tris::Vector{Float64}, doc, buffers, node_index::Int, parent::SMatrix{4, 4, Float64})
    node = doc["nodes"][node_index + 1]
    M = parent * _gltf_node_matrix(node)
    if haskey(node, "mesh")
        mesh = doc["meshes"][Int(node["mesh"]) + 1]
        for prim in mesh["primitives"]
            Int(get(prim, "mode", 4)) == 4 || continue  # triangles only
            haskey(prim["attributes"], "POSITION") || continue
            P = _gltf_accessor(doc, buffers, Int(prim["attributes"]["POSITION"]))
            idx = haskey(prim, "indices") ? Int.(vec(_gltf_accessor(doc, buffers, Int(prim["indices"])))) .+ 1 : collect(1:size(P, 2))
            for k in 1:3:(length(idx) - 2)
                for v in (idx[k], idx[k + 1], idx[k + 2])
                    p = M * SVector{4, Float64}(P[1, v], P[2, v], P[3, v], 1.0)
                    append!(tris, (p[1], p[2], p[3]))
                end
            end
        end
    end
    for child in get(node, "children", Any[])
        _gltf_collect_triangles!(tris, doc, buffers, Int(child), M)
    end
    return tris
end

function _load_gltf_triangles(path::AbstractString)::Matrix{Float64}
    doc, buffers = _gltf_document(path)
    required = String[String(x) for x in get(doc, "extensionsRequired", Any[])]
    bad = [e for e in required if e in GLTF_UNSUPPORTED_REQUIRED]
    isempty(bad) || throw(ArgumentError("$(basename(path)) requires $(join(bad, ", ")); decompress it first (see `gltf_required_extensions`)."))
    tris = Float64[]
    scene_index = Int(get(doc, "scene", 0))
    roots = haskey(doc, "scenes") && !isempty(doc["scenes"]) ? Int.(get(doc["scenes"][scene_index + 1], "nodes", Any[])) : collect(0:(length(get(doc, "nodes", Any[])) - 1))
    identity = SMatrix{4, 4, Float64}(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    for r in roots
        _gltf_collect_triangles!(tris, doc, buffers, r, identity)
    end
    return reshape(tris, 3, length(tris) ÷ 3)
end

# Rodrigues rotation matrix about a unit axis.
@inline function _axis_angle_matrix(axis, angle_deg::Real)::SMatrix{3, 3, Float64}
    a = SVector{3, Float64}(axis)
    n = norm(a)
    n > 0 || throw(ArgumentError("articulation axis must be non-zero."))
    a = a / n
    θ = deg2rad(Float64(angle_deg))
    c, s = cos(θ), sin(θ)
    K = @SMatrix [0.0 -a[3] a[2]; a[3] 0.0 -a[1]; -a[2] a[1] 0.0]
    return SMatrix{3, 3, Float64}(I) + s * K + (1 - c) * (K * K)
end

@inline _articulation_get(spec, key, default) = spec isa AbstractDict ? get(spec, key, get(spec, String(key), default)) : (hasproperty(spec, key) ? getproperty(spec, key) : default)

# The region test of an articulation: an axis-aligned box in model units, any
# bound left out being unbounded.
@inline function _articulation_bounds(spec)
    region = _articulation_get(spec, :region, nothing)
    region === nothing && throw(ArgumentError("an articulation needs a region (x_min/x_max/y_min/y_max/z_min/z_max in model units)."))
    g = (k) -> Float64(_articulation_get(region, k, k in (:x_min, :y_min, :z_min) ? -Inf : Inf))
    return SVector{3, Float64}(g(:x_min), g(:y_min), g(:z_min)), SVector{3, Float64}(g(:x_max), g(:y_max), g(:z_max))
end

"""
    articulate_triangles(tris, articulations) -> Matrix{Float64}

Rotate parts of a triangle soup (3 x 3N, model units) about an axis: each
articulation is a NamedTuple or Dict with `region` (an axis-aligned box in
model units, `x_min`/`x_max`/`y_min`/`y_max`/`z_min`/`z_max`, missing bounds
unbounded), `axis` (rotation axis in model axes), `angle_deg`, and `pivot`
(a point on the axis in model units, or `:centroid` (default) for the center
of the selected vertices' bounding box). Vertices inside the region move;
used to pose parts of a CAD model the file holds in another position, such
as solar wings turned broadside for aerobraking. The viewer applies the same
articulations to the drawn model (`model_articulations`), so the picture and
the aerodynamic mesh agree.
"""
function articulate_triangles(tris::AbstractMatrix{<:Real}, articulations)::Matrix{Float64}
    out = Matrix{Float64}(tris)
    for (k, spec) in enumerate(articulations)
        lo, hi = _articulation_bounds(spec)
        R = _axis_angle_matrix(_articulation_get(spec, :axis, (1.0, 0.0, 0.0)), _articulation_get(spec, :angle_deg, 0.0))
        selected = Int[]
        @inbounds for i in 1:size(out, 2)
            p = SVector{3, Float64}(out[1, i], out[2, i], out[3, i])
            all(lo .<= p .<= hi) && push!(selected, i)
        end
        isempty(selected) && throw(ArgumentError("articulation $(k): no vertices inside its region."))
        pivot_spec = _articulation_get(spec, :pivot, :centroid)
        pivot = if pivot_spec === :centroid || pivot_spec == "centroid"
            mn = SVector{3, Float64}(minimum(out[1, selected]), minimum(out[2, selected]), minimum(out[3, selected]))
            mx = SVector{3, Float64}(maximum(out[1, selected]), maximum(out[2, selected]), maximum(out[3, selected]))
            0.5 * (mn + mx)
        else
            SVector{3, Float64}(pivot_spec)
        end
        @inbounds for i in selected
            p = SVector{3, Float64}(out[1, i], out[2, i], out[3, i])
            q = R * (p - pivot) + pivot
            out[1, i] = q[1]; out[2, i] = q[2]; out[3, i] = q[3]
        end
    end
    return out
end

"""
    articulation_payload(articulations) -> Vector{Dict{String, Any}}

The articulations in the form the viewer's model loader reads: explicit
bounds (infinite bounds as `nothing`), unit axis, angle and pivot (the
centroid resolved against `tris` when given as `:centroid`; pass `tris` to
resolve it).
"""
function articulation_payload(articulations, tris::Union{Nothing, AbstractMatrix{<:Real}}=nothing)::Vector{Dict{String, Any}}
    out = Dict{String, Any}[]
    for spec in articulations
        lo, hi = _articulation_bounds(spec)
        axis = normalize(SVector{3, Float64}(_articulation_get(spec, :axis, (1.0, 0.0, 0.0))))
        pivot_spec = _articulation_get(spec, :pivot, :centroid)
        pivot = if pivot_spec === :centroid || pivot_spec == "centroid"
            tris === nothing && throw(ArgumentError("articulation_payload needs the triangles to resolve a :centroid pivot."))
            sel = [i for i in 1:size(tris, 2) if all(lo .<= SVector{3, Float64}(tris[1, i], tris[2, i], tris[3, i]) .<= hi)]
            isempty(sel) && throw(ArgumentError("articulation: no vertices inside its region."))
            0.5 * (SVector{3, Float64}(minimum(tris[1, sel]), minimum(tris[2, sel]), minimum(tris[3, sel])) + SVector{3, Float64}(maximum(tris[1, sel]), maximum(tris[2, sel]), maximum(tris[3, sel])))
        else
            SVector{3, Float64}(pivot_spec)
        end
        push!(out, Dict{String, Any}(
            "region" => Dict{String, Any}("min" => [isfinite(v) ? v : nothing for v in lo], "max" => [isfinite(v) ? v : nothing for v in hi]),
            "axis" => collect(axis), "angle_deg" => Float64(_articulation_get(spec, :angle_deg, 0.0)), "pivot" => collect(pivot),
        ))
    end
    return out
end

"""
    load_model_triangles(path; scale=1.0, rotation_deg=(0, 0, 0), articulations=()) -> Matrix{Float64}

Triangle vertices (3 x 3N, meters) of an STL, OBJ, glTF or GLB file in model
axes, with `articulations` (see [`articulate_triangles`](@ref)) applied in
model units first, then scaled by `scale` (meters per model unit) and
rotated by XYZ Euler angles in degrees, the same transform the viewer
applies through `model_scale`, `model_rotation_deg` and `model_articulations`.
"""
function load_model_triangles(path::AbstractString; scale::Real=1.0, rotation_deg=(0.0, 0.0, 0.0), articulations=())::Matrix{Float64}
    format, _ = model_format(path)
    tris = format == "stl" ? _load_stl_triangles(path) : format == "obj" ? _load_obj_triangles(path) : _load_gltf_triangles(path)
    isempty(articulations) || (tris = articulate_triangles(tris, articulations))
    R = _rotation_xyz_deg(rotation_deg)
    out = Matrix{Float64}(undef, 3, size(tris, 2))
    @inbounds for i in 1:size(tris, 2)
        p = R * (Float64(scale) * SVector{3, Float64}(tris[1, i], tris[2, i], tris[3, i]))
        out[1, i] = p[1]; out[2, i] = p[2]; out[3, i] = p[3]
    end
    return out
end

"""
    model_bounding_box(path) -> (min::SVector{3}, max::SVector{3})

Axis-aligned bounds of a model in its own units and axes.
"""
function model_bounding_box(path::AbstractString)
    tris = load_model_triangles(path)
    size(tris, 2) > 0 || throw(ArgumentError("$(basename(path)) has no triangles."))
    lo = SVector{3, Float64}(minimum(tris[1, :]), minimum(tris[2, :]), minimum(tris[3, :]))
    hi = SVector{3, Float64}(maximum(tris[1, :]), maximum(tris[2, :]), maximum(tris[3, :]))
    return lo, hi
end

"""
    model_bounding_box_center(path) -> SVector{3, Float64}
"""
model_bounding_box_center(path::AbstractString) = (b = model_bounding_box(path); 0.5 * (b[1] + b[2]))

"""
    sample_model_pointcloud(path; n_points=10000, rng=Random.default_rng(), scale=1.0, rotation_deg=(0,0,0), center=true) -> Matrix{Float64}

Area-weighted random surface samples (3 x n_points, meters) of a model after
the viewer's scale and rotation, optionally translated so its bounding-box
center is at the origin, which is where the viewer puts a centerd model.
"""
function sample_model_pointcloud(path::AbstractString; n_points::Integer=10000, rng=Random.default_rng(), scale::Real=1.0, rotation_deg=(0.0, 0.0, 0.0), center::Bool=true)::Matrix{Float64}
    tris = load_model_triangles(path; scale=scale, rotation_deg=rotation_deg)
    ntri = size(tris, 2) ÷ 3
    ntri > 0 || throw(ArgumentError("$(basename(path)) has no triangles."))
    if center
        c = SVector{3, Float64}(0.5 * (minimum(tris[1, :]) + maximum(tris[1, :])), 0.5 * (minimum(tris[2, :]) + maximum(tris[2, :])), 0.5 * (minimum(tris[3, :]) + maximum(tris[3, :])))
        tris = tris .- c
    end
    areas = Vector{Float64}(undef, ntri)
    @inbounds for t in 1:ntri
        a = SVector{3, Float64}(tris[:, 3t - 2]); b = SVector{3, Float64}(tris[:, 3t - 1]); c3 = SVector{3, Float64}(tris[:, 3t])
        areas[t] = 0.5 * norm(cross(b - a, c3 - a))
    end
    total = sum(areas)
    total > 0 || throw(ArgumentError("$(basename(path)) has zero surface area."))
    cdf = cumsum(areas) ./ total
    out = Matrix{Float64}(undef, 3, Int(n_points))
    @inbounds for k in 1:Int(n_points)
        t = clamp(searchsortedfirst(cdf, rand(rng)), 1, ntri)
        a = SVector{3, Float64}(tris[:, 3t - 2]); b = SVector{3, Float64}(tris[:, 3t - 1]); c3 = SVector{3, Float64}(tris[:, 3t])
        u, v = rand(rng), rand(rng)
        if u + v > 1.0
            u, v = 1.0 - u, 1.0 - v
        end
        p = a + u * (b - a) + v * (c3 - a)
        out[1, k] = p[1]; out[2, k] = p[2]; out[3, k] = p[3]
    end
    return out
end
