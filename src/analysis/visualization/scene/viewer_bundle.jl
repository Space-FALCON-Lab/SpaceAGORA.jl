# Single-file HTML bundler for the three.js viewer under top-level `viewer/`.
#
# The page carries everything: the vendored three.js modules and our own ES
# modules through an import map of `data:` URLs (so the browser resolves
# `import ... from 'three'` with no network and no server), the planet texture
# as a data URI, the scene sidecar as JSON, and the decimated trajectory as
# base64 little-endian Float32 blocks. Opening the file from disk is enough.

const VIEWER_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "viewer"))
const TEXTURES_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "data", "textures"))
const VIEWER_MODULES = ("data.js", "colormaps.js", "globe.js", "atmosphere.js", "spacecraft.js", "lod.js", "ensemble.js", "paths.js", "references.js", "video.js", "plots.js", "terrain.js", "plumes.js", "timeline.js", "ui.js", "main.js")
const VIEWER_VENDOR = (
    "three" => joinpath("vendor", "three.module.js"),
    "three/addons/controls/OrbitControls.js" => joinpath("vendor", "OrbitControls.js"),
    "three/addons/loaders/STLLoader.js" => joinpath("vendor", "STLLoader.js"),
    "three/addons/loaders/OBJLoader.js" => joinpath("vendor", "OBJLoader.js"),
    "three/addons/loaders/GLTFLoader.js" => joinpath("vendor", "GLTFLoader.js"),
    "three/addons/utils/BufferGeometryUtils.js" => joinpath("vendor", "BufferGeometryUtils.js"),
    "mp4-muxer" => joinpath("vendor", "mp4-muxer.mjs"),
)
# Runs up to this many spacecraft embed Float64 positions so a 3 m assembly
# does not jitter at planetary distances; larger runs keep Float32.
const FLOAT64_POSITION_MAX_SPACECRAFT = 64
const DEFAULT_MAX_FRAMES = 2000
const DEFAULT_DATA_BUDGET_MB = 150.0

# ---------------------------------------------------------------------------
# Textures
# ---------------------------------------------------------------------------

"""
    texture_manifest(dir=TEXTURES_DIR) -> Dict{String, Dict{String, Dict{String, Any}}}

Entries of `data/textures/manifest.toml` keyed by body (lowercase planet
name) and then by resolution tier (`"4k"`, `"8k"`, ...; an entry without a
`resolution` field is treated as `"4k"`), each with an absolute `path` added.
Empty when the manifest is absent.
"""
function texture_manifest(dir::AbstractString=TEXTURES_DIR)::Dict{String, Dict{String, Dict{String, Any}}}
    manifest_path = joinpath(dir, "manifest.toml")
    out = Dict{String, Dict{String, Dict{String, Any}}}()
    isfile(manifest_path) || return out
    parsed = TOML.parsefile(manifest_path)
    for entry in get(parsed, "texture", Any[])
        body = lowercase(String(entry["body"]))
        tier = lowercase(String(get(entry, "resolution", "4k")))
        e = Dict{String, Any}(entry)
        e["resolution"] = tier
        e["path"] = joinpath(dir, String(entry["file"]))
        get!(out, body, Dict{String, Dict{String, Any}}())[tier] = e
    end
    return out
end

@inline _tier_pixels(tier::AbstractString)::Int = something(tryparse(Int, replace(lowercase(String(tier)), "k" => "")), 0) * 1024

"""
    texture_entry(body_key; resolution=:best, dir=TEXTURES_DIR) -> Union{Nothing, Dict}

The manifest entry for a body: the largest tier whose file exists for
`:best`, the exact tier for a string or symbol such as `"8k"` (falling back
to the largest available tier at or below it), or `nothing`.
"""
function texture_entry(body_key::AbstractString; resolution=:best, dir::AbstractString=TEXTURES_DIR)
    tiers = get(texture_manifest(dir), lowercase(String(body_key)), nothing)
    tiers === nothing && return nothing
    present = [e for (_, e) in tiers if isfile(e["path"])]
    isempty(present) && return nothing
    sort!(present; by=e -> _tier_pixels(e["resolution"]))
    resolution === :best && return present[end]
    wanted = _tier_pixels(String(resolution))
    wanted > 0 || throw(ArgumentError("texture resolution must be :best or a tier such as \"4k\" or \"8k\", got $(repr(resolution))."))
    at_or_below = [e for e in present if _tier_pixels(e["resolution"]) <= wanted]
    return isempty(at_or_below) ? present[1] : at_or_below[end]
end

@inline function _mime_for(path::AbstractString)::String
    ext = lowercase(splitext(path)[2])
    ext in (".jpg", ".jpeg") && return "image/jpeg"
    ext == ".png" && return "image/png"
    ext == ".webp" && return "image/webp"
    throw(ArgumentError("Unsupported texture type $(ext) for $(path); use jpg, png or webp."))
end

@inline _data_url(bytes::Vector{UInt8}, mime::AbstractString)::String = string("data:", mime, ";base64,", base64encode(bytes))
@inline _js_data_url(source::AbstractString)::String = _data_url(Vector{UInt8}(codeunits(String(source))), "text/javascript")

"""
    texture_payload(planet_key; resolution=:best, dir=TEXTURES_DIR) -> Union{Nothing, Dict}

`{"url" => data URI, "lon_left_deg" => ..., "resolution" => ..., "width" => ...}`
for the body at the chosen tier (see `texture_entry`), or `nothing` when no
texture is registered (the viewer then draws a flat color).
"""
function texture_payload(planet_key::AbstractString; resolution=:best, dir::AbstractString=TEXTURES_DIR)
    entry = texture_entry(planet_key; resolution=resolution, dir=dir)
    entry === nothing && return nothing
    path = entry["path"]
    return Dict{String, Any}(
        "url" => _data_url(read(path), _mime_for(path)),
        "lon_left_deg" => Float64(get(entry, "lon_left_deg", -180.0)),
        "resolution" => String(entry["resolution"]),
        "width" => Int(get(entry, "width", _tier_pixels(entry["resolution"]))),
        "height" => Int(get(entry, "height", _tier_pixels(entry["resolution"]) ÷ 2)),
        "source" => String(get(entry, "source", "")),
        "license" => String(get(entry, "license", "")),
    )
end

# ---------------------------------------------------------------------------
# Frame payload
# ---------------------------------------------------------------------------

@inline function _float32_base64(values::AbstractVector{<:Real})::String
    buf = Vector{Float32}(undef, length(values))
    @inbounds for i in eachindex(values)
        buf[i] = Float32(values[i])
    end
    ENDIAN_BOM == 0x04030201 || (buf .= bswap.(buf))
    return base64encode(reinterpret(UInt8, buf))
end

@inline function _float64_base64(values::AbstractVector{<:Real})::String
    buf = Vector{Float64}(undef, length(values))
    @inbounds for i in eachindex(values)
        buf[i] = Float64(values[i])
    end
    ENDIAN_BOM == 0x04030201 || (buf .= bswap.(buf))
    return base64encode(reinterpret(UInt8, buf))
end

@inline _has_columns(df, cols) = all(c -> c in names(df), cols)

"""
    kept_row_indices(n_rows, stride) -> Vector{Int}

Every `stride`-th row from the first, always including the last row.
"""
function kept_row_indices(n_rows::Integer, stride::Integer)::Vector{Int}
    n_rows <= 0 && return Int[]
    idx = collect(1:Int(stride):Int(n_rows))
    idx[end] == n_rows || push!(idx, Int(n_rows))
    return idx
end

"""
    build_viewer_frames(df, scene; max_frames=2000, data_budget_mb=150.0) -> Dict{String, Any}

Decimate the results table to the playback budget and encode positions (km),
velocities (km/s), attitude quaternions and link poses as base64 Float32
blocks, row-major over (frame, spacecraft, component).
"""
function build_viewer_frames(
    df::DataFrame,
    scene::VisualizationScene;
    max_frames::Integer=DEFAULT_MAX_FRAMES,
    data_budget_mb::Real=DEFAULT_DATA_BUDGET_MB
)::Dict{String, Any}
    S = length(scene.spacecraft)
    S >= 1 || throw(ArgumentError("The scene has no spacecraft."))
    n_rows = nrow(df)
    "time" in names(df) || throw(ArgumentError("Results table has no `time` column."))
    for i in 1:S
        _has_columns(df, ("sc$(i)_pos_1", "sc$(i)_pos_2", "sc$(i)_pos_3")) ||
            throw(ArgumentError("Results table has no position columns for spacecraft $(i) (sc$(i)_pos_1..3)."))
    end
    has_vel = all(i -> _has_columns(df, ("sc$(i)_vel_1", "sc$(i)_vel_2", "sc$(i)_vel_3")), 1:S)
    has_q = all(i -> _has_columns(df, ("sc$(i)_q_1", "sc$(i)_q_2", "sc$(i)_q_3", "sc$(i)_q_4")), 1:S)
    has_mass = all(i -> "sc$(i)_mass" in names(df), 1:S)
    has_density = all(i -> "sc$(i)_density" in names(df), 1:S)
    has_heat = all(i -> "sc$(i)_heat_rate" in names(df), 1:S)
    has_drag = all(i -> _has_columns(df, ("sc$(i)_drag_1", "sc$(i)_drag_2", "sc$(i)_drag_3")), 1:S)
    has_wind = all(i -> _has_columns(df, ("sc$(i)_wind_1", "sc$(i)_wind_2", "sc$(i)_wind_3")), 1:S)
    pos_f64 = S <= FLOAT64_POSITION_MAX_SPACECRAFT
    stride_lp = scene.link_pose_stride
    counts = Int[max(0, length(sc.links) - 1) for sc in scene.spacecraft]
    has_lp = any(>(0), counts) && all(1:S) do i
        counts[i] == 0 || _has_columns(df, ["sc$(i)_$(scene.link_pose_field)_$(k)" for k in 1:(stride_lp * counts[i])])
    end
    lp_offsets = Int[]
    lp_total = 0
    for i in 1:S
        push!(lp_offsets, lp_total)
        lp_total += has_lp ? stride_lp * counts[i] : 0
    end
    # Thruster firing levels (0 to 1), one value per thruster in
    # `scene.spacecraft[i].thrusters` order; absent when the run saved none.
    # A spacecraft whose control effectors report no levels has no columns and
    # counts zero, so one vehicle without them does not drop the whole block.
    thruster_counts = zeros(Int, S)
    for i in 1:S
        n = length(scene.spacecraft[i].thrusters)
        n > 0 || continue
        _has_columns(df, ["sc$(i)_thruster_level_$(k)" for k in 1:n]) || continue
        thruster_counts[i] = n
    end
    thr_total = sum(thruster_counts)
    thr_offsets = Int[]
    let acc = 0
        for n in thruster_counts
            push!(thr_offsets, acc)
            acc += n
        end
    end
    arm_counts = Int[sc.arm === nothing ? 0 : length(sc.arm.links) for sc in scene.spacecraft]
    has_arm = any(>(0), arm_counts) && all(1:S) do i
        arm_counts[i] == 0 || _has_columns(df, ["sc$(i)_arm_pose_$(k)" for k in 1:(stride_lp * arm_counts[i])])
    end
    arm_offsets = Int[]
    arm_total = 0
    for i in 1:S
        push!(arm_offsets, arm_total)
        arm_total += has_arm ? stride_lp * arm_counts[i] : 0
    end

    bytes_per_frame = S * ((pos_f64 ? 24 : 12) + (has_vel ? 12 : 0) + (has_q ? 16 : 0) + (has_mass ? 4 : 0) +
                           (has_density ? 4 : 0) + (has_heat ? 4 : 0) + (has_drag ? 4 : 0) + (has_wind ? 12 : 0)) + 4 * lp_total + 4 * arm_total + 4 * thr_total + 8
    budget = visualization_frame_budget(n_rows, S; max_frames=max_frames, data_budget_mb=data_budget_mb,
                                        bytes_per_sat_frame=cld(bytes_per_frame, S))
    rows = kept_row_indices(n_rows, budget.stride)
    N = length(rows)

    times = Float64[Float64(df[r, "time"]) for r in rows]
    pos = Vector{Float64}(undef, N * S * 3)
    vel = has_vel ? Vector{Float64}(undef, N * S * 3) : Float64[]
    q = has_q ? Vector{Float64}(undef, N * S * 4) : Float64[]
    mass = has_mass ? Vector{Float64}(undef, N * S) : Float64[]
    density = has_density ? Vector{Float64}(undef, N * S) : Float64[]
    heat = has_heat ? Vector{Float64}(undef, N * S) : Float64[]
    drag = has_drag ? Vector{Float64}(undef, N * S) : Float64[]
    wind = has_wind ? Vector{Float64}(undef, N * S * 3) : Float64[]
    lp = has_lp ? Vector{Float64}(undef, N * lp_total) : Float64[]
    ap = has_arm ? Vector{Float64}(undef, N * arm_total) : Float64[]
    thr = thr_total > 0 ? Vector{Float64}(undef, N * thr_total) : Float64[]
    for i in 1:S
        pcols = [df[!, "sc$(i)_pos_$(c)"] for c in 1:3]
        acols = (has_arm && arm_counts[i] > 0) ? [df[!, "sc$(i)_arm_pose_$(k)"] for k in 1:(stride_lp * arm_counts[i])] : nothing
        vcols = has_vel ? [df[!, "sc$(i)_vel_$(c)"] for c in 1:3] : nothing
        qcols = has_q ? [df[!, "sc$(i)_q_$(c)"] for c in 1:4] : nothing
        mcol = has_mass ? df[!, "sc$(i)_mass"] : nothing
        dcol = has_density ? df[!, "sc$(i)_density"] : nothing
        hcol = has_heat ? df[!, "sc$(i)_heat_rate"] : nothing
        fcols = has_drag ? [df[!, "sc$(i)_drag_$(c)"] for c in 1:3] : nothing
        wcols = has_wind ? [df[!, "sc$(i)_wind_$(c)"] for c in 1:3] : nothing
        lcols = (has_lp && counts[i] > 0) ? [df[!, "sc$(i)_$(scene.link_pose_field)_$(k)"] for k in 1:(stride_lp * counts[i])] : nothing
        tcols = thruster_counts[i] > 0 ? [df[!, "sc$(i)_thruster_level_$(k)"] for k in 1:thruster_counts[i]] : nothing
        @inbounds for (f, r) in enumerate(rows)
            base3 = ((f - 1) * S + (i - 1)) * 3
            for c in 1:3
                pos[base3 + c] = Float64(pcols[c][r]) / 1000.0
            end
            if vcols !== nothing
                for c in 1:3
                    vel[base3 + c] = Float64(vcols[c][r]) / 1000.0
                end
            end
            if qcols !== nothing
                base4 = ((f - 1) * S + (i - 1)) * 4
                for c in 1:4
                    q[base4 + c] = Float64(qcols[c][r])
                end
            end
            if mcol !== nothing
                mass[(f - 1) * S + i] = Float64(mcol[r])
            end
            if dcol !== nothing
                density[(f - 1) * S + i] = Float64(dcol[r])
            end
            if hcol !== nothing
                heat[(f - 1) * S + i] = Float64(hcol[r])
            end
            if fcols !== nothing
                drag[(f - 1) * S + i] = sqrt(Float64(fcols[1][r])^2 + Float64(fcols[2][r])^2 + Float64(fcols[3][r])^2)
            end
            if wcols !== nothing
                for c in 1:3
                    wind[base3 + c] = Float64(wcols[c][r])
                end
            end
            if lcols !== nothing
                basel = (f - 1) * lp_total + lp_offsets[i]
                for k in eachindex(lcols)
                    lp[basel + k] = Float64(lcols[k][r])
                end
            end
            if acols !== nothing
                basea = (f - 1) * arm_total + arm_offsets[i]
                for k in eachindex(acols)
                    ap[basea + k] = Float64(acols[k][r])
                end
            end
            if tcols !== nothing
                baset = (f - 1) * thr_total + thr_offsets[i]
                for k in eachindex(tcols)
                    thr[baset + k] = clamp(Float64(tcols[k][r]), 0.0, 1.0)
                end
            end
        end
    end

    frames = Dict{String, Any}(
        "count" => N,
        "sats" => S,
        "source_rows" => n_rows,
        "stride_rows" => budget.stride,
        "t_dtype" => "f64",
        "t_s" => _float64_base64(times),
        "pos_dtype" => pos_f64 ? "f64" : "f32",
        "pos_km" => pos_f64 ? _float64_base64(pos) : _float32_base64(pos),
        "vel_kms" => has_vel ? _float32_base64(vel) : nothing,
        "q" => has_q ? _float32_base64(q) : nothing,
        "mass_kg" => has_mass ? _float32_base64(mass) : nothing,
        "density_kg_m3" => has_density ? _float32_base64(density) : nothing,
        "heat_rate_w_m2" => has_heat ? _float32_base64(heat) : nothing,
        "drag_n" => has_drag ? _float32_base64(drag) : nothing,
        "wind_ms" => has_wind ? _float32_base64(wind) : nothing,
        "link_pose" => has_lp ? Dict{String, Any}(
            "stride" => stride_lp, "counts" => counts, "offsets" => lp_offsets, "total" => lp_total,
            "data" => _float32_base64(lp)
        ) : nothing,
        "thruster_level" => thr_total > 0 ? _float32_base64(thr) : nothing,
        "thruster_counts" => thr_total > 0 ? thruster_counts : nothing,
        "arm_pose" => has_arm ? Dict{String, Any}(
            "stride" => stride_lp, "counts" => arm_counts, "offsets" => arm_offsets, "total" => arm_total,
            "data" => _float32_base64(ap)
        ) : nothing,
    )
    return frames
end

# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------

"""
    viewer_import_map(viewer_dir=VIEWER_DIR) -> Dict{String, Any}

Import map whose entries are `data:` URLs of the vendored three.js modules
and the viewer's own modules, so the bundled page resolves every import
offline.
"""
function viewer_import_map(viewer_dir::AbstractString=VIEWER_DIR)::Dict{String, Any}
    imports = Dict{String, String}()
    for (specifier, rel) in VIEWER_VENDOR
        path = joinpath(viewer_dir, rel)
        isfile(path) || throw(ArgumentError("Missing vendored viewer module $(path)."))
        imports[specifier] = _js_data_url(read(path, String))
    end
    for name in VIEWER_MODULES
        path = joinpath(viewer_dir, "src", name)
        isfile(path) || throw(ArgumentError("Missing viewer module $(path)."))
        imports["viewer/" * name] = _js_data_url(read(path, String))
    end
    return Dict{String, Any}("imports" => imports)
end

# JSON embedded in a <script> must never contain the sequence "</" (it would
# close the tag); "<\/" is the same JSON string.
@inline _script_safe_json(value)::String = replace(JSON.json(value), "</" => "<\\/")

@inline _per_id(value, id::Int, default)::Float64 = value isa AbstractDict ? Float64(get(value, id, default)) : Float64(value)

function _check_gltf_supported(path::AbstractString)
    required = gltf_required_extensions(path)
    bad = [e for e in required if e in GLTF_UNSUPPORTED_REQUIRED]
    isempty(bad) || throw(ArgumentError(
        "$(basename(path)) requires $(join(bad, ", ")), which the viewer cannot decode. " *
        "Re-export it without compression, e.g. `npx @gltf-transform/cli copy in.glb out.glb` " *
        "(decodes Draco on read) or Blender's glTF export with compression off."))
    return nothing
end

"""
    model_payloads(scene; models=Dict(), model_scale=1.0, model_rotation_deg=Dict(), stl=Dict(), stl_scale=1.0) -> Dict{String, Any}

3D model overrides keyed by spacecraft id. `models` maps an id to an STL,
OBJ, glTF or GLB file and wins over the `stl_path` recorded in the scene
(`stl` is an older spelling of the same mapping). `model_scale` is meters per
model unit, a number for all or a `Dict` per id (`stl_scale` applies when the
id is absent). `model_rotation_deg` maps an id to XYZ Euler angles in degrees
applied to the model in the body frame. `model_center` (a Bool or a per-id
`Dict`, default true) shifts the model so its bounding-box center sits on
the spacecraft. Each entry is
`{"url" => data URI, "format" => ..., "scale" => ..., "rotation_deg" => [rx, ry, rz], "center" => [cx, cy, cz] (model units), "source" => file name}`.
A `.gltf` file must embed its buffers; external files are not carried along.
"""
function model_payloads(
    scene::VisualizationScene;
    models::AbstractDict=Dict{Int, String}(),
    model_scale=nothing,
    model_rotation_deg::AbstractDict=Dict{Int, Any}(),
    model_center=true,
    model_articulations::AbstractDict=Dict{Int, Any}(),
    stl::AbstractDict=Dict{Int, String}(),
    stl_scale::Real=1.0
)::Dict{String, Any}
    out = Dict{String, Any}()
    for sc in scene.spacecraft
        path = haskey(models, sc.id) ? String(models[sc.id]) : (haskey(stl, sc.id) ? String(stl[sc.id]) : sc.stl_path)
        path === nothing && continue
        isfile(path) || throw(ArgumentError("3D model for spacecraft $(sc.id) not found: $(path)"))
        format, mime = model_format(path)
        format in ("glb", "gltf") && _check_gltf_supported(path)
        scale = model_scale === nothing ? Float64(stl_scale) : _per_id(model_scale, sc.id, stl_scale)
        rot = get(model_rotation_deg, sc.id, (0.0, 0.0, 0.0))
        length(rot) == 3 || throw(ArgumentError("model_rotation_deg entries must be three angles (rx, ry, rz) in degrees."))
        centerd = model_center isa AbstractDict ? Bool(get(model_center, sc.id, true)) : Bool(model_center)
        articulations = get(model_articulations, sc.id, ())
        center = [0.0, 0.0, 0.0]
        articulation_dicts = Dict{String, Any}[]
        if centerd || !isempty(articulations)
            try
                raw = load_model_triangles(path)
                articulation_dicts = articulation_payload(articulations, raw)
                posed = isempty(articulations) ? raw : articulate_triangles(raw, articulations)
                centerd && (center = [0.5 * (minimum(posed[c, :]) + maximum(posed[c, :])) for c in 1:3])
            catch err
                isempty(articulations) && @warn "Could not read $(basename(path)) to center it; the viewer will use the file's own origin." exception=(err, catch_backtrace())
                isempty(articulations) || rethrow()
            end
        end
        out[string(sc.id)] = Dict{String, Any}(
            "url" => _data_url(read(path), mime),
            "format" => format,
            "scale" => scale,
            "rotation_deg" => Float64[Float64(r) for r in rot],
            "center" => center,
            "articulations" => articulation_dicts,
            "source" => basename(path),
        )
    end
    return out
end

"""
    viewer_payload(scene, df; textures_dir=TEXTURES_DIR, options=Dict(), max_frames=2000, data_budget_mb=150.0, stl=Dict(), stl_scale=1.0) -> Dict{String, Any}

The `window.SPACEAGORA_VIEWER` object: sidecar, encoded frames, the planet's
texture (if registered), STL models (if any) and viewer options.
"""
function viewer_payload(
    scene::VisualizationScene,
    df::DataFrame;
    textures_dir::AbstractString=TEXTURES_DIR,
    include_textures::Bool=true,
    texture_resolution=:best,
    options::AbstractDict=Dict{String, Any}(),
    max_frames::Integer=DEFAULT_MAX_FRAMES,
    data_budget_mb::Real=DEFAULT_DATA_BUDGET_MB,
    models::AbstractDict=Dict{Int, String}(),
    model_scale=nothing,
    model_rotation_deg::AbstractDict=Dict{Int, Any}(),
    model_center=true,
    model_articulations::AbstractDict=Dict{Int, Any}(),
    stl::AbstractDict=Dict{Int, String}(),
    stl_scale::Real=1.0,
    paths=(),
    references=(),
    terrain=nothing
)::Dict{String, Any}
    textures = Dict{String, Any}()
    if include_textures
        entry = texture_payload(scene.planet.texture; resolution=texture_resolution, dir=textures_dir)
        entry === nothing || (textures[scene.planet.texture] = entry)
    end
    return Dict{String, Any}(
        "scene" => scene_dict(scene),
        "frames" => build_viewer_frames(df, scene; max_frames=max_frames, data_budget_mb=data_budget_mb),
        "textures" => textures,
        "models" => model_payloads(scene; models=models, model_scale=model_scale, model_rotation_deg=model_rotation_deg, model_center=model_center, model_articulations=model_articulations, stl=stl, stl_scale=stl_scale),
        "paths" => path_payloads(paths),
        "references" => reference_payloads(references, scene),
        "terrain" => terrain === nothing ? nothing : terrain_payload(terrain),
        "options" => Dict{String, Any}(options),
    )
end

"""
    terrain_payload(site_json; max_grid=512) -> Dict{String, Any}

Site terrain for the page: the DEM grids of a site directory written by
`scripts/dev/terrain/fetch_moon_site.py` (finest first, each subsampled to
at most `max_grid` samples per side, heights as base64 Float32) and its
imagery levels (JPEG data URLs with their latitude/longitude boxes). The
viewer drapes the imagery over the displaced grids and cuts the globe open
under the outermost level.
"""
function terrain_payload(site_json::AbstractString; max_grid::Integer=512)::Dict{String, Any}
    isfile(site_json) || throw(ArgumentError("terrain site file not found: $(site_json)"))
    meta = JSON.parsefile(String(site_json))
    dir = dirname(String(site_json))
    grids = Dict{String, Any}[]
    radius = 0.0
    for d in meta["dem"]
        g = TerrainModels.load_dem_grid(joinpath(dir, String(d["name"]) * ".json"))
        radius = Float64(get(d, "reference_radius_m", 1737400.0))
        rows, cols = size(g.heights)
        stride = max(1, ceil(Int, max(rows, cols) / Int(max_grid)))
        sub = g.heights[1:stride:end, 1:stride:end]
        # the subsampled grid keeps the same outer edges only approximately; state the edges it does cover
        r2, c2 = size(sub)
        dlat = (g.lat_max - g.lat_min) / rows; dlon = (g.lon_max - g.lon_min) / cols
        push!(grids, Dict{String, Any}(
            "name" => String(d["name"]), "rows" => r2, "cols" => c2,
            "lat_max" => g.lat_max, "lat_min" => g.lat_max - r2 * stride * dlat,
            "lon_min" => g.lon_min, "lon_max" => g.lon_min + c2 * stride * dlon,
            "heights" => _float32_base64(vec(permutedims(sub))),
            "source" => g.source,
        ))
    end
    levels = Dict{String, Any}[]
    imagery_rel = get(meta, "imagery", nothing)
    if imagery_rel !== nothing && isfile(joinpath(dir, String(imagery_rel)))
        im = JSON.parsefile(joinpath(dir, String(imagery_rel)))
        for lvl in im["levels"]
            path = joinpath(dir, dirname(String(imagery_rel)), String(lvl["file"]))
            isfile(path) || continue
            push!(levels, Dict{String, Any}(
                "lat_min" => lvl["lat_min"], "lat_max" => lvl["lat_max"], "lon_min" => lvl["lon_min"], "lon_max" => lvl["lon_max"],
                "width" => lvl["width"], "height" => lvl["height"], "m_per_px" => lvl["m_per_px"],
                "url" => _data_url(read(path), "image/jpeg"),
            ))
        end
    end
    site = meta["site"]
    model, info = TerrainModels.load_site_terrain(String(site_json))
    return Dict{String, Any}(
        "site" => Dict{String, Any}("lat_deg" => site["lat_deg"], "lon_deg" => site["lon_deg"], "name" => get(site, "name", "site"), "height_m" => info.height_m),
        "reference_radius_m" => radius, "grids" => grids, "imagery" => levels,
    )
end

"""
    path_payloads(paths) -> Vector{Dict{String, Any}}

Reference polylines drawn beside the flown trajectories. Each entry of
`paths` is a NamedTuple or Dict with `name`, `points_m` (3 x N, meters),
`frame` (`:inertial`, or `:rtn` for the radial/transverse/normal frame of
spacecraft `target`, 1-based index into the run's spacecraft list, or
`:body` for that spacecraft's body frame) and optionally `color` (hex
string) and `dashed` (Bool).
"""
function path_payloads(paths)::Vector{Dict{String, Any}}
    out = Dict{String, Any}[]
    for (k, path) in enumerate(paths)
        get_ = (key, default) -> path isa AbstractDict ? get(path, key, get(path, String(key), default)) : (hasproperty(path, key) ? getproperty(path, key) : default)
        pts = get_(:points_m, nothing)
        pts === nothing && throw(ArgumentError("path $(k) needs points_m (3 x N, meters)."))
        M = Matrix{Float64}(pts)
        size(M, 1) == 3 || throw(ArgumentError("path $(k): points_m must be 3 x N."))
        frame = Symbol(get_(:frame, :inertial))
        frame in (:inertial, :rtn, :body) || throw(ArgumentError("path $(k): frame must be :inertial, :rtn or :body."))
        push!(out, Dict{String, Any}(
            "name" => String(get_(:name, "path $(k)")),
            "frame" => String(frame),
            "target" => Int(get_(:target, 1)),
            "color" => String(get_(:color, "#7fe0ff")),
            "dashed" => Bool(get_(:dashed, true)),
            "points_km" => _float32_base64(vec(M) ./ 1000.0),
            "count" => size(M, 2),
        ))
    end
    return out
end

"""
    reference_payloads(references, scene) -> Vector{Dict{String, Any}}

Reference trajectories drawn as translucent ghosts of a spacecraft on the
run's own timeline: a SPICE reconstruction, a telemetry record or a plan.
Each entry of `references` is a NamedTuple or Dict with `name`, `t_s`
(elapsed seconds from the run epoch, length N, non-decreasing), `pos_m`
(3 x N, inertial meters) and optionally `vel_mps` (3 x N), `q` (4 x N,
scalar-last body-to-inertial attitude; velocity-aligned when absent),
`target` (1-based index of the spacecraft whose geometry and 3D model the
ghost copies, default 1), `color` (hex string), `opacity` (0..1, default
0.45) and `trail` (draw the whole reference line, default true). Times and
positions are embedded as Float64 so a ghost sits within meters of the flown
spacecraft when the two agree.
"""
function reference_payloads(references, scene::VisualizationScene)::Vector{Dict{String, Any}}
    out = Dict{String, Any}[]
    n_sc = length(scene.spacecraft)
    for (k, ref) in enumerate(references)
        get_ = (key, default) -> ref isa AbstractDict ? get(ref, key, get(ref, String(key), default)) : (hasproperty(ref, key) ? getproperty(ref, key) : default)
        t = get_(:t_s, nothing)
        t === nothing && throw(ArgumentError("reference $(k) needs t_s (elapsed seconds from the run epoch)."))
        times = Float64[Float64(x) for x in vec(collect(t))]
        N = length(times)
        N >= 1 || throw(ArgumentError("reference $(k): t_s is empty."))
        all(isfinite, times) && issorted(times) || throw(ArgumentError("reference $(k): t_s must be finite and non-decreasing."))
        pts = get_(:pos_m, nothing)
        pts === nothing && throw(ArgumentError("reference $(k) needs pos_m (3 x N, meters)."))
        P = Matrix{Float64}(pts)
        size(P) == (3, N) || throw(ArgumentError("reference $(k): pos_m must be 3 x $(N) to match t_s, got $(size(P))."))
        all(isfinite, P) || throw(ArgumentError("reference $(k): pos_m has non-finite entries."))
        vel = get_(:vel_mps, nothing)
        V = vel === nothing ? nothing : Matrix{Float64}(vel)
        V === nothing || size(V) == (3, N) || throw(ArgumentError("reference $(k): vel_mps must be 3 x $(N)."))
        qs = get_(:q, nothing)
        Q = qs === nothing ? nothing : Matrix{Float64}(qs)
        Q === nothing || size(Q) == (4, N) || throw(ArgumentError("reference $(k): q must be 4 x $(N) (scalar-last)."))
        target = Int(get_(:target, 1))
        1 <= target <= n_sc || throw(ArgumentError("reference $(k): target must be a spacecraft index in 1:$(n_sc), got $(target)."))
        opacity = Float64(get_(:opacity, 0.45))
        0.0 <= opacity <= 1.0 || throw(ArgumentError("reference $(k): opacity must be within 0..1."))
        push!(out, Dict{String, Any}(
            "name" => String(get_(:name, "reference $(k)")),
            "target" => target,
            "count" => N,
            "t_s" => _float64_base64(times),
            "pos_km" => _float64_base64(vec(P) ./ 1000.0),
            "vel_kms" => V === nothing ? nothing : _float32_base64(vec(V) ./ 1000.0),
            "q" => Q === nothing ? nothing : _float32_base64(vec(Q)),
            "color" => String(get_(:color, "#ff8c69")),
            "opacity" => opacity,
            "trail" => Bool(get_(:trail, true)),
        ))
    end
    return out
end

"""
    render_viewer_html(payload; viewer_dir=VIEWER_DIR, title="SpaceAGORA viewer") -> String

Fill `viewer/template.html` with the import map and the payload.
"""
function render_viewer_html(payload::AbstractDict; viewer_dir::AbstractString=VIEWER_DIR, title::AbstractString="SpaceAGORA viewer")::String
    template_path = joinpath(viewer_dir, "template.html")
    isfile(template_path) || throw(ArgumentError("Missing viewer template $(template_path)."))
    html = read(template_path, String)
    for token in ("__TITLE__", "__IMPORTMAP__", "__PAYLOAD__")
        occursin(token, html) || throw(ArgumentError("Viewer template is missing the $(token) placeholder."))
    end
    safe_title = replace(String(title), "<" => "&lt;", ">" => "&gt;", "&" => "&amp;")
    html = replace(html, "__TITLE__" => safe_title)
    html = replace(html, "__IMPORTMAP__" => _script_safe_json(viewer_import_map(viewer_dir)))
    return replace(html, "__PAYLOAD__" => _script_safe_json(payload))
end

@inline function _viewer_options(; trail_s, trail_orbits, frame, speed, title)
    frame in (:inertial, :planet_fixed) || throw(ArgumentError("frame must be :inertial or :planet_fixed, got $(frame)."))
    options = Dict{String, Any}("frame" => String(frame))
    trail_s === nothing || (options["trail_s"] = Float64(trail_s))
    trail_orbits === nothing || (options["trail_orbits"] = Float64(trail_orbits))
    speed === nothing || (options["speed"] = Float64(speed))
    title === nothing || (options["title"] = String(title))
    return options
end

"""
    export_visualization(prefix; out=prefix * "_viewer.html", kwargs...) -> String

Build the self-contained viewer page for a run written with
`save_visualization_scene = true`. `prefix` is the results bundle prefix
(`<results_directory>/simulation_results`); the sidecar `<prefix>_scene.json`
and the Arrow file `<prefix>.feather` must exist.

Keywords: `max_frames` (default 2000) and `data_budget_mb` (default 150)
bound the embedded trajectory; `trail_orbits` (default 3, estimated from the
trajectory's periapsis passages) or `trail_s` set the trail window; `frame`
is `:inertial` (default) or `:planet_fixed`; `speed` is the initial playback
rate in simulated seconds per wall second; `title` names the page;
`textures=false` skips the surface texture and `texture_resolution` picks a
tier (`:best`, the default, takes the largest registered, e.g. 8k for Earth;
`"4k"` keeps the page small); `models` maps spacecraft ids to STL, OBJ, glTF
or GLB files drawn instead of the link boxes, at `model_scale` meters per
model unit (a number or a per-id `Dict`) and rotated by `model_rotation_deg`
(per-id XYZ Euler angles), posed by `model_articulations` (per-id list of
parts rotated about an axis, see `articulate_triangles`) and centerd on the
spacecraft unless `model_center=false`; `stl`/`stl_scale` are the older spelling for STL only; `references` draws
reference trajectories as translucent ghosts of a spacecraft (see
`reference_payloads`); `paths` overlays
reference polylines (see `path_payloads`), e.g. a planned RPO path in the
target's RTN frame; `viewer_dir` and `textures_dir` override the repository locations.
"""
function export_visualization(
    prefix::AbstractString;
    out::Union{Nothing, AbstractString}=nothing,
    max_frames::Integer=DEFAULT_MAX_FRAMES,
    data_budget_mb::Real=DEFAULT_DATA_BUDGET_MB,
    trail_s::Union{Nothing, Real}=nothing,
    trail_orbits::Union{Nothing, Real}=nothing,
    frame::Symbol=:inertial,
    speed::Union{Nothing, Real}=nothing,
    title::Union{Nothing, AbstractString}=nothing,
    textures::Bool=true,
    texture_resolution=:best,
    models::AbstractDict=Dict{Int, String}(),
    model_scale=nothing,
    model_rotation_deg::AbstractDict=Dict{Int, Any}(),
    model_center=true,
    model_articulations::AbstractDict=Dict{Int, Any}(),
    stl::AbstractDict=Dict{Int, String}(),
    stl_scale::Real=1.0,
    paths=(),
    references=(),
    terrain=nothing,
    viewer_dir::AbstractString=VIEWER_DIR,
    textures_dir::AbstractString=TEXTURES_DIR
)::String
    prefix = String(prefix)
    scene_path = prefix * "_scene.json"
    isfile(scene_path) || throw(ArgumentError("No visualization sidecar at $(scene_path); run with simulation_settings.save_visualization_scene = true."))
    scene = read_visualization_scene(scene_path)
    feather_path = joinpath(dirname(scene_path), scene.results_feather)
    isfile(feather_path) || throw(ArgumentError("Results file $(feather_path) referenced by the sidecar does not exist."))
    df = DataFrame(Arrow.Table(feather_path))
    payload = viewer_payload(
        scene, df;
        textures_dir=textures_dir, include_textures=textures, texture_resolution=texture_resolution,
        options=_viewer_options(; trail_s=trail_s, trail_orbits=trail_orbits, frame=frame, speed=speed, title=title),
        max_frames=max_frames, data_budget_mb=data_budget_mb,
        models=models, model_scale=model_scale, model_rotation_deg=model_rotation_deg, model_center=model_center, model_articulations=model_articulations, stl=stl, stl_scale=stl_scale,
        paths=paths, references=references, terrain=terrain
    )
    page_title = title === nothing ? "SpaceAGORA · $(scene.planet.name) · $(basename(prefix))" : String(title)
    html = render_viewer_html(payload; viewer_dir=viewer_dir, title=page_title)
    out_path = out === nothing ? prefix * "_viewer.html" : String(out)
    return IOSerialization._atomic_write_file(out_path, tmp -> write(tmp, html))
end

"""
    export_visualization(args::SimulationConfiguration; kwargs...) -> String

Same as the prefix form, for the results bundle of `args`.
"""
function export_visualization(args::SimulationConfiguration; kwargs...)::String
    return export_visualization(IOConfig._results_bundle_prefix(args); kwargs...)
end

"""
    write_viewer_dev_payload(prefix, out_js; kwargs...) -> String

Write `window.SPACEAGORA_VIEWER = {...};` for `viewer/dev.html`, so the
viewer modules can be edited and reloaded without rebundling.
"""
function write_viewer_dev_payload(prefix::AbstractString, out_js::AbstractString; textures::Bool=true, textures_dir::AbstractString=TEXTURES_DIR, kwargs...)::String
    scene_path = String(prefix) * "_scene.json"
    scene = read_visualization_scene(scene_path)
    df = DataFrame(Arrow.Table(joinpath(dirname(scene_path), scene.results_feather)))
    payload = viewer_payload(scene, df; textures_dir=textures_dir, include_textures=textures, kwargs...)
    return IOSerialization._atomic_write_file(String(out_js), tmp -> write(tmp, "window.SPACEAGORA_VIEWER = " * _script_safe_json(payload) * ";\n"))
end
