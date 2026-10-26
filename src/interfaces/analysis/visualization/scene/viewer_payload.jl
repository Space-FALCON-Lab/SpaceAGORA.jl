# Renderer-independent saved-data conversion. The optional visualizer owns
# browser modules and HTML assembly; the data contract stays with SpaceAGORA.

# Mission history is optional plain data shared with the independent viewer.
# Validate structure without changing times, names or the underlying results.
function _validated_mission_history(mission)
    mission === nothing && return nothing
    mission isa AbstractDict || throw(ArgumentError("mission must be a JSON object."))
    schema = get(mission, "schema", nothing)
    schema isa Integer && !(schema isa Bool) && schema == 1 ||
        throw(ArgumentError("Unsupported mission history schema; expected 1."))
    name_ok(x) = x isa AbstractString && !isempty(strip(x))
    time_ok(x) = x isa Real && !(x isa Bool) && isfinite(x) && x >= 0
    name_ok(get(mission, "name", nothing)) || throw(ArgumentError("mission.name must be a nonempty string."))
    segments = get(mission, "segments", nothing)
    transitions = get(mission, "transitions", nothing)
    segments isa AbstractVector || throw(ArgumentError("mission.segments must be an array."))
    transitions isa AbstractVector || throw(ArgumentError("mission.transitions must be an array."))
    names = Set{String}()
    for (i, segment) in enumerate(segments)
        segment isa AbstractDict || throw(ArgumentError("mission segment must be an object."))
        name = get(segment, "name", nothing)
        name_ok(name) && !(name in names) || throw(ArgumentError("Mission segment names must be nonempty and unique."))
        push!(names, name)
        start_s, end_s = get(segment, "start_s", nothing), get(segment, "end_s", nothing)
        time_ok(start_s) && time_ok(end_s) && start_s <= end_s ||
            throw(ArgumentError("Mission segment times must be finite and satisfy 0 <= start_s <= end_s."))
        i == 1 || segments[i - 1]["end_s"] == start_s ||
            throw(ArgumentError("Mission segment intervals must be contiguous and ordered."))
    end
    n = length(segments)
    (n == 0 ? isempty(transitions) : length(transitions) in (n - 1, n)) ||
        throw(ArgumentError("Mission history needs one transition per adjacent segment and an optional final exit."))
    for (i, transition) in enumerate(transitions)
        transition isa AbstractDict || throw(ArgumentError("mission transition must be an object."))
        t = get(transition, "time_s", nothing)
        time_ok(t) && t == segments[i]["end_s"] ||
            throw(ArgumentError("Mission transition time must match its segment boundary."))
        get(transition, "from", nothing) == segments[i]["name"] ||
            throw(ArgumentError("Mission transition.from must name its outgoing segment."))
        haskey(transition, "to") && transition["to"] == (i < n ? segments[i + 1]["name"] : nothing) ||
            throw(ArgumentError("Mission transition.to must name the next segment, or be null at final exit."))
        name_ok(get(transition, "reason", nothing)) ||
            throw(ArgumentError("Mission transition.reason must be a nonempty string."))
    end
    return mission
end

const TEXTURES_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "..", "data", "scene", "textures"))
# Plume-surface interaction block of the frames payload: one Float32 array per
# quantity, frame-major then spacecraft, from the `sc{i}_plume_*` result columns.
const PLUME_FRAME_FIELDS = ("height_m", "shear_pa", "pressure_pa", "erosion_kg_s", "eroded_kg", "ejecta_mps", "ground_effect_n")

# Runs up to this many spacecraft embed Float64 positions so a 3 m assembly
# does not jitter at planetary distances; larger runs do too when that costs
# no frames (`_float64_positions`), and keep Float32 otherwise.
const FLOAT64_POSITION_MAX_SPACECRAFT = 64
const DEFAULT_MAX_FRAMES = 2000
const DEFAULT_DATA_BUDGET_MB = 150.0
# Saved `sc{i}_heat_rate` columns are W/cm²; the frames field `heat_rate_w_m2` is W/m².
const HEAT_RATE_W_CM2_TO_W_M2 = 1.0e4

# ---------------------------------------------------------------------------
# Textures
# ---------------------------------------------------------------------------

"""
    texture_manifest(dir=TEXTURES_DIR) -> Dict{String, Dict{String, Dict{String, Any}}}

Entries of `data/scene/textures/manifest.toml` keyed by body (lowercase planet
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
    texture_payload(planet_key; resolution="4k", dir=TEXTURES_DIR) -> Union{Nothing, Dict}

`{"url" => data URI, "lon_left_deg" => ..., "resolution" => ..., "width" => ...}`
for the body at the chosen tier (see `texture_entry`), or `nothing` when no
texture is registered (the viewer then draws a flat color).
"""
function texture_payload(planet_key::AbstractString; resolution="4k", dir::AbstractString=TEXTURES_DIR)
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

# Keep exact saved boundaries while staying inside the requested frame budget.
# A transition without a saved row remains metadata; no state is invented here.
# `hint` ends the refusal's advice.
function _mission_kept_rows(df, mission, regular_rows, limit::Int; hint::AbstractString="")
    mission === nothing && return regular_rows
    isempty(mission["segments"]) && return regular_rows
    nrow(df) == 0 && return Int[]
    boundaries = Set(v for segment in mission["segments"] for v in (segment["start_s"], segment["end_s"]))
    required = unique!([1; findall(t -> t in boundaries, df[!, "time"]); nrow(df)])
    length(required) <= limit || throw(ArgumentError(
        "Mission boundary samples exceed the viewer frame budget; increase max_frames or data_budget_mb$(hint)."))
    available = setdiff(collect(1:nrow(df)), required)
    count = min(limit - length(required), length(available))
    extra = count == 0 ? Int[] : count == 1 ? [available[cld(length(available), 2)]] :
        available[unique(round.(Int, range(1, length(available); length=count)))]
    return sort!([required; extra])
end

# `position_precision` as given, refused unless it is one of the three values;
# the export entry points check it before reading any results.
function _checked_position_precision(position_precision)::Symbol
    position_precision isa Symbol && position_precision in (:auto, :f64, :f32) && return position_precision
    throw(ArgumentError("position_precision must be :auto, :f64 or :f32, got $(repr(position_precision))."))
end

# Whether `S` spacecraft get Float64 positions: always with `:f64`, never with
# `:f32`. With `:auto`, at `FLOAT64_POSITION_MAX_SPACECRAFT` or fewer, and above
# that when Float64 costs no frames: `same_frames()` says whether the frames kept
# at 24 bytes per position are those kept at 12, within the budget.
function _float64_positions(position_precision, S::Integer, same_frames)::Bool
    precision = _checked_position_precision(position_precision)
    precision === :f64 && return true
    precision === :f32 && return false
    return S <= FLOAT64_POSITION_MAX_SPACECRAFT || same_frames()
end

# What decides the frames kept at `bytes_per_frame`: the budget, whose regular
# stride picks the rows, and the row limit that bounds a mission's rows.
function _viewer_frame_plan(n_rows::Integer, S::Integer, bytes_per_frame::Integer, max_frames::Integer, data_budget_mb::Real)
    budget = visualization_frame_budget(n_rows, S; max_frames=max_frames, data_budget_mb=data_budget_mb,
                                        bytes_per_sat_frame=cld(bytes_per_frame, S))
    limit = min(n_rows, Int(max_frames), max(min(n_rows, 2), floor(Int, Float64(data_budget_mb) * 1.0e6 / bytes_per_frame)))
    return (budget=budget, limit=limit)
end

"""
    build_viewer_frames(df, scene; max_frames=2000, data_budget_mb=150.0, position_precision=:auto) -> Dict{String, Any}

Decimate the results table to the playback budget and encode positions (km),
velocities (km/s), attitude quaternions and link poses as base64 Float32
blocks, row-major over (frame, spacecraft, component). Positions are Float64
for 64 spacecraft or fewer, and above that when Float64 keeps the same frames
as Float32 would; `position_precision=:f64` or `:f32` forces either.
"""
function build_viewer_frames(
    df::DataFrame,
    scene::VisualizationScene;
    max_frames::Integer=DEFAULT_MAX_FRAMES,
    data_budget_mb::Real=DEFAULT_DATA_BUDGET_MB,
    channels=(),
    mission=nothing,
    position_precision=:auto
)::Dict{String, Any}
    _checked_position_precision(position_precision)
    mission = _validated_mission_history(mission)
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
    has_plume = all(i -> _has_columns(df, ["sc$(i)_plume_$(f)" for f in PLUME_FRAME_FIELDS]), 1:S)
    # Sun direction: one unit vector per row for the whole scene (the viewer's
    # sun lighting), absent when the run's ephemerides could not resolve the Sun.
    has_sun = _has_columns(df, ("sun_dir_1", "sun_dir_2", "sun_dir_3"))
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
    # Thruster exhaust directions (three per thruster, in its link's frame, in
    # the same order), saved when a run turns a thruster with a gimbal; the
    # frames carry the block only then.
    dir_counts = zeros(Int, S)
    for i in 1:S
        n = length(scene.spacecraft[i].thrusters)
        n > 0 || continue
        _has_columns(df, ["sc$(i)_thruster_direction_$(k)" for k in 1:(3 * n)]) || continue
        dir_counts[i] = n
    end
    dir_total = 3 * sum(dir_counts)
    dir_offsets = Int[]
    let acc = 0
        for n in dir_counts
            push!(dir_offsets, acc)
            acc += 3 * n
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

    channel_specs = _resolved_channels(df, channels, S)
    # Bytes per frame besides the positions: per spacecraft, and for the frame.
    other_sat_bytes = (has_vel ? 12 : 0) + (has_q ? 16 : 0) + (has_mass ? 4 : 0) +
                      (has_density ? 4 : 0) + (has_heat ? 4 : 0) + (has_drag ? 4 : 0) + (has_wind ? 12 : 0) +
                      (has_plume ? 4 * length(PLUME_FRAME_FIELDS) : 0) + 4 * length(channel_specs)
    other_frame_bytes = 4 * lp_total + 4 * arm_total + 4 * thr_total + 4 * dir_total + (has_sun ? 12 : 0) + 8
    plan(position_bytes) = _viewer_frame_plan(n_rows, S, S * (position_bytes + other_sat_bytes) + other_frame_bytes,
        max_frames, data_budget_mb)
    # The kept rows follow the stride, and with mission segments the row limit.
    # Below two frames the budget's floor keeps two at either precision, and the
    # Float64 frames then exceed the budget: that costs, so Float32.
    mission_rows = mission !== nothing && !isempty(mission["segments"])
    same_frames(a, b) = a.budget.stride == b.budget.stride && (!mission_rows || a.limit == b.limit) &&
        a.budget.bytes <= Float64(data_budget_mb) * 1.0e6
    pos_f64 = _float64_positions(position_precision, S, () -> same_frames(plan(24), plan(12)))
    (; budget, limit) = plan(pos_f64 ? 24 : 12)
    regular_rows = kept_row_indices(n_rows, budget.stride)
    forced_f64 = position_precision === :f64 && S > FLOAT64_POSITION_MAX_SPACECRAFT
    rows = _mission_kept_rows(df, mission, regular_rows, limit;
        hint=forced_f64 ? ", or use position_precision=:auto, which can keep Float32 positions to leave room for more rows" : "")
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
    sun = has_sun ? Vector{Float64}(undef, N * 3) : Float64[]
    lp = has_lp ? Vector{Float64}(undef, N * lp_total) : Float64[]
    ap = has_arm ? Vector{Float64}(undef, N * arm_total) : Float64[]
    plume = has_plume ? [Vector{Float64}(undef, N * S) for _ in PLUME_FRAME_FIELDS] : Vector{Float64}[]
    thr = thr_total > 0 ? Vector{Float64}(undef, N * thr_total) : Float64[]
    tdir = dir_total > 0 ? Vector{Float64}(undef, N * dir_total) : Float64[]
    if has_sun
        scols = [df[!, "sun_dir_$(c)"] for c in 1:3]
        @inbounds for (f, r) in enumerate(rows)
            for c in 1:3
                sun[(f - 1) * 3 + c] = Float64(scols[c][r])
            end
        end
    end
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
        plcols = has_plume ? [df[!, "sc$(i)_plume_$(f)"] for f in PLUME_FRAME_FIELDS] : nothing
        lcols = (has_lp && counts[i] > 0) ? [df[!, "sc$(i)_$(scene.link_pose_field)_$(k)"] for k in 1:(stride_lp * counts[i])] : nothing
        tcols = thruster_counts[i] > 0 ? [df[!, "sc$(i)_thruster_level_$(k)"] for k in 1:thruster_counts[i]] : nothing
        dcols = dir_counts[i] > 0 ? [df[!, "sc$(i)_thruster_direction_$(k)"] for k in 1:(3 * dir_counts[i])] : nothing
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
                heat[(f - 1) * S + i] = Float64(hcol[r]) * HEAT_RATE_W_CM2_TO_W_M2
            end
            if fcols !== nothing
                drag[(f - 1) * S + i] = sqrt(Float64(fcols[1][r])^2 + Float64(fcols[2][r])^2 + Float64(fcols[3][r])^2)
            end
            if wcols !== nothing
                for c in 1:3
                    wind[base3 + c] = Float64(wcols[c][r])
                end
            end
            if plcols !== nothing
                for k in eachindex(plcols)
                    plume[k][(f - 1) * S + i] = Float64(plcols[k][r])
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
            if dcols !== nothing
                based = (f - 1) * dir_total + dir_offsets[i]
                for k in eachindex(dcols)
                    tdir[based + k] = Float64(dcols[k][r])
                end
            end
        end
    end

    channel_values = [Vector{Float32}(undef, N * S) for _ in channel_specs]
    for (c, spec) in enumerate(channel_specs)
        cols = [df[!, "sc$(i)_$(spec.column)"] for i in 1:S]
        @inbounds for (f, r) in enumerate(rows), i in 1:S
            channel_values[c][(f - 1) * S + i] = _channel_float32(cols[i][r], spec.column)
        end
    end
    frames = Dict{String, Any}(
        "channels" => isempty(channel_specs) ? nothing : [Dict{String, Any}(
            "name" => spec.column, "label" => spec.label, "unit" => spec.unit,
            "digits" => spec.digits, "log" => spec.log,
            "data" => _float32_base64(channel_values[k])
        ) for (k, spec) in enumerate(channel_specs)],
        "count" => N,
        "sats" => S,
        "source_rows" => n_rows,
        "stride_rows" => rows == regular_rows ? budget.stride : nothing,
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
        "sun_dir" => has_sun ? _float32_base64(sun) : nothing,
        "link_pose" => has_lp ? Dict{String, Any}(
            "stride" => stride_lp, "counts" => counts, "offsets" => lp_offsets, "total" => lp_total,
            "data" => _float32_base64(lp)
        ) : nothing,
        "plume" => has_plume ? Dict{String, Any}(
            String(PLUME_FRAME_FIELDS[k]) => _float32_base64(plume[k]) for k in eachindex(PLUME_FRAME_FIELDS)
        ) : nothing,
        "thruster_level" => thr_total > 0 ? _float32_base64(thr) : nothing,
        "thruster_counts" => thr_total > 0 ? thruster_counts : nothing,
        "arm_pose" => has_arm ? Dict{String, Any}(
            "stride" => stride_lp, "counts" => arm_counts, "offsets" => arm_offsets, "total" => arm_total,
            "data" => _float32_base64(ap)
        ) : nothing,
    )
    if dir_total > 0
        frames["thruster_direction"] = Dict{String, Any}(
            "stride" => 3, "counts" => dir_counts, "offsets" => dir_offsets, "total" => dir_total,
            "data" => _float32_base64(tdir))
    end
    return frames
end

# Missing samples and NaN are plot gaps; infinities and Float32 overflow are errors.
function _channel_float32(value, column::AbstractString)::Float32
    value === missing && return NaN32
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("Channel $(column) samples must be real numbers or missing."))
    x = try
        Float32(value)
    catch
        throw(ArgumentError("Channel $(column) sample cannot be represented as Float32."))
    end
    isinf(x) && throw(ArgumentError("Channel $(column) samples must not be infinite or overflow Float32."))
    return x
end

function _resolved_channels(df::DataFrame, channels, S::Integer)
    channels isa Union{Tuple, AbstractVector} || throw(ArgumentError("channels must be a tuple or vector of specifications."))
    specs = NamedTuple[]
    seen = Set{String}()
    for spec in channels
        spec isa Union{NamedTuple, AbstractDict} || throw(ArgumentError("Each channel must be a NamedTuple or Dict."))
        field(key, default) = spec isa AbstractDict ? get(spec, key, get(spec, String(key), default)) : get(spec, key, default)
        column = field(:column, nothing)
        column isa Union{AbstractString, Symbol} || throw(ArgumentError("Each channel requires a string or symbol column suffix."))
        column = String(column)
        isempty(column) && throw(ArgumentError("Channel column suffix must not be empty."))
        column in seen && throw(ArgumentError("Duplicate channel column $(column)."))
        push!(seen, column)
        label = field(:label, replace(column, '_' => ' '))
        unit = field(:unit, "")
        label isa AbstractString && unit isa AbstractString || throw(ArgumentError("Channel label and unit must be strings."))
        digits = field(:digits, 3)
        digits isa Integer && !(digits isa Bool) && 0 <= digits <= 100 ||
            throw(ArgumentError("Channel digits must be an integer from 0 through 100."))
        log = field(:log, false)
        log isa Bool || log === :auto || (log isa AbstractString && log == "auto") ||
            throw(ArgumentError("Channel log must be true, false, :auto or \"auto\"."))
        log = log isa Bool ? log : "auto"
        # Omit a channel unless every spacecraft has it, keeping spacecraft indices aligned.
        all(i -> "sc$(i)_$(column)" in names(df), 1:S) || continue
        for i in 1:S, value in df[!, "sc$(i)_$(column)"]
            _channel_float32(value, column)
        end
        push!(specs, (column=column, label=String(label), unit=String(unit), digits=Int(digits), log=log))
    end
    return specs
end

# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------

"""
    viewer_import_map(viewer_dir=nothing) -> Dict{String, Any}

Import map whose entries are `data:` URLs of the vendored three.js modules
and the viewer's own modules, so the bundled page resolves every import
offline.
"""
function viewer_import_map(viewer_dir::Union{Nothing, AbstractString}=nothing)::Dict{String, Any}
    return VisualizerServices.viewer_import_map(viewer_dir)
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
Identical bytes with the same format are embedded once; later entries carry
`"url_from" => owning spacecraft id` instead of `"url"`. Transforms and
articulations remain per spacecraft. A `.gltf` file must embed its buffers; external files are not carried along.
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
    bytes_by_path = Dict{String, Vector{UInt8}}()
    owner_by_bytes = Dict{Tuple{String, Vector{UInt8}}, String}()
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
        bytes = get!(() -> read(path), bytes_by_path, path)
        entry = Dict{String, Any}(
            "format" => format,
            "scale" => scale,
            "rotation_deg" => Float64[Float64(r) for r in rot],
            "center" => center,
            "articulations" => articulation_dicts,
            "source" => basename(path),
        )
        # Content and parser format define identity, independent of file name and transforms.
        key = (String(format), bytes)
        owner = get(owner_by_bytes, key, nothing)
        if owner === nothing
            owner_by_bytes[key] = string(sc.id)
            entry["url"] = _data_url(bytes, mime)
        else
            entry["url_from"] = owner
        end
        out[string(sc.id)] = entry
    end
    return out
end

"""
    viewer_payload(scene, df; textures_dir=TEXTURES_DIR, options=Dict(), max_frames=2000, data_budget_mb=150.0, position_precision=:auto, stl=Dict(), stl_scale=1.0) -> Dict{String, Any}

The `window.SPACEAGORA_VIEWER` object: sidecar, encoded frames, the planet's
texture (if registered), STL models (if any) and viewer options.
"""
function viewer_payload(
    scene::VisualizationScene,
    df::DataFrame;
    textures_dir::AbstractString=TEXTURES_DIR,
    include_textures::Bool=true,
    texture_resolution="4k",
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
    channels=(),
    terrain=nothing,
    terrain_max_grid=512,
    mission=nothing,
    position_precision=:auto,
    link_schedule=nothing
)::Dict{String, Any}
    _checked_position_precision(position_precision)
    mission = _validated_mission_history(mission)
    # The links' windows share the data budget with the frames.
    interlinks = link_schedule === nothing ? nothing : interlink_payload(link_schedule, length(scene.spacecraft); ground_stations=scene.ground_stations)
    if interlinks !== nothing
        link_mb = _interlink_payload_bytes(interlinks) / 1.0e6
        link_mb < data_budget_mb || throw(ArgumentError("The interlink windows ($(round(link_mb; digits=3)) MB) " *
            "exceed the viewer's data budget; increase data_budget_mb, or pass interlinks=false to leave the links out."))
        data_budget_mb -= link_mb
    end
    textures = Dict{String, Any}()
    if include_textures
        entry = texture_payload(scene.planet.texture; resolution=texture_resolution, dir=textures_dir)
        entry === nothing || (textures[scene.planet.texture] = entry)
    end
    payload = Dict{String, Any}(
        "viewer_schema" => 1,
        "scene" => scene_dict(scene),
        "frames" => build_viewer_frames(df, scene; max_frames=max_frames, data_budget_mb=data_budget_mb, channels=channels,
            mission=mission, position_precision=position_precision),
        "textures" => textures,
        "models" => model_payloads(scene; models=models, model_scale=model_scale, model_rotation_deg=model_rotation_deg, model_center=model_center, model_articulations=model_articulations, stl=stl, stl_scale=stl_scale),
        "paths" => path_payloads(paths),
        "references" => reference_payloads(references, scene),
        "terrain" => terrain === nothing ? nothing : terrain_payload(terrain; max_grid=terrain_max_grid),
        "options" => Dict{String, Any}(options),
    )
    mission === nothing || (payload["mission"] = deepcopy(mission))
    interlinks === nothing || (payload["interlinks"] = interlinks)
    return payload
end

"""
    interlink_payload(schedule::DataFrame, spacecraft_count) -> Dict{String, Any}

The `interlinks` block of the viewer payload, from a run's link schedule (see
`link_schedule`): `schema` 1; `links`, one entry per link that appears in a
window, in link order, with `link`, `kind`, `directed`, the spacecraft `a` and
`b` (numbered from 1, as `paths[].target`) and their terminals; `active`, the
spans over which each link was selected, as `count`, `index` (0-based positions
in `links`) and Float64 `start_s` and `end_s` on the clock of `frames.t_s`; and
`contact`, the same for the spans over which a link was available, or `nothing`
when the run did not record them. `history_start_s` and `history_end_s` give the recorded interval when known;
times outside it are unknown, not evidence of no links. Valid checkpoint history
preserves the earlier interval after resume.
A link is drawn at playback time `t` when
`start_s <= t < end_s`, and at `t == end_s` when that is its schedule's last time.
Ground station ends retain their negative node number. Supply `ground_stations`
from the scene to include their names and planet-fixed positions in meters.
"""
function interlink_payload(schedule::DataFrame, spacecraft_count::Integer; ground_stations=GroundStationGeometry[])::Dict{String, Any}
    for column in ("window", "link", "kind", "directed", "a_spacecraft", "a_terminal", "b_spacecraft", "b_terminal",
                   "start_s", "end_s")
        column in names(schedule) || throw(ArgumentError("The link schedule has no `$(column)` column."))
    end
    station_nodes = Set(station.node for station in ground_stations)
    length(station_nodes) == length(ground_stations) && all(<(0), station_nodes) ||
        throw(ArgumentError("Ground station nodes must be unique negative integers."))
    all(station -> !isempty(station.name) && all(isfinite, station.position_m), ground_stations) ||
        throw(ArgumentError("Ground stations need names and finite planet-fixed positions."))
    for spacecraft in Iterators.flatten((schedule.a_spacecraft, schedule.b_spacecraft))
        (1 <= spacecraft <= spacecraft_count || spacecraft in station_nodes) || throw(ArgumentError(
            "The link schedule names node $(spacecraft), which is absent from the scene."))
    end
    numbers = sort!(unique(schedule.link))
    position = Dict(number => k - 1 for (k, number) in enumerate(numbers))
    first_row = Dict(number => findfirst(==(number), schedule.link) for number in numbers)
    links = Any[Dict{String, Any}("link" => number, "kind" => String(schedule.kind[row]),
        "directed" => Bool(schedule.directed[row]), "a" => Int(schedule.a_spacecraft[row]),
        "a_terminal" => Int(schedule.a_terminal[row]), "b" => Int(schedule.b_spacecraft[row]),
        "b_terminal" => Int(schedule.b_terminal[row])) for (number, row) in ((n, first_row[n]) for n in numbers)]
    windows(kind) = (rows = findall(==(kind), schedule.window); Dict{String, Any}("count" => length(rows),
        "index" => Int[position[schedule.link[row]] for row in rows],
        "start_s" => _float64_base64(Float64[schedule.start_s[row] for row in rows]),
        "end_s" => _float64_base64(Float64[schedule.end_s[row] for row in rows])))
    recorded = metadata(schedule, "record_available", nothing)
    contact = recorded === true ? windows("contact") : recorded === false ? nothing :
        any(==("contact"), schedule.window) ? windows("contact") : nothing
    finite_time(key) = begin
        value = metadata(schedule, key, NaN)
        value isa Real && isfinite(value) ? Float64(value) : nothing
    end
    block = Dict{String, Any}("schema" => 1, "links" => links, "active" => windows("active"), "contact" => contact,
        "history_start_s" => finite_time("history_start_s"), "history_end_s" => finite_time("history_end_s"))
    isempty(ground_stations) || (block["stations"] = [_scene_dict(station) for station in ground_stations])
    return block
end

# The bytes the block adds to the payload, as JSON text.
_interlink_payload_bytes(block) = sizeof(JSON.json(block))

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
    render_viewer_html(payload; viewer_dir=nothing, title="SpaceAGORA viewer") -> String

Render through the loaded SpaceAGORAVisualizer companion. `viewer_dir` optionally selects a custom renderer asset directory.
"""
function render_viewer_html(payload::AbstractDict; viewer_dir::Union{Nothing, AbstractString}=nothing, title::AbstractString="SpaceAGORA viewer")::String
    return VisualizerServices.render_viewer_html(payload; viewer_dir=viewer_dir, title=title)
end

@inline function _viewer_options(; trail_s, trail_orbits, frame, speed, title, ground_tracks=false)
    frame in (:inertial, :planet_fixed) || throw(ArgumentError("frame must be :inertial or :planet_fixed, got $(frame)."))
    options = Dict{String, Any}("frame" => String(frame))
    trail_s === nothing || (options["trail_s"] = Float64(trail_s))
    trail_orbits === nothing || (options["trail_orbits"] = Float64(trail_orbits))
    speed === nothing || (options["speed"] = Float64(speed))
    title === nothing || (options["title"] = String(title))
    options["ground_tracks"] = Bool(ground_tracks)
    return options
end

# The result manifest owns the association, so a normal result overwrite drops
# it even when the new trajectory happens to contain identical bytes. Hashes
# detect separately replaced history, scene or Arrow files. Relative names keep
# the association valid when the complete saved folder is moved.
function _mission_file_sha256(path::AbstractString)
    return open(path, "r") do io
        bytes2hex(SHA.sha256(io))
    end
end

# The manifest's entry for the results feather a scene names.
_manifest_feather_key(scene::VisualizationScene) = scene.results_layout === :long ? "spacecraft_feather" : "feather"

# The results a saved scene's viewer reads, as the wide table. A wide run's
# feather is read as it is. A long run's spacecraft table is opened memory-mapped
# and its columns become the wide `sc{i}_<column>` columns as views of the rows
# of spacecraft `i` (and the run table's columns are added), so no row is copied:
# the viewer then reads only the rows it keeps.
function _saved_results_table(scene::VisualizationScene, feather_path::AbstractString)::DataFrame
    scene.results_layout === :long || return DataFrame(Arrow.Table(feather_path))
    run_path = replace(String(feather_path), r"_spacecraft\.feather$" => "_run.feather")
    run = run_path != feather_path && isfile(run_path) ? DataFrame(Arrow.Table(run_path); copycols=false) : nothing
    arrow_table = Arrow.Table(feather_path)
    long = DataFrame(arrow_table; copycols=false)
    arrow_metadata = Arrow.getmetadata(arrow_table)
    key = "spaceagora.long_column_spacecraft.v1"
    if arrow_metadata !== nothing && haskey(arrow_metadata, key)
        metadata!(long, key, arrow_metadata[key]; style=:note)
    end
    return _wide_view(long, run, length(scene.spacecraft))
end

const _LONG_IDENTITY_COLUMNS = ("spacecraft", "spacecraft_id", "plane", "slot", "time")

"""
    _wide_view(spacecraft_table, run_table, num_sats) -> DataFrame

The wide table of a long spacecraft table (one row per spacecraft and time), as
views: `time`, then `sc{i}_<column>` for each spacecraft `i` and column, then the
run table's columns. The rows may come in any order but must hold every
spacecraft at every time once.
"""
function _wide_view(long::DataFrame, run::Union{Nothing, DataFrame}, num_sats::Integer)::DataFrame
    raw_presence = metadata(long, "spaceagora.long_column_spacecraft.v1", nothing)
    presence = if raw_presence === nothing
        nothing
    else
        parsed = try
            JSON.parse(raw_presence)
        catch
            throw(ArgumentError("Invalid long-results column-presence metadata."))
        end
        parsed isa AbstractDict || throw(ArgumentError("Invalid long-results column-presence metadata."))
        for name in names(long)
            name in _LONG_IDENTITY_COLUMNS && continue
            members = get(parsed, name, nothing)
            members isa AbstractVector && all(i -> i isa Integer && !(i isa Bool) && 1 <= i <= num_sats, members) &&
                length(unique(members)) == length(members) ||
                throw(ArgumentError("Invalid long-results column-presence metadata for column $(name)."))
        end
        parsed
    end
    spacecraft, times = long[!, "spacecraft"], long[!, "time"]
    n = length(spacecraft)
    n % num_sats == 0 || throw(ArgumentError("The long results table does not hold every spacecraft at every time."))
    in_order = all(r -> spacecraft[r] == mod1(r, num_sats) && times[r] == times[r - mod1(r, num_sats) + 1], 1:n)
    order = in_order ? collect(1:n) : sortperm(1:n; by=r -> (times[r], spacecraft[r]))
    all(r -> spacecraft[order[r]] == mod1(r, num_sats) && times[order[r]] == times[order[r - mod1(r, num_sats) + 1]], 1:n) ||
        throw(ArgumentError("The long results table does not hold every spacecraft at every time."))
    rows = [view(order, i:num_sats:n) for i in 1:num_sats]
    wide = DataFrame("time" => view(times, rows[1]); copycols=false)
    for name in names(long)
        name in _LONG_IDENTITY_COLUMNS && continue
        column = long[!, name]
        for i in 1:num_sats
            spacecraft_column = view(column, rows[i])
            if presence === nothing
                # Compatibility with older draft bundles: their all-missing columns
                # cannot distinguish absent fields from present measurement gaps.
                Missing <: eltype(column) && all(ismissing, spacecraft_column) && continue
            else
                i in presence[name] || continue
            end
            wide[!, "sc$(i)_$(name)"] = spacecraft_column
        end
    end
    if run !== nothing
        nrow(run) == nrow(wide) && run[!, "time"] == wide[!, "time"] ||
            throw(ArgumentError("The run table's times differ from the spacecraft table's."))
        for name in names(run)
            name == "time" || (wide[!, name] = run[!, name])
        end
    end
    return wide
end

function _bind_mission_history!(prefix::AbstractString)
    manifest_path = String(prefix) * ".manifest.toml"
    scene_path = String(prefix) * "_scene.json"
    history_path = String(prefix) * "_mission.json"
    scene = read_visualization_scene(scene_path)
    feather_path = joinpath(dirname(scene_path), scene.results_feather)
    _validated_mission_history(JSON.parsefile(history_path))
    manifest = TOML.parsefile(manifest_path)
    feather_hash = _mission_file_sha256(feather_path)
    get(get(get(manifest, "files", Dict()), _manifest_feather_key(scene), Dict()), "sha256", nothing) == feather_hash ||
        throw(ArgumentError("Cannot associate mission history with a different saved results bundle."))
    manifest["mission_history"] = Dict(
        "schema" => 1, "path" => basename(history_path),
        "sha256" => _mission_file_sha256(history_path),
        "feather_sha256" => feather_hash,
        "scene_sha256" => _mission_file_sha256(scene_path))
    IOSerialization._atomic_write_file(manifest_path, tmp -> open(tmp, "w") do io
        TOML.print(io, manifest)
    end)
    return nothing
end

function _saved_mission_history(prefix::AbstractString, scene_path::String, feather_path::String)
    manifest_path = String(prefix) * ".manifest.toml"
    isfile(manifest_path) || return nothing
    manifest = try
        TOML.parsefile(manifest_path)
    catch err
        err isa TOML.ParserError || rethrow()
        return nothing
    end
    association = get(manifest, "mission_history", nothing)
    association === nothing && return nothing
    history_path = String(prefix) * "_mission.json"
    association isa AbstractDict && get(association, "schema", nothing) === 1 &&
        get(association, "path", nothing) == basename(history_path) ||
        throw(ArgumentError("Unsupported saved mission history association in $manifest_path."))
    for (path, key) in ((history_path, "sha256"), (feather_path, "feather_sha256"), (scene_path, "scene_sha256"))
        isfile(path) && get(association, key, nothing) == _mission_file_sha256(path) ||
            throw(ArgumentError("Saved mission history no longer matches $path; export the matching saved bundle or supply mission explicitly."))
    end
    return JSON.parsefile(history_path)
end

function _saved_viewer_payload(
    prefix::AbstractString;
    max_frames::Integer=DEFAULT_MAX_FRAMES,
    data_budget_mb::Real=DEFAULT_DATA_BUDGET_MB,
    trail_s::Union{Nothing, Real}=nothing,
    trail_orbits::Union{Nothing, Real}=nothing,
    frame::Symbol=:inertial,
    speed::Union{Nothing, Real}=nothing,
    title::Union{Nothing, AbstractString}=nothing,
    ground_tracks::Bool=false,
    textures::Bool=true,
    texture_resolution="4k",
    models::AbstractDict=Dict{Int, String}(),
    model_scale=nothing,
    model_rotation_deg::AbstractDict=Dict{Int, Any}(),
    model_center=true,
    model_articulations::AbstractDict=Dict{Int, Any}(),
    stl::AbstractDict=Dict{Int, String}(),
    stl_scale::Real=1.0,
    paths=(),
    references=(),
    channels=(),
    terrain=nothing,
    terrain_max_grid=512,
    textures_dir::AbstractString=TEXTURES_DIR,
    mission=nothing,
    position_precision=:auto,
    interlinks::Bool=true
)
    _checked_position_precision(position_precision)
    prefix = String(prefix)
    scene_path = prefix * "_scene.json"
    isfile(scene_path) || throw(ArgumentError("No visualization sidecar at $(scene_path); run with simulation_settings.save_visualization_scene = true."))
    scene = read_visualization_scene(scene_path)
    feather_path = joinpath(dirname(scene_path), scene.results_feather)
    isfile(feather_path) || throw(ArgumentError("Results file $(feather_path) referenced by the sidecar does not exist."))
    df = _saved_results_table(scene, feather_path)
    mission === nothing && (mission = _saved_mission_history(prefix, scene_path, feather_path))
    link_schedule = interlinks ? _saved_link_schedule(prefix, scene) : nothing
    payload = viewer_payload(
        scene, df; link_schedule=link_schedule,
        textures_dir=textures_dir, include_textures=textures, texture_resolution=texture_resolution,
        options=_viewer_options(; trail_s=trail_s, trail_orbits=trail_orbits, frame=frame, speed=speed, title=title, ground_tracks=ground_tracks),
        max_frames=max_frames, data_budget_mb=data_budget_mb, position_precision=position_precision,
        models=models, model_scale=model_scale, model_rotation_deg=model_rotation_deg, model_center=model_center, model_articulations=model_articulations, stl=stl, stl_scale=stl_scale,
        paths=paths, references=references, channels=channels, terrain=terrain, terrain_max_grid=terrain_max_grid, mission=mission
    )
    page_title = title === nothing ? "SpaceAGORA · $(scene.planet.name) · $(basename(prefix))" : String(title)
    return payload, page_title
end

# The link schedule the scene names, checked against the hash the results manifest
# recorded for it, so a schedule written by another run is never drawn; `nothing`
# when the scene names none.
function _saved_link_schedule(prefix::AbstractString, scene::VisualizationScene)
    scene.link_schedule === nothing && return nothing
    path = joinpath(dirname(String(prefix)), scene.link_schedule)
    manifest_path = String(prefix) * ".manifest.toml"
    manifest = isfile(manifest_path) ? TOML.parsefile(manifest_path) : Dict{String, Any}()
    files = get(manifest, "files", Dict())
    trajectory_path = joinpath(dirname(String(prefix)), scene.results_feather)
    for (file, entry) in ((path, "link_schedule_feather"), (trajectory_path, _manifest_feather_key(scene)))
        recorded = get(get(files, entry, Dict()), "sha256", nothing)
        isfile(file) && recorded !== nothing && recorded == IOSerialization._sha256_hex(file) ||
            throw(ArgumentError("The interlink input $(file) is missing or no longer matches its results manifest; export the " *
                "matching saved bundle, or pass interlinks=false to leave the links out."))
    end
    schedule = DataFrame(Arrow.Table(path))
    entry = get(manifest, "interlinks", Dict{String, Any}())
    for key in ("history_start_s", "history_end_s", "record_available")
        haskey(entry, key) && metadata!(schedule, key, entry[key]; style=:note)
    end
    return schedule
end

function _check_viewer_output(path::AbstractString, protected)
    isdir(path) && throw(ArgumentError("Viewer output must be a file, not a directory: $path"))
    for input in protected
        (abspath(path) == abspath(input) || (isfile(path) && isfile(input) && samefile(path, input))) &&
            throw(ArgumentError("Viewer output cannot replace saved input data: $input"))
    end
    return String(path)
end

function _viewer_output_path(prefix::AbstractString, payload::AbstractDict, out, suffix::String)
    path = out === nothing ? String(prefix) * suffix : String(out)
    schedule = get(payload["scene"]["results"], "link_schedule", nothing)
    protected = (String(prefix) * "_scene.json", String(prefix) * "_mission.json", String(prefix) * ".manifest.toml",
        joinpath(dirname(String(prefix)), String(payload["scene"]["results"]["feather"])),
        (schedule === nothing ? () : (joinpath(dirname(String(prefix)), String(schedule["file"])),))...)
    return _check_viewer_output(path, protected)
end

"""
    export_viewer_payload(prefix; out=nothing, kwargs...) -> String
    export_viewer_payload(configuration; kwargs...) -> String

Convert a saved scene and results bundle to version-1 viewer JSON. The default
file is `<prefix>_viewer.json`. Frame conventions, precision, decimation and
options match `export_visualization`; `textures=false` omits local textures.
The result manifest links `<prefix>_mission.json` to its own saved run. Only a
verified association is carried into the optional `mission` block; an unbound
sidecar is never attached automatically. Move the complete saved folder to
preserve that association, or supply a deliberate `mission` keyword value.
Exact saved boundary rows are retained within the frame budget; missing boundary
states are never reconstructed. A run with interlinks adds its link schedule's
windows as the `interlinks` block (see `interlink_payload`), checked against the
results manifest; `interlinks=false` leaves them out. No renderer package, browser, download or simulation is needed. Load
SpaceAGORAVisualizer to turn this portable payload into an HTML page.
"""
function export_viewer_payload(prefix::AbstractString; out::Union{Nothing, AbstractString}=nothing, kwargs...)::String
    payload, _ = _saved_viewer_payload(prefix; kwargs...)
    path = _viewer_output_path(prefix, payload, out, "_viewer.json")
    return IOSerialization._atomic_write_file(path, tmp -> write(tmp, JSON.json(payload)))
end

function export_viewer_payload(args::SimulationConfiguration; kwargs...)::String
    return export_viewer_payload(IOConfig._results_bundle_prefix(args); kwargs...)
end

"""
    export_visualization(prefix; out=prefix * "_viewer.html", kwargs...) -> String
    export_visualization(args::SimulationConfiguration; kwargs...) -> String

Pass a results prefix or the configuration that owns the results bundle.
Build the self-contained viewer page for a run written with
`save_visualization_scene = true`. `prefix` is the results bundle prefix
(`<results_directory>/simulation_results`); the sidecar `<prefix>_scene.json`
and the Arrow file `<prefix>.feather` must exist.

Keywords: `max_frames` (default 2000) and `data_budget_mb` (default 150)
bound the embedded trajectory. Positions are Float64 for 64 spacecraft or
fewer, and above that when Float64 keeps the frames Float32 would; otherwise
Float32. `position_precision=:f64` keeps Float64 at the cost of frames (with
mission segments, it can need a larger budget for their boundary rows), and
`:f32` forces Float32 (`:auto` by default). `trail_orbits` (default 3, estimated from the
trajectory's periapsis passages) or `trail_s` set the trail window; `frame`
is `:inertial` (default) or `:planet_fixed`; `speed` is the initial playback
rate in simulated seconds per wall second; `title` names the page;
`textures=false` skips the surface texture and `texture_resolution` picks a
tier (`"4k"` by default; `"8k"` opts into more detail for Earth, Mars and Moon;
`:best` takes the largest registered tier); `ground_tracks=true` starts the page with the
sub-satellite tracks drawn on the body's surface (the toolbar toggles them
either way); `models` maps spacecraft ids to STL, OBJ, glTF
or GLB files drawn instead of the link boxes, at `model_scale` meters per
model unit (a number or a per-id `Dict`) and rotated by `model_rotation_deg`
(per-id XYZ Euler angles), posed by `model_articulations` (per-id list of
parts rotated about an axis, see `articulate_triangles`) and centered on the
spacecraft unless `model_center=false`; `stl`/`stl_scale` are the older spelling for STL only; `references` draws
reference trajectories as translucent ghosts of a spacecraft (see
`reference_payloads`); `paths` overlays
reference polylines (see `path_payloads`), e.g. a planned RPO path in the
target's RTN frame; `channels` is a tuple or vector of NamedTuples/Dicts with a
`column` suffix (`sc<i>_<column>`), optional string `label`/`unit`, integer
`digits` in 0:100 (default 3), and `log` true, false or `:auto`/`"auto"`.
A channel missing from any spacecraft is omitted. Samples are Float32;
`missing`/NaN are gaps and infinities/overflow are rejected. Channels use the
same decimated rows and budget as positions. `terrain` names a local site.json;
`terrain_max_grid` (default 512) limits its exported DEM dimensions (see
`terrain_payload`). A run with interlinks draws its active links from its link
schedule, which shares `data_budget_mb` with the frames; `interlinks=false`
leaves them out. `viewer_dir` overrides the optional visualizer package assets; `textures_dir`
overrides the local texture directory. Load `SpaceAGORAVisualizer` before
requesting an HTML page, or use `export_viewer_payload` to save portable JSON.
"""
function export_visualization(prefix::AbstractString; out::Union{Nothing, AbstractString}=nothing,
    viewer_dir::Union{Nothing, AbstractString}=nothing, kwargs...)::String
    VisualizerServices.require_available()
    payload, page_title = _saved_viewer_payload(prefix; kwargs...)
    path = _viewer_output_path(prefix, payload, out, "_viewer.html")
    payload_path = String(prefix) * "_viewer.json"
    (abspath(path) == abspath(payload_path) ||
        (isfile(path) && isfile(payload_path) && samefile(path, payload_path))) &&
        throw(ArgumentError("HTML output cannot replace the viewer JSON input: $payload_path"))
    html = render_viewer_html(payload; viewer_dir=viewer_dir, title=page_title)
    return IOSerialization._atomic_write_file(path, tmp -> write(tmp, html))
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

Write `window.SPACEAGORA_VIEWER = {...};` for `viewer/dev.html` in the
SpaceAGORAVisualizer.jl companion repository, so its viewer modules can be
edited and reloaded without rebundling.
"""
function write_viewer_dev_payload(prefix::AbstractString, out_js::AbstractString; textures::Bool=true, textures_dir::AbstractString=TEXTURES_DIR, mission=nothing, interlinks::Bool=true, kwargs...)::String
    _checked_position_precision(get(kwargs, :position_precision, :auto))
    scene_path = String(prefix) * "_scene.json"
    scene = read_visualization_scene(scene_path)
    feather_path = joinpath(dirname(scene_path), scene.results_feather)
    df = _saved_results_table(scene, feather_path)
    mission === nothing && (mission = _saved_mission_history(prefix, scene_path, feather_path))
    payload = viewer_payload(scene, df; textures_dir=textures_dir, include_textures=textures, mission=mission,
        link_schedule=interlinks ? _saved_link_schedule(prefix, scene) : nothing, kwargs...)
    path = _viewer_output_path(prefix, payload, out_js, "_viewer.js")
    _check_viewer_output(path, (String(prefix) * "_viewer.json",))
    return IOSerialization._atomic_write_file(path, tmp -> write(tmp, "window.SPACEAGORA_VIEWER = " * _script_safe_json(payload) * ";\n"))
end
