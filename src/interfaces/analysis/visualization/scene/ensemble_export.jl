# Ensembles: many runs of the same configuration (Monte Carlo samples, or the
# per-member bundles a constellation ensemble writes) merged into one viewer
# page. Every sample is a subdirectory of a campaign directory holding its own
# results bundle and scene sidecar; they are resampled onto one time axis and
# become S pseudo-spacecraft in the frame payload, coloured by a per-sample
# scalar in the viewer.

using ..SimConfig: _with_fields

const ENSEMBLE_MANIFEST = "ensemble_manifest.json"
const ENSEMBLE_SAMPLE_DIR_RE = r"^(sample_\d+|sat_\d+_id_\d+)$"

"""
    EnsembleSample

One member of an ensemble: `index` (1-based), the `seed` it was run with (kept
as a string for the manifest), whether the run finished, a scalar the viewer
colors it by (`NaN` when none), a display label and its results directory.
"""
struct EnsembleSample
    index::Int
    seed::String
    success::Bool
    scalar::Float64
    label::String
    directory::String
end

"""
    sample_results_directory(campaign_dir, index) -> String

`<campaign_dir>/sample_0001`-style directory for one ensemble member.
"""
@inline sample_results_directory(campaign_dir::AbstractString, index::Integer)::String =
    joinpath(String(campaign_dir), "sample_" * lpad(string(Int(index)), 4, '0'))

"""
    with_results_directory(args, dir; visualization=true) -> SimulationConfiguration

Copy of `args` writing its bundle to `dir`, with the scene sidecar flag set
unless `visualization=false`.
"""
function with_results_directory(args::SimulationConfiguration, dir::AbstractString; visualization::Bool=true)::SimulationConfiguration
    settings = _with_fields(args.simulation_settings;
        results=true, results_directory=String(dir), save_visualization_scene=visualization)
    return _with_configuration(args;
        simulation_settings=settings,
    )
end

"""
    default_sample_scalar(df) -> Float64

The last spherical `sc1_periapsis_altitude` value in kilometers, or `NaN` when the
column is absent; the viewer's default colouring for an ensemble.
"""
function default_sample_scalar(df::DataFrame)::Float64
    "sc1_periapsis_altitude" in names(df) || return NaN
    nrow(df) == 0 && return NaN
    return Float64(df[end, "sc1_periapsis_altitude"]) / 1000.0
end

@inline _sample_dict(s::EnsembleSample) = Dict{String, Any}(
    "index" => s.index, "seed" => s.seed, "success" => s.success,
    "scalar" => isfinite(s.scalar) ? s.scalar : nothing, "label" => s.label, "directory" => s.directory
)

function _sample_from(d, campaign_dir::AbstractString)::EnsembleSample
    scalar = get(d, "scalar", nothing)
    dir = String(get(d, "directory", ""))
    isabspath(dir) || (dir = joinpath(campaign_dir, dir))
    return EnsembleSample(Int(d["index"]), String(get(d, "seed", "")), Bool(get(d, "success", true)),
                          scalar === nothing ? NaN : Float64(scalar), String(get(d, "label", "sample $(d["index"])")), dir)
end

"""
    write_ensemble_manifest(campaign_dir, samples; scalar_name="", nominal=nothing) -> String

Write `ensemble_manifest.json` listing the samples and the name of the scalar
they carry. Directories are stored relative to `campaign_dir` when possible.
`nominal` is the `index` of the unperturbed sample, which the viewer draws
distinctly and centers the 3-sigma tube on.
"""
function write_ensemble_manifest(campaign_dir::AbstractString, samples::AbstractVector{EnsembleSample}; scalar_name::AbstractString="", nominal::Union{Nothing, Integer}=nothing)::String
    nominal === nothing || any(s -> s.index == nominal, samples) || throw(ArgumentError("nominal sample index $(nominal) is not among the samples."))
    entries = Any[]
    for s in samples
        d = _sample_dict(s)
        d["directory"] = startswith(s.directory, String(campaign_dir)) ? relpath(s.directory, String(campaign_dir)) : s.directory
        push!(entries, d)
    end
    payload = Dict{String, Any}("schema" => 1, "scalar_name" => String(scalar_name), "samples" => entries, "nominal" => nominal)
    path = joinpath(String(campaign_dir), ENSEMBLE_MANIFEST)
    return IOSerialization._atomic_write_file(path, tmp -> open(io -> JSON.print(io, payload), tmp, "w"))
end

"""
    read_ensemble_manifest(campaign_dir) -> (samples::Vector{EnsembleSample}, scalar_name::String)
"""
function read_ensemble_manifest(campaign_dir::AbstractString)
    path = joinpath(String(campaign_dir), ENSEMBLE_MANIFEST)
    d = JSON.parsefile(path)
    samples = EnsembleSample[_sample_from(x, String(campaign_dir)) for x in d["samples"]]
    return samples, String(get(d, "scalar_name", ""))
end

"""
    read_ensemble_nominal(campaign_dir) -> Union{Nothing, Int}

The `index` of the nominal sample recorded in the manifest, or `nothing`.
"""
function read_ensemble_nominal(campaign_dir::AbstractString)::Union{Nothing, Int}
    path = joinpath(String(campaign_dir), ENSEMBLE_MANIFEST)
    isfile(path) || return nothing
    v = get(JSON.parsefile(path), "nominal", nothing)
    return v === nothing ? nothing : Int(v)
end

"""
    discover_ensemble_samples(campaign_dir) -> (samples, scalar_name)

The manifest when present; otherwise every `sample_NNNN` or
`sat_<i>_id_<id>` subdirectory holding a scene sidecar, in name order. A
constellation ensemble's folder with the `constellation_summary.csv` that
`run_constellation_ensemble(args; merge_results=true)` writes labels each member
by its id, plane and slot (`"id 7, plane 2, slot 3"`) and gives it its plane as
the scalar, `"plane"`, so the page colours members by plane. Otherwise samples
are labelled by folder name, with no scalar.
"""
function discover_ensemble_samples(campaign_dir::AbstractString)
    isfile(joinpath(String(campaign_dir), ENSEMBLE_MANIFEST)) && return read_ensemble_manifest(campaign_dir)
    identities = _constellation_member_identities(String(campaign_dir))
    samples = EnsembleSample[]
    for name in sort(readdir(String(campaign_dir)))
        occursin(ENSEMBLE_SAMPLE_DIR_RE, name) || continue
        label, plane, destination = get(identities, name, (name, NaN, missing))
        dir = IOSerialization._member_results_directory(joinpath(String(campaign_dir), name), destination)
        isfile(joinpath(dir, "simulation_results_scene.json")) || continue
        push!(samples, EnsembleSample(length(samples) + 1, "", true, plane, label, dir))
    end
    return samples, any(sample -> !isnan(sample.scalar), samples) ? "plane" : ""
end

# The members of a constellation ensemble by folder name, from its
# `constellation_summary.csv`: a label from the member's id, plane and slot, and
# its plane (`NaN` without one). Empty without the summary.
function _constellation_member_identities(campaign_dir::String)
    path = joinpath(campaign_dir, "constellation_summary.csv")
    isfile(path) || return Dict{String, Tuple{String, Float64, Union{Missing, String}}}()
    summary = CSV.read(path, DataFrame)
    identities = Dict{String, Tuple{String, Float64, Union{Missing, String}}}()
    for row in eachrow(summary)
        located = !ismissing(row.plane) && !ismissing(row.slot)
        label = located ? "id $(row.spacecraft_id), plane $(row.plane), slot $(row.slot)" : "id $(row.spacecraft_id)"
        destination = :results_subdirectory in propertynames(row) ? row.results_subdirectory : missing
        identities[String(row.member_tag)] = (label, located ? Float64(row.plane) : NaN, destination)
    end
    return identities
end

# ---------------------------------------------------------------------------
# Resampling onto a common time axis
# ---------------------------------------------------------------------------

# Linear interpolation of `y` (sampled at strictly increasing `t`) at `tq`; NaN outside [t[1], t[end]].
function _resample_linear!(out::AbstractVector{Float64}, t::AbstractVector{<:Real}, y::AbstractVector{<:Real}, tq::AbstractVector{<:Real})
    n = length(t)
    @inbounds for (k, q) in enumerate(tq)
        if n == 0 || q < t[1] || q > t[n]
            out[k] = NaN
        elseif n == 1
            out[k] = Float64(y[1])
        else
            hi = searchsortedfirst(t, q)
            hi = clamp(hi, 2, n)
            lo = hi - 1
            f = (Float64(q) - Float64(t[lo])) / max(Float64(t[hi]) - Float64(t[lo]), eps())
            out[k] = Float64(y[lo]) + clamp(f, 0.0, 1.0) * (Float64(y[hi]) - Float64(y[lo]))
        end
    end
    return out
end

@inline function _slerp(a::SVector{4, Float64}, b::SVector{4, Float64}, f::Float64)::SVector{4, Float64}
    d = dot(a, b)
    if d < 0.0
        b = -b
        d = -d
    end
    if d > 0.9995
        q = a + f * (b - a)
        return q / norm(q)
    end
    θ = acos(clamp(d, -1.0, 1.0))
    return (sin((1 - f) * θ) * a + sin(f * θ) * b) / sin(θ)
end

function _resample_quaternions!(out::AbstractMatrix{Float64}, t::AbstractVector{<:Real}, q::AbstractMatrix{<:Real}, tq::AbstractVector{<:Real})
    # `q` is 4 x n, `out` is 4 x length(tq)
    n = length(t)
    @inbounds for (k, s) in enumerate(tq)
        if n == 0 || s < t[1] || s > t[n]
            out[:, k] .= NaN
        elseif n == 1
            out[:, k] .= q[:, 1]
        else
            hi = clamp(searchsortedfirst(t, s), 2, n)
            lo = hi - 1
            f = clamp((Float64(s) - Float64(t[lo])) / max(Float64(t[hi]) - Float64(t[lo]), eps()), 0.0, 1.0)
            a = SVector{4, Float64}(q[1, lo], q[2, lo], q[3, lo], q[4, lo])
            b = SVector{4, Float64}(q[1, hi], q[2, hi], q[3, hi], q[4, hi])
            out[:, k] .= _slerp(a, b, f)
        end
    end
    return out
end

"""
    ensemble_time_axis(ends_s, cadences_s; max_frames=2000) -> Vector{Float64}

Common playback axis from 0 to the longest sample, at the finest sample
cadence or coarser so it holds at most `max_frames` rows.
"""
function ensemble_time_axis(ends_s::AbstractVector{<:Real}, cadences_s::AbstractVector{<:Real}; max_frames::Integer=DEFAULT_MAX_FRAMES)::Vector{Float64}
    isempty(ends_s) && return Float64[0.0]
    t_end = maximum(Float64, ends_s)
    t_end > 0.0 || return Float64[0.0]
    cadence = max(minimum(Float64, cadences_s), t_end / (max(Int(max_frames), 2) - 1))
    return rotation_sample_times(t_end, cadence; max_samples=max(Int(max_frames), 2))
end

@inline function _table_cadence(t::AbstractVector{<:Real})::Float64
    length(t) < 2 && return 1.0
    d = Float64(t[end] - t[1]) / (length(t) - 1)
    return d > 0 ? d : 1.0
end

# A sample's results table: a `DataFrame`, or a function that reads it, so that
# the export holds one sample's table at a time.
_ensemble_table(table::DataFrame) = table
_ensemble_table(read_table) = read_table()::DataFrame

# What the frames need of a sample's table before its rows: its column names,
# last saved time, cadence and row count.
function _ensemble_table_outline(table)
    df = _ensemble_table(table)
    return (names=Set{String}(names(df)), last_s=Float64(df.time[end]), cadence=_table_cadence(df.time), rows=nrow(df))
end

"""
    build_ensemble_frames(samples, scenes, tables; max_frames=2000, data_budget_mb=150.0, position_precision=:auto) -> (frames::Dict, scene::VisualizationScene, per_sample::Int)

Merge the samples' results tables onto one time axis; each entry of `tables`
is a `DataFrame` or a function that reads one. Every sample must have
the same number of spacecraft; sample `i`, spacecraft `k` becomes
pseudo-spacecraft `(i-1)*per_sample + k`. Rows outside a sample's own span are
`NaN`, which the viewer treats as "not present". Positions follow the rule of
`build_viewer_frames` for the pseudo-spacecraft count.
"""
function build_ensemble_frames(
    samples::AbstractVector{EnsembleSample},
    scenes::AbstractVector{VisualizationScene},
    tables::AbstractVector;
    max_frames::Integer=DEFAULT_MAX_FRAMES,
    data_budget_mb::Real=DEFAULT_DATA_BUDGET_MB,
    position_precision=:auto
)
    _checked_position_precision(position_precision)
    n = length(samples)
    n >= 1 || throw(ArgumentError("An ensemble needs at least one sample."))
    length(scenes) == n && length(tables) == n || throw(ArgumentError("samples, scenes and tables must have the same length."))
    per_sample = length(scenes[1].spacecraft)
    all(sc -> length(sc.spacecraft) == per_sample, scenes) || throw(ArgumentError("Every ensemble sample must have the same number of spacecraft."))
    planet = scenes[1].planet.name
    all(sc -> sc.planet.name == planet, scenes) || throw(ArgumentError("Every ensemble sample must orbit the same body."))
    S = n * per_sample
    # A first pass reads each table for its outline; the second, below, reads
    # each again for its rows.
    outlines = [_ensemble_table_outline(table) for table in tables]
    has(outline, columns) = all(column -> column in outline.names, columns)
    has_vel = all(o -> all(k -> has(o, ("sc$(k)_vel_1", "sc$(k)_vel_2", "sc$(k)_vel_3")), 1:per_sample), outlines)
    has_q = all(o -> all(k -> has(o, ("sc$(k)_q_1", "sc$(k)_q_2", "sc$(k)_q_3", "sc$(k)_q_4")), 1:per_sample), outlines)
    has_mass = all(o -> all(k -> "sc$(k)_mass" in o.names, 1:per_sample), outlines)
    has_density = all(o -> all(k -> "sc$(k)_density" in o.names, 1:per_sample), outlines)
    has_heat = all(o -> all(k -> "sc$(k)_heat_rate" in o.names, 1:per_sample), outlines)
    stride_lp = scenes[1].link_pose_stride
    counts = Int[max(0, length(sc.links) - 1) for sc in scenes[1].spacecraft]
    has_lp = any(>(0), counts) && all(outlines) do o
        all(1:per_sample) do k
            counts[k] == 0 || has(o, ["sc$(k)_link_pose_$(j)" for j in 1:(stride_lp * counts[k])])
        end
    end
    lp_counts = repeat(counts, n)
    lp_offsets = Int[]
    lp_total = 0
    for c in lp_counts
        push!(lp_offsets, lp_total)
        lp_total += has_lp ? stride_lp * c : 0
    end
    # The time axis within the budget at `position_bytes` per position.
    other_sat_bytes = (has_vel ? 12 : 0) + (has_q ? 16 : 0) + (has_mass ? 4 : 0)
    other_frame_bytes = 4 * lp_total + 8
    ends = [o.last_s for o in outlines]
    cadences = [o.cadence for o in outlines]
    budget_frames(position_bytes) =
        floor(Int, Float64(data_budget_mb) * 1.0e6 / (S * (position_bytes + other_sat_bytes) + other_frame_bytes))
    time_axis(position_bytes) =
        ensemble_time_axis(ends, cadences; max_frames=max(2, min(Int(max_frames), budget_frames(position_bytes))))
    # Below two frames the floor keeps two at either precision, over the budget:
    # that costs, so Float32.
    pos_f64 = _float64_positions(position_precision, S, () -> budget_frames(24) >= 2 && time_axis(24) == time_axis(12))
    times = time_axis(pos_f64 ? 24 : 12)
    N = length(times)

    pos = fill(NaN, N * S * 3)
    vel = has_vel ? fill(NaN, N * S * 3) : Float64[]
    q = has_q ? fill(NaN, N * S * 4) : Float64[]
    mass = has_mass ? fill(NaN, N * S) : Float64[]
    density = has_density ? fill(NaN, N * S) : Float64[]
    heat = has_heat ? fill(NaN, N * S) : Float64[]
    lp = has_lp ? fill(NaN, N * lp_total) : Float64[]
    col = Vector{Float64}(undef, N)
    for i in 1:n
        df = _ensemble_table(tables[i])
        t = Float64.(df.time)
        for k in 1:per_sample
            s = (i - 1) * per_sample + k
            for c in 1:3
                _resample_linear!(col, t, df[!, "sc$(k)_pos_$(c)"], times)
                @inbounds for f in 1:N
                    pos[((f - 1) * S + (s - 1)) * 3 + c] = col[f] / 1000.0
                end
                if has_vel
                    _resample_linear!(col, t, df[!, "sc$(k)_vel_$(c)"], times)
                    @inbounds for f in 1:N
                        vel[((f - 1) * S + (s - 1)) * 3 + c] = col[f] / 1000.0
                    end
                end
            end
            if has_mass
                _resample_linear!(col, t, df[!, "sc$(k)_mass"], times)
                @inbounds for f in 1:N
                    mass[(f - 1) * S + s] = col[f]
                end
            end
            if has_density
                _resample_linear!(col, t, df[!, "sc$(k)_density"], times)
                @inbounds for f in 1:N
                    density[(f - 1) * S + s] = col[f]
                end
            end
            if has_heat
                _resample_linear!(col, t, df[!, "sc$(k)_heat_rate"], times)
                @inbounds for f in 1:N
                    heat[(f - 1) * S + s] = col[f] * HEAT_RATE_W_CM2_TO_W_M2
                end
            end
            if has_q
                qsrc = Matrix{Float64}(undef, 4, nrow(df))
                for c in 1:4
                    qsrc[c, :] .= Float64.(df[!, "sc$(k)_q_$(c)"])
                end
                qdst = Matrix{Float64}(undef, 4, N)
                _resample_quaternions!(qdst, t, qsrc, times)
                @inbounds for f in 1:N, c in 1:4
                    q[((f - 1) * S + (s - 1)) * 4 + c] = qdst[c, f]
                end
            end
            if has_lp && counts[k] > 0
                for link in 1:counts[k]
                    base = stride_lp * (link - 1)
                    for c in 1:3
                        _resample_linear!(col, t, df[!, "sc$(k)_link_pose_$(base + c)"], times)
                        @inbounds for f in 1:N
                            lp[(f - 1) * lp_total + lp_offsets[s] + base + c] = col[f]
                        end
                    end
                    qsrc = Matrix{Float64}(undef, 4, nrow(df))
                    for c in 1:4
                        qsrc[c, :] .= Float64.(df[!, "sc$(k)_link_pose_$(base + 3 + c)"])
                    end
                    qdst = Matrix{Float64}(undef, 4, N)
                    _resample_quaternions!(qdst, t, qsrc, times)
                    @inbounds for f in 1:N, c in 1:4
                        lp[(f - 1) * lp_total + lp_offsets[s] + base + 3 + c] = qdst[c, f]
                    end
                end
            end
        end
    end

    frames = Dict{String, Any}(
        "count" => N,
        "sats" => S,
        "source_rows" => sum(o -> o.rows, outlines),
        "stride_rows" => 1,
        "t_dtype" => "f64",
        "t_s" => _float64_base64(times),
        "pos_dtype" => pos_f64 ? "f64" : "f32",
        "pos_km" => pos_f64 ? _float64_base64(pos) : _float32_base64(pos),
        "vel_kms" => has_vel ? _float32_base64(vel) : nothing,
        "q" => has_q ? _float32_base64(q) : nothing,
        "mass_kg" => has_mass ? _float32_base64(mass) : nothing,
        "density_kg_m3" => has_density ? _float32_base64(density) : nothing,
        "heat_rate_w_m2" => has_heat ? _float32_base64(heat) : nothing,
        "drag_n" => nothing,
        "wind_ms" => nothing,
        "arm_pose" => nothing,
        "link_pose" => has_lp ? Dict{String, Any}(
            "stride" => stride_lp, "counts" => lp_counts, "offsets" => lp_offsets, "total" => lp_total,
            "data" => _float32_base64(lp)
        ) : nothing,
    )

    # One scene with the geometry repeated per sample so labels, LOD and the
    # info panel address pseudo-spacecraft by index.
    base = scenes[1]
    spacecraft = SpacecraftGeometry[]
    for i in 1:n, k in 1:per_sample
        g = base.spacecraft[k]
        push!(spacecraft, SpacecraftGeometry((i - 1) * per_sample + k, "$(samples[i].label)/$(g.name)", g.links, g.thrusters, g.facets, g.joints, g.bounding_radius_m, g.stl_path, g.arm))
    end
    scene = VisualizationScene(base.schema, base.epoch_et_start_s, base.epoch_utc, base.planet, spacecraft, base.orientation_sim,
                               base.results_feather, base.link_pose_field, base.link_pose_stride, base.atmosphere,
                               base.results_layout)
    return frames, scene, per_sample
end

"""
    export_ensemble_visualization(campaign_dir; out=joinpath(campaign_dir, "ensemble_viewer.html"), kwargs...) -> String

Build one viewer page for every sample found in `campaign_dir` (see
`discover_ensemble_samples`). Samples are coloured by the manifest scalar
(`scalar_name` overrides its label); the page adds a sample selector and the
full-history spaghetti view. Accepts the same page keywords as
`export_visualization` except `stl`.
"""
function export_ensemble_visualization(
    campaign_dir::AbstractString;
    out::Union{Nothing, AbstractString}=nothing,
    max_frames::Integer=DEFAULT_MAX_FRAMES,
    data_budget_mb::Real=DEFAULT_DATA_BUDGET_MB,
    position_precision=:auto,
    scalar_name::Union{Nothing, AbstractString}=nothing,
    trail_s::Union{Nothing, Real}=nothing,
    trail_orbits::Union{Nothing, Real}=nothing,
    frame::Symbol=:inertial,
    speed::Union{Nothing, Real}=nothing,
    title::Union{Nothing, AbstractString}=nothing,
    textures::Bool=true,
    texture_resolution="4k",
    viewer_dir::Union{Nothing, AbstractString}=nothing,
    textures_dir::AbstractString=TEXTURES_DIR
)::String
    _checked_position_precision(position_precision)
    VisualizerServices.require_available()
    campaign_dir = String(campaign_dir)
    samples, manifest_scalar = discover_ensemble_samples(campaign_dir)
    isempty(samples) && throw(ArgumentError("No ensemble samples with a scene sidecar under $(campaign_dir)."))
    finished = EnsembleSample[s for s in samples if isfile(joinpath(s.directory, "simulation_results_scene.json"))]
    isempty(finished) && throw(ArgumentError("None of the ensemble samples under $(campaign_dir) wrote a scene sidecar."))
    scenes = VisualizationScene[read_visualization_scene(joinpath(s.directory, "simulation_results_scene.json")) for s in finished]
    # Each sample's table is read when the frames need it, one at a time.
    tables = [() -> _saved_results_table(scenes[i], joinpath(s.directory, scenes[i].results_feather))
              for (i, s) in enumerate(finished)]
    frames, scene, per_sample = build_ensemble_frames(finished, scenes, tables; max_frames=max_frames, data_budget_mb=data_budget_mb,
        position_precision=position_precision)

    textures_payload = Dict{String, Any}()
    if textures
        entry = texture_payload(scene.planet.texture; resolution=texture_resolution, dir=textures_dir)
        entry === nothing || (textures_payload[scene.planet.texture] = entry)
    end
    name = scalar_name === nothing ? manifest_scalar : String(scalar_name)
    nominal_index = read_ensemble_nominal(campaign_dir)
    nominal_position = nominal_index === nothing ? nothing : findfirst(s -> s.index == nominal_index, finished)
    payload = Dict{String, Any}(
        "viewer_schema" => 1,
        "scene" => scene_dict(scene),
        "frames" => frames,
        "textures" => textures_payload,
        "models" => Dict{String, Any}(),
        "options" => _viewer_options(; trail_s=trail_s, trail_orbits=trail_orbits, frame=frame, speed=speed, title=title),
        "ensemble" => Dict{String, Any}(
            "count" => length(finished),
            "spacecraft_per_sample" => per_sample,
            "scalar_name" => name,
            "samples" => Any[_sample_dict(s) for s in finished],
            "nominal" => nominal_position === nothing ? nothing : nominal_position - 1,
        ),
    )
    page_title = title === nothing ? "SpaceAGORA · $(scene.planet.name) · ensemble of $(length(finished))" : String(title)
    html = render_viewer_html(payload; viewer_dir=viewer_dir, title=page_title)
    out_path = out === nothing ? joinpath(campaign_dir, "ensemble_viewer.html") : String(out)
    protected = String[joinpath(campaign_dir, ENSEMBLE_MANIFEST)]
    for (sample, scene) in zip(finished, scenes)
        push!(protected, joinpath(sample.directory, "simulation_results_scene.json"))
        push!(protected, joinpath(sample.directory, scene.results_feather))
        push!(protected, joinpath(sample.directory, "simulation_results_viewer.json"))
        push!(protected, joinpath(sample.directory, "simulation_results_mission.json"))
        push!(protected, joinpath(sample.directory, "simulation_results.manifest.toml"))
    end
    _check_viewer_output(out_path, protected)
    return IOSerialization._atomic_write_file(out_path, tmp -> write(tmp, html))
end
