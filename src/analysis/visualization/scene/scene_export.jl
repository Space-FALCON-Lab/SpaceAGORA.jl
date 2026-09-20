# Scene assembly from a SimulationConfiguration, the JSON sidecar format, and
# the playback frame budget the exporter will use.

"""
    visualization_scene_path(args) -> String

Sidecar path next to the results bundle: `<results_directory>/simulation_results_scene.json`.
"""
@inline function visualization_scene_path(args::SimulationConfiguration)::String
    return IOConfig._results_bundle_prefix(args) * "_scene.json"
end

@inline function _epoch_utc_string(initial_time)::String
    return Dates.format(_initial_time_datetime(initial_time), dateformat"yyyy-mm-dd\THH:MM:SS.sss") * "Z"
end

"""
    build_visualization_scene(args; rotation_max_samples=4096, stl_paths=Dict()) -> VisualizationScene

Assemble the scene for a configuration: planet with rotation table over the
mission span, one `SpacecraftGeometry` per spacecraft, and the epoch.
`stl_paths` maps a spacecraft id to a CAD file recorded as that spacecraft's
`stl_path` (the viewer bundler embeds it when present).
"""
function build_visualization_scene(
    args::SimulationConfiguration;
    rotation_max_samples::Integer=ROTATION_MAX_SAMPLES,
    stl_paths::AbstractDict=Dict{Int, String}(),
    density_params=nothing
)::VisualizationScene
    ephemerides_model = args.environment_model.ephemerides_model
    planet = args.environment_model.planet
    et_start = ephemerides_time_seconds(args.initial_time, ephemerides_model)
    mission = args.mission_configuration
    sample_times = rotation_sample_times(mission.mission_time, mission.data_rate; max_samples=rotation_max_samples)
    spacecraft = SpacecraftGeometry[
        spacecraft_geometry(model; stl_path=get(stl_paths, Int(model.id), nothing), arm=robot_arm_plan_for(args, i))
        for (i, model) in enumerate(args.dynamics_model.spacecraft)
    ]
    return VisualizationScene(
        SCENE_SCHEMA_VERSION,
        et_start,
        _epoch_utc_string(args.initial_time),
        planet_spec(planet, ephemerides_model, et_start, sample_times),
        spacecraft,
        mission.orientation_sim,
        basename(IOConfig._results_bundle_prefix(args)) * ".feather",
        String(LINK_POSE_FIELD),
        LINK_POSE_STRIDE,
        atmosphere_spec(args; density_params=density_params)
    )
end

# ---------------------------------------------------------------------------
# JSON encoding
# ---------------------------------------------------------------------------

@inline _plain(v::SVector) = collect(Float64, v)

function _scene_dict(box::LinkBox)
    return Dict{String, Any}(
        "name" => box.name,
        "root" => box.root,
        "dims_m" => _plain(box.dims_m),
        "r_m" => _plain(box.r_m),
        "q" => _plain(box.q),
        "mass_kg" => box.mass_kg
    )
end

_scene_dict(g::ThrusterGlyph) = Dict{String, Any}(
    "link" => g.link, "location_m" => _plain(g.location_m), "direction" => _plain(g.direction), "max_thrust_n" => g.max_thrust_n
)
_scene_dict(g::FacetGlyph) = Dict{String, Any}(
    "link" => g.link, "name" => g.name, "normal" => _plain(g.normal), "cp_m" => _plain(g.cp_m), "area_m2" => g.area_m2
)
_scene_dict(g::JointGlyph) = Dict{String, Any}(
    "link1" => g.link1, "link2" => g.link2, "p1_m" => _plain(g.p1_m), "p2_m" => _plain(g.p2_m)
)

function _scene_dict(sc::SpacecraftGeometry)
    return Dict{String, Any}(
        "id" => sc.id,
        "name" => sc.name,
        "links" => Any[_scene_dict(b) for b in sc.links],
        "thrusters" => Any[_scene_dict(g) for g in sc.thrusters],
        "facets" => Any[_scene_dict(g) for g in sc.facets],
        "joints" => Any[_scene_dict(g) for g in sc.joints],
        "bounding_radius_m" => sc.bounding_radius_m,
        "stl_path" => sc.stl_path,
        "arm" => sc.arm === nothing ? nothing : _scene_dict(sc.arm)
    )
end

function _scene_dict(p::PlanetSpec)
    return Dict{String, Any}(
        "name" => p.name,
        "equatorial_radius_m" => p.equatorial_radius_m,
        "polar_radius_m" => p.polar_radius_m,
        "spin_rad_s" => _plain(p.spin_rad_s),
        "inertial_frame" => p.inertial_frame,
        "texture" => p.texture,
        "rotation" => Dict{String, Any}(
            "t_s" => copy(p.rotation_t_s),
            "q_pi" => Any[_plain(q) for q in p.rotation_q_pi]
        )
    )
end

"""
    scene_dict(scene) -> Dict{String, Any}

JSON-ready form of a `VisualizationScene`; `read_visualization_scene` inverts it.
"""
function scene_dict(scene::VisualizationScene)::Dict{String, Any}
    return Dict{String, Any}(
        "schema" => scene.schema,
        "epoch" => Dict{String, Any}("et_start_s" => scene.epoch_et_start_s, "utc" => scene.epoch_utc),
        "planet" => _scene_dict(scene.planet),
        "spacecraft" => Any[_scene_dict(sc) for sc in scene.spacecraft],
        "orientation_sim" => scene.orientation_sim,
        "atmosphere" => scene.atmosphere === nothing ? nothing : _scene_dict(scene.atmosphere),
        "results" => Dict{String, Any}(
            "feather" => scene.results_feather,
            "link_pose" => Dict{String, Any}(
                "field" => scene.link_pose_field,
                "stride" => scene.link_pose_stride,
                "layout" => collect(String, LINK_POSE_LAYOUT),
                "columns" => "sc{i}_$(scene.link_pose_field)_{1..$(scene.link_pose_stride)n}"
            ),
            "arm_pose" => Dict{String, Any}(
                "field" => String(ARM_POSE_FIELD),
                "stride" => LINK_POSE_STRIDE,
                "layout" => collect(String, LINK_POSE_LAYOUT),
                "frame" => "inertial, position relative to the spacecraft in meters",
                "columns" => "sc{i}_$(ARM_POSE_FIELD)_{1..$(LINK_POSE_STRIDE)n}"
            )
        )
    )
end

_link_box_from(d) = LinkBox(String(d["name"]), Bool(d["root"]), _svec3(d["dims_m"]), _svec3(d["r_m"]), _svec4(d["q"]), Float64(d["mass_kg"]))
_thruster_from(d) = ThrusterGlyph(Int(d["link"]), _svec3(d["location_m"]), _svec3(d["direction"]), Float64(d["max_thrust_n"]))
_facet_from(d) = FacetGlyph(Int(d["link"]), String(d["name"]), _svec3(d["normal"]), _svec3(d["cp_m"]), Float64(d["area_m2"]))
_joint_from(d) = JointGlyph(Int(d["link1"]), Int(d["link2"]), _svec3(d["p1_m"]), _svec3(d["p2_m"]))

_scene_dict(l::ArmLinkGeometry) = Dict{String, Any}(
    "name" => l.name, "vector_m" => _plain(l.vector_m), "com_offset_m" => _plain(l.com_offset_m), "radius_m" => l.radius_m, "mass_kg" => l.mass_kg
)
_scene_dict(a::ArmGeometry) = Dict{String, Any}(
    "links" => Any[_scene_dict(l) for l in a.links], "mount_offset_body_m" => _plain(a.mount_offset_body_m), "reach_m" => a.reach_m
)
_arm_link_from(d) = ArmLinkGeometry(String(d["name"]), _svec3(d["vector_m"]), _svec3(d["com_offset_m"]), Float64(d["radius_m"]), Float64(d["mass_kg"]))
_arm_from(d) = d === nothing ? nothing : ArmGeometry(ArmLinkGeometry[_arm_link_from(x) for x in d["links"]], _svec3(d["mount_offset_body_m"]), Float64(d["reach_m"]))

function _spacecraft_from(d)::SpacecraftGeometry
    stl = get(d, "stl_path", nothing)
    return SpacecraftGeometry(
        Int(d["id"]),
        String(d["name"]),
        LinkBox[_link_box_from(x) for x in d["links"]],
        ThrusterGlyph[_thruster_from(x) for x in d["thrusters"]],
        FacetGlyph[_facet_from(x) for x in d["facets"]],
        JointGlyph[_joint_from(x) for x in d["joints"]],
        Float64(d["bounding_radius_m"]),
        stl === nothing ? nothing : String(stl),
        _arm_from(get(d, "arm", nothing))
    )
end

function _planet_from(d)::PlanetSpec
    rotation = d["rotation"]
    return PlanetSpec(
        String(d["name"]),
        Float64(d["equatorial_radius_m"]),
        Float64(d["polar_radius_m"]),
        _svec3(d["spin_rad_s"]),
        String(d["inertial_frame"]),
        String(d["texture"]),
        Float64[Float64(t) for t in rotation["t_s"]],
        SVector{4, Float64}[_svec4(q) for q in rotation["q_pi"]]
    )
end

function _scene_from_dict(d)::VisualizationScene
    schema = Int(d["schema"])
    schema == SCENE_SCHEMA_VERSION || throw(ArgumentError("Unsupported visualization scene schema $(schema); expected $(SCENE_SCHEMA_VERSION)."))
    results = d["results"]
    link_pose = results["link_pose"]
    return VisualizationScene(
        schema,
        Float64(d["epoch"]["et_start_s"]),
        String(d["epoch"]["utc"]),
        _planet_from(d["planet"]),
        SpacecraftGeometry[_spacecraft_from(x) for x in d["spacecraft"]],
        Bool(d["orientation_sim"]),
        String(results["feather"]),
        String(link_pose["field"]),
        Int(link_pose["stride"]),
        _atmosphere_from(get(d, "atmosphere", nothing))
    )
end

"""
    write_visualization_scene(path, scene) -> String

Atomically write the scene sidecar as JSON and return `path`.
"""
function write_visualization_scene(path::AbstractString, scene::VisualizationScene)::String
    payload = scene_dict(scene)
    return IOSerialization._atomic_write_file(String(path), tmp -> open(io -> JSON.print(io, payload), tmp, "w"))
end

"""
    read_visualization_scene(path) -> VisualizationScene

Parse a sidecar written by `write_visualization_scene`.
"""
function read_visualization_scene(path::AbstractString)::VisualizationScene
    return _scene_from_dict(JSON.parsefile(String(path)))
end

"""
    with_visualization_scene(args, flag=true) -> SimulationConfiguration

Copy of `args` with `simulation_settings.save_visualization_scene` set to
`flag`; every other field is carried over unchanged.
"""
function with_visualization_scene(args::SimulationConfiguration, flag::Bool=true)::SimulationConfiguration
    s = args.simulation_settings
    s.save_visualization_scene == flag && return args
    settings = SimulationSettings(
        results=s.results,
        verbose=s.verbose,
        results_directory=s.results_directory,
        generate_plots=s.generate_plots,
        generate_filenames=s.generate_filenames,
        normalize=s.normalize,
        save_csv=s.save_csv,
        save_visualization_scene=flag,
        checkpoint_enabled=s.checkpoint_enabled,
        checkpoint_interval_s=s.checkpoint_interval_s,
        checkpoint_directory=s.checkpoint_directory,
        resume_from_checkpoint=s.resume_from_checkpoint
    )
    return _with_configuration(args;
        simulation_settings=settings,
    )
end

"""
    write_visualization_scene!(args; density_params=nothing) -> Union{Nothing, String}

Write the sidecar for `args` when `simulation_settings.save_visualization_scene`
is set; otherwise do nothing and return `nothing`. Called by the engine after
the results bundle is written, with its integrator parameters as
`density_params`. Verified pure analytic atmosphere models and fixed grid
snapshots are sampled automatically within their coverage; native and user
models keep their shell and saved trajectory density without additional calls.
"""
function write_visualization_scene!(args::SimulationConfiguration; density_params=nothing)::Union{Nothing, String}
    args.simulation_settings.save_visualization_scene || return nothing
    return write_visualization_scene(visualization_scene_path(args), build_visualization_scene(args; density_params=density_params))
end

# ---------------------------------------------------------------------------
# Playback budget (used by the viewer exporter; pure arithmetic, tested here)
# ---------------------------------------------------------------------------

"""
    visualization_frame_budget(n_rows, num_sats; max_frames=2000, data_budget_mb=150.0, bytes_per_sat_frame=12) -> (frames, stride, bytes)

How many result rows the exporter keeps for playback: the smaller of
`max_frames` and what `data_budget_mb` allows at `bytes_per_sat_frame` per
spacecraft per frame, never fewer than two when two rows exist. `stride` is
the row step that realizes it and `bytes` the resulting payload size.
"""
function visualization_frame_budget(
    n_rows::Integer,
    num_sats::Integer;
    max_frames::Integer=2000,
    data_budget_mb::Real=150.0,
    bytes_per_sat_frame::Integer=12
)
    n_rows >= 0 || throw(ArgumentError("n_rows must be non-negative, got $(n_rows)."))
    num_sats >= 1 || throw(ArgumentError("num_sats must be at least 1, got $(num_sats)."))
    max_frames >= 2 || throw(ArgumentError("max_frames must be at least 2, got $(max_frames)."))
    data_budget_mb > 0 || throw(ArgumentError("data_budget_mb must be positive, got $(data_budget_mb)."))
    bytes_per_sat_frame >= 1 || throw(ArgumentError("bytes_per_sat_frame must be at least 1, got $(bytes_per_sat_frame)."))
    n_rows == 0 && return (frames=0, stride=1, bytes=0)
    per_frame = Int(num_sats) * Int(bytes_per_sat_frame)
    budget_frames = floor(Int, Float64(data_budget_mb) * 1.0e6 / per_frame)
    frames = min(Int(n_rows), Int(max_frames), budget_frames)
    frames = max(frames, min(Int(n_rows), 2))
    stride = cld(Int(n_rows), frames)
    frames = cld(Int(n_rows), stride)
    return (frames=frames, stride=stride, bytes=frames * per_frame)
end
