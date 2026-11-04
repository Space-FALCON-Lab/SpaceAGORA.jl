module MissionWorkflow

using ..SimulationModel: SimulationConfiguration, SimulationSettings
import ..SimulationModel
import ..SimulationEngine
import SciMLBase
import JSON
import ..SimulationModel.SceneVisualization: export_visualization

export MissionPlan, MissionPhase, MissionEvent, MissionRun, run_mission
export MissionSegment, MissionDeadline, MissionDuration, MissionCondition

function _name(name::Union{Symbol, AbstractString}, kind::AbstractString)
    value = strip(String(name))
    isempty(value) && throw(ArgumentError("$kind needs a nonempty name."))
    return value
end

include(joinpath(@__DIR__, "mission_segments.jl"))

"""
    MissionPhase(name; window_s)

A named observation window `(start_s, end_s)` in elapsed simulation seconds.
Windows must be finite, nonnegative, increasing and nonoverlapping in a
[`MissionPlan`](@ref). They describe portions of one solve; they do not switch
models, reset controllers, insert solver stops or execute separate integrations.
The run may stop before a window ends because of an existing stop condition.
"""
struct MissionPhase
    name::Symbol
    window_s::Tuple{Float64, Float64}
    function MissionPhase(name::Union{Symbol, AbstractString}; window_s)
        window_s isa Tuple && length(window_s) == 2 && all(x -> x isa Real, window_s) ||
            throw(ArgumentError("MissionPhase.window_s must be a pair of real elapsed times."))
        start_s, end_s = Float64.(window_s)
        isfinite(start_s) && isfinite(end_s) && 0.0 <= start_s < end_s ||
            throw(ArgumentError("MissionPhase.window_s must satisfy 0 <= start_s < end_s with finite times."))
        return new(Symbol(_name(name, "MissionPhase")), (start_s, end_s))
    end
end

"""
    MissionEvent(name, build_callback)

A named factory `build_callback(configuration)` that returns a fresh solver
callback (or callback set) for each run. It receives the configuration after
SpaceAGORA has aligned epochs, flattened spacecraft GNC and isolated run state.
Allocate captured mutable callback state inside the factory. Do not return a
shared callback, mutate the configuration during construction, or perform I/O.

The callback keeps its own trigger and action semantics. A state anchor replaces
state, a maneuver changes it according to its existing model, and a discrete
stop checks accepted steps. A declared event name does not assert that it fired.
"""
struct MissionEvent{F <: Function}
    name::Symbol
    build_callback::F
    function MissionEvent(name::Union{Symbol, AbstractString}, build_callback::F) where {F <: Function}
        return new{F}(Symbol(_name(name, "MissionEvent")), build_callback)
    end
end

"""
    MissionPlan(name, build_configuration; phases=(), events=(), segments=())

A reusable mission definition. `build_configuration(; parameters...)` must
construct and return a `SimulationConfiguration`; `run_mission` calls it once
per run. Define spacecraft, controllers and other mutable model state inside
that builder. The engine retains its normal isolation as an additional boundary.

`phases` contains ordered [`MissionPhase`](@ref) observation windows; `events`
contains [`MissionEvent`](@ref) callback factories. Names must be unique within
each collection. This plan neither changes scientific configuration defaults
nor replaces the existing `MissionConfiguration` duration and stop settings.
`segments` contains ordered operational `MissionSegment` values. These use one
continuous solve and explicit command actions at timed or condition boundaries.
"""
struct MissionPlan{F <: Function, P <: Tuple, E <: Tuple, S <: Tuple}
    name::String
    build_configuration::F
    phases::P
    events::E
    segments::S
    function MissionPlan(name::AbstractString, build_configuration::F; phases=(), events=(), segments=()) where {F <: Function}
        phases isa Union{Tuple, AbstractVector} || throw(ArgumentError("MissionPlan.phases must be a tuple or vector of MissionPhase values."))
        events isa Union{Tuple, AbstractVector} || throw(ArgumentError("MissionPlan.events must be a tuple or vector of MissionEvent values."))
        segments isa Union{Tuple, AbstractVector} || throw(ArgumentError("MissionPlan.segments must be a tuple or vector of MissionSegment values."))
        phase_tuple, event_tuple, segment_tuple = Tuple(phases), Tuple(events), Tuple(segments)
        all(s -> s isa MissionSegment, segment_tuple) || throw(ArgumentError("MissionPlan.segments must contain only MissionSegment values."))
        length(unique(map(s -> s.name, segment_tuple))) == length(segment_tuple) || throw(ArgumentError("MissionSegment names must be unique."))
        all(s -> s.until !== nothing, segment_tuple[1:max(0, length(segment_tuple) - 1)]) || throw(ArgumentError("Only the final MissionSegment may omit until."))
        all(p -> p isa MissionPhase, phase_tuple) || throw(ArgumentError("MissionPlan.phases must contain only MissionPhase values."))
        all(e -> e isa MissionEvent, event_tuple) || throw(ArgumentError("MissionPlan.events must contain only MissionEvent values."))
        phase_names = map(p -> p.name, phase_tuple)
        event_names = map(e -> e.name, event_tuple)
        length(unique(phase_names)) == length(phase_names) || throw(ArgumentError("MissionPhase names must be unique."))
        length(unique(event_names)) == length(event_names) || throw(ArgumentError("MissionEvent names must be unique."))
        for i in 2:length(phase_tuple)
            phase_tuple[i - 1].window_s[2] <= phase_tuple[i].window_s[1] ||
                throw(ArgumentError("MissionPhase windows must be ordered and nonoverlapping."))
        end
        return new{F, typeof(phase_tuple), typeof(event_tuple), typeof(segment_tuple)}(
            _name(name, "MissionPlan"), build_configuration, phase_tuple, event_tuple, segment_tuple)
    end
end

"""
    MissionRun

The result of [`run_mission`](@ref). `results` is the existing `SimulationResults`
with its table, isolated configuration and optional solution. `directory` and
`prefix` identify the saved run; `payload_path` identifies the standalone viewer
JSON when requested. `files` lists the files written by the simulation and its
payload export. `phases` contains descriptive windows, and `event_names` lists
configured event factories rather than an event-occurrence log. `history` holds
actual operational segment intervals and transitions; `history_path` names their
JSON sidecar, or is `nothing` for a plan without segments.

Paths refer to the original destination. After moving a run directory, pass its
new path to the standalone viewer. No live simulation or Julia object is required
to view the saved payload. `export_visualization(run)` instead builds the existing
HTML viewer from the bundle at the handle's original prefix.
"""
struct MissionRun{R <: SimulationEngine.SimulationResults, P <: Tuple, E <: Tuple}
    name::String
    results::R
    directory::String
    prefix::String
    payload_path::Union{Nothing, String}
    phases::P
    event_names::E
    files::Vector{String}
    history::typeof(_empty_segment_history())
    history_path::Union{Nothing, String}
end

function _destination(output_dir::AbstractString)
    isempty(strip(output_dir)) && throw(ArgumentError("run_mission needs a nonempty output_dir."))
    dir = abspath(output_dir)
    islink(dir) && throw(ArgumentError("run_mission output_dir must not be a symbolic link: $dir"))
    if ispath(dir)
        isdir(dir) || throw(ArgumentError("run_mission output_dir must name a directory: $dir"))
        isempty(readdir(dir)) || throw(ArgumentError("run_mission needs an absent or empty output directory; destination is populated: $dir"))
    end
    return dir
end

function _configuration(plan::MissionPlan, parameters::NamedTuple, destination::String, save_scene::Bool)
    config = plan.build_configuration(; parameters...)
    config isa SimulationConfiguration || throw(ArgumentError("MissionPlan builder must return a SimulationConfiguration; got $(typeof(config))."))
    settings = config.simulation_settings
    (settings.checkpoint_enabled || settings.resume_from_checkpoint) && throw(ArgumentError(
        "run_mission does not support checkpoint writing or resume yet; use run_simulation for checkpoint workflows."))
    horizon = config.mission_configuration.mission_time
    isfinite(horizon) && horizon > 0 || throw(ArgumentError("run_mission needs a finite positive mission_time."))
    all(phase -> phase.window_s[2] <= horizon, plan.phases) || throw(ArgumentError(
        "MissionPhase observation windows must fit within mission_configuration.mission_time."))
    _validate_segments(plan.segments, config)
    # These three output choices are the complete wrapper override. Carry every
    # other setting by field name, including future scientific settings.
    names = fieldnames(typeof(settings))
    values = NamedTuple{names}(map(n -> getfield(settings, n), names))
    output_settings = SimulationSettings(; merge(values, (
        results=true, results_directory=destination, save_visualization_scene=save_scene))...)
    return SimulationModel.SimConfig._with_configuration(config; simulation_settings=output_settings)
end

"""
    run_mission(plan; output_dir, parameters=(;), save_scene=true,
                save_payload=true, return_solution=false, save_fields=nothing)

Build and run a mission once, returning a [`MissionRun`](@ref). The destination
must be absent or empty. This output-oriented workflow enables saved results and
selects `output_dir`; `save_scene` also selects saved scene/pose output. All
physical models, epoch, solver settings, tolerances and sampling settings come
from the builder unchanged. The normal engine state isolation remains enabled.

By default a scene and renderer-independent viewer JSON without textures are written.
Texture-rich payloads can be exported later without rerunning the simulation. Set both
`save_scene=false` and `save_payload=false` for table/bundle output only.
`save_payload=true` requires `save_scene=true`; scene output also requires result
bundles to be enabled. The workflow does not launch a browser or produce HTML.
Use `export_visualization(run)` explicitly for the existing HTML export.

Plan, builder, destination and factory errors are checked before SpaceAGORA
creates run outputs. User builders/factories must themselves avoid I/O. Runtime
failures retain the engine's existing partial-output behavior and are rethrown;
a failed run is not returned as a successful handle. Checkpoint writing and
resume are intentionally unsupported here; use `run_simulation` for them.
"""
function run_mission(plan::MissionPlan; output_dir::AbstractString,
    parameters::NamedTuple=(;), save_scene::Bool=true, save_payload::Bool=true,
    return_solution::Bool=false, save_fields=nothing)
    save_payload && !save_scene && throw(ArgumentError("save_payload=true requires save_scene=true."))
    save_scene && !SimulationEngine._typed_save_bundle_enabled() && throw(ArgumentError(
        "Mission scene output requires result bundles; enable SPACEAGORA_SAVE_BUNDLE or disable save_scene and save_payload."))
    destination = _destination(output_dir)
    config = _configuration(plan, parameters, destination, save_scene)
    _destination(destination)
    factories = map(event -> event.build_callback, plan.events)
    runtime_ref = Ref{Any}(nothing)
    priority_factories = isempty(plan.segments) ? () :
        (configuration -> _segment_callbacks(plan.segments, configuration, runtime_ref),)
    results = SimulationEngine.run_simulation(config; return_results=true,
        return_solution=return_solution, save_fields=save_fields,
        extra_callback_factories=factories, priority_callback_factories=priority_factories, visualization=false)
    prefix = SimulationModel.IOConfig._results_bundle_prefix(results.configuration)
    history = runtime_ref[] === nothing ? _empty_segment_history() : runtime_ref[].history
    history_path = isempty(plan.segments) ? nothing : _write_mission_history(prefix, plan.name, history)
    if history_path !== nothing && save_scene
        # Associate only outputs this run actually produced, never an old scene
        # or manifest surviving in a reused output directory.
        all(path -> path in results.files, (prefix * ".manifest.toml", prefix * "_scene.json")) ||
            throw(ArgumentError("Mission history needs this run's saved scene and result manifest."))
        SimulationModel.SceneVisualization._bind_mission_history!(prefix)
    end
    payload_path = save_payload ? SimulationModel.SceneVisualization.export_viewer_payload(prefix; textures=false) : nothing
    files = copy(results.files)
    payload_path === nothing || push!(files, payload_path)
    history_path === nothing || push!(files, history_path)
    return MissionRun(plan.name, results, destination, prefix, payload_path,
        plan.phases, map(event -> event.name, plan.events), unique(files), history, history_path)
end

"""
    export_visualization(run::MissionRun; kwargs...) -> String

Build the existing HTML viewer from this run's saved scene and results bundle.
This is an explicit post-processing step; it never reruns the mission.
The standalone viewer can instead read `run.payload_path` or `run.directory`.
"""
function export_visualization(run::MissionRun; kwargs...)
    return export_visualization(run.prefix; kwargs...)
end

end # module MissionWorkflow
