# The run state a checkpoint holds next to the integrated state, written with every
# checkpoint (`_write_checkpoint!`) and read back when a run resumes. Its entries,
# in the order a resumed run uses them (see "Run state: owners, resets and
# restores" in docs/contributing/architecture.md):
#
# - `activity`, `interlinks`, `callbacks`, `maneuver_buffers`, `orbits` and `entries`: the
#   engine's and the callbacks' run state, restored by `_resume_engine_state!`
#   (execution.jl) before the models start and before the callbacks initialize;
# - `models`: the models' own run state, checked when the checkpoint is read and
#   restored by `_start_run_models!` after the models' initializers;
# - `gnc_scopes`: the spacecraft each guidance, navigation and control model acts
#   on, compared by `_start_run_models!` after the initializers and before the
#   models' run state is restored;
# - `periodic_ticks`: the periodic ticks not yet taken at the checkpoint's time,
#   which a resumed run takes at its start (`_start_run_clock!`);
# - `results_rows`: how much of the results rows kept beside the checkpoint
#   belongs to it, read before the run's first segment (`_start_results_rows!`);
# - `interlink_history`: likewise, how much of the interlink history kept beside
#   the checkpoint belongs to it (`_start_interlink_history!`).
#
# The definitions below follow that order, then the segment ends and the file.

@inline _checkpoint_directory(args) = SimulationModel.IOConfig._checkpoint_directory(args)
@inline _checkpoint_paths(args) = SimulationModel.IOConfig._checkpoint_paths(args)

# Which spacecraft are active: an impact or a touchdown deactivates one for the
# rest of the run. Written with every checkpoint and restored before a resumed
# run's initializers, so a spacecraft inactive at the checkpoint stays inactive:
# its state is not integrated, the GNC models that skip inactive spacecraft do
# not tick for it, and its impact and touchdown do not fire again (both act on
# active spacecraft only, and the touchdown's start check skips inactive ones).
# An interlink checkpoint also carries the flags in its own entry; the two must
# agree. A checkpoint without the entry, written by an earlier SpaceAGORA,
# resumes with every spacecraft active, except in an interlinked run, which
# restores its entry's flags.
_activity_checkpoint_state(is_active) = (version=1, is_active=Vector{Bool}(is_active))

function _restore_activity_checkpoint_state!(args, checkpoint, is_active)
    runtime = get(checkpoint, :runtime_state, nothing)
    entry = runtime !== nothing && hasproperty(runtime, :activity) ? runtime.activity : nothing
    entry === nothing && return nothing
    hasproperty(entry, :version) && entry.version === 1 ||
        throw(ArgumentError("Spacecraft activity checkpoint state version 1 is required to resume this run."))
    count = length(args.dynamics_model.spacecraft)
    hasproperty(entry, :is_active) && entry.is_active isa AbstractVector{Bool} &&
        length(entry.is_active) == count && length(is_active) == count ||
        throw(ArgumentError("The checkpoint's spacecraft activity flags must match the spacecraft count."))
    interlinks = hasproperty(runtime, :interlinks) ? runtime.interlinks : nothing
    interlinks === nothing || !hasproperty(interlinks, :is_active) || interlinks.is_active == entry.is_active ||
        throw(ArgumentError("The checkpoint's spacecraft activity flags differ from its interlink entry's."))
    copyto!(is_active, entry.is_active)
    return nothing
end

# The burn scheduler's command channel (`SharedBuffers`): the maneuver command each
# spacecraft's maneuver guidance last issued, and the burn plan installed for it.
# A run carries both from one checkpoint segment to the next, so a resumed run
# takes them from its checkpoint: a command issued before the checkpoint and not
# yet planned, or a plan not yet flown, is then planned and flown as in the run
# without the break. Written only when some entry differs from its default, so the
# checkpoints of other runs do not change, and a checkpoint without the entry,
# written by an earlier SpaceAGORA, resumes with the defaults, as before.
function _maneuver_buffer_checkpoint_state(p)
    buffers = p.shared_buffers
    commands, plans = buffers.maneuver_commands, buffers.maneuver_burn_plans
    all(c -> isequal(c, SimulationModel.PropulsiveManeuverCommand()), commands) &&
        all(b -> isequal(b, SimulationModel.PropulsiveBurnPlan()), plans) && return nothing
    return (version=1, commands=copy(commands), plans=copy(plans))
end

function _restore_maneuver_buffer_checkpoint_state!(checkpoint, p)
    runtime = get(checkpoint, :runtime_state, nothing)
    entry = runtime !== nothing && hasproperty(runtime, :maneuver_buffers) ? runtime.maneuver_buffers : nothing
    entry === nothing && return nothing
    hasproperty(entry, :version) && entry.version === 1 ||
        throw(ArgumentError("Maneuver command checkpoint state version 1 is required to resume this run."))
    buffers = p.shared_buffers
    length(entry.commands) == length(buffers.maneuver_commands) &&
        length(entry.plans) == length(buffers.maneuver_burn_plans) ||
        throw(ArgumentError("The checkpoint's maneuver commands and burn plans must match the spacecraft count."))
    copyto!(buffers.maneuver_commands, entry.commands)
    copyto!(buffers.maneuver_burn_plans, entry.plans)
    return nothing
end

# The orbit each spacecraft is on (`p.orbit_counter`, from 1; the apoapsis callback
# adds one at each apoapsis). The orbit-count mission end, the campaign maneuver
# guidance, pass-indexed density models and periapsis pulses read it, so a resumed
# run continues the count. Written once some spacecraft is past its first orbit, so
# the checkpoints of shorter runs do not change, and a checkpoint without the
# entry, written by an earlier SpaceAGORA, resumes from the first orbit, as before.
function _orbit_counter_checkpoint_state(p)
    all(==(1), p.orbit_counter) && return nothing
    return (version=1, orbit_counter=copy(p.orbit_counter))
end

function _restore_orbit_counter_checkpoint_state!(checkpoint, p)
    runtime = get(checkpoint, :runtime_state, nothing)
    entry = runtime !== nothing && hasproperty(runtime, :orbits) ? runtime.orbits : nothing
    entry === nothing && return nothing
    hasproperty(entry, :version) && entry.version === 1 ||
        throw(ArgumentError("Orbit counter checkpoint state version 1 is required to resume this run."))
    length(entry.orbit_counter) == length(p.orbit_counter) ||
        throw(ArgumentError("The checkpoint's orbit counters must match the spacecraft count."))
    copyto!(p.orbit_counter, entry.orbit_counter)
    return nothing
end

# Atmosphere entries counted by the entry-end callback. These used to live only
# in its closure, which survived solve segments but was rebuilt on resume. Keep
# them with the run buffers so a checkpoint can restore them before callbacks run.
# An older checkpoint without this entry retains the historical zero-count start.
function _entry_counter_checkpoint_state(p)
    counts = p.shared_buffers.entry_counter
    all(iszero, counts) && return nothing
    return (version=1, entry_counter=copy(counts))
end

function _restore_entry_counter_checkpoint_state!(checkpoint, p)
    runtime = get(checkpoint, :runtime_state, nothing)
    entry = runtime !== nothing && hasproperty(runtime, :entries) ? runtime.entries : nothing
    entry === nothing && return nothing
    hasproperty(entry, :version) && entry.version === 1 ||
        throw(ArgumentError("Atmosphere entry counter checkpoint state version 1 is required to resume this run."))
    counts = p.shared_buffers.entry_counter
    hasproperty(entry, :entry_counter) && entry.entry_counter isa AbstractVector{<:Integer} &&
        length(entry.entry_counter) == length(counts) && all(>=(0), entry.entry_counter) ||
        throw(ArgumentError("The checkpoint's atmosphere entry counters must be nonnegative and match the spacecraft count."))
    copyto!(counts, entry.entry_counter)
    return nothing
end

# Keep feature identity and lifecycle knowledge in the simulation engine. The
# file layer stores this plain-data snapshot without importing interlink models.
# A link's identity is its kind and its fields. The kind is named by the full
# path of the module that defines it and its name, so same-named kinds of two
# packages, or of two modules of one package, differ; a kind moved to another
# module needs its old checkpoints migrated. Its concrete type parameters are
# part of it, since they can select different hooks (`ModeLink{1}` and
# `ModeLink{-1}`).
function _interlink_link_identity(link)
    kind = typeof(link)
    names = fieldnames(kind)
    return (kind=_interlink_kind_identity(kind),
        fields=NamedTuple{names}(map(name -> getfield(link, name), names)))
end

# A type as plain data that the checkpoint file serializes: the full path of
# its module, its name and its parameters. A parameter that is a type is
# described the same way, and a tuple parameter element by element. A value
# parameter (a number, symbol or other bits value, the only values a type
# parameter can hold) is recorded with its own type, since `1`, `1.0`, `true`
# and `Int32(1)` select different methods yet compare equal as numbers. A
# floating-point value also keeps its bits, so that -0.0 and 0.0, and NaNs with
# different payloads, stay apart as they do in types. A composite bits value (a
# `ComplexF64`, a struct of numbers) is recorded field by field, the same way:
# its own `isequal` can merge values that are distinct types (two NaN payloads
# in one `Complex`), and its padding bytes carry no meaning.
_interlink_kind_identity(kind::DataType) = (module_path=fullname(parentmodule(kind)), name=nameof(kind),
    parameters=map(_interlink_kind_parameter, Tuple(kind.parameters)))
_interlink_kind_parameter(parameter::DataType) = _interlink_kind_identity(parameter)
_interlink_kind_parameter(parameter::Union) = (union=map(_interlink_kind_parameter, Tuple(Base.uniontypes(parameter))),)
_interlink_kind_parameter(parameter::Type) = string(parameter)   # a UnionAll parameter: its printed form
_interlink_kind_parameter(parameter::Tuple) = map(_interlink_kind_parameter, parameter)
function _interlink_kind_parameter(parameter)
    kind = typeof(parameter)
    isbitstype(kind) && !isprimitivetype(kind) && return (type=_interlink_kind_identity(kind),
        fields=ntuple(index -> _interlink_kind_parameter(getfield(parameter, index)), fieldcount(kind)))
    return (type=_interlink_kind_identity(kind), value=parameter)
end
_interlink_kind_parameter(parameter::Base.IEEEFloat) = (type=_interlink_kind_identity(typeof(parameter)),
    value=parameter, bits=reinterpret(Base.uinttype(typeof(parameter)), parameter))

function _interlink_checkpoint_configuration(args)
    model = args.interlink_model
    graph = model.graph
    configuration = (;
        vertices=copy(graph.vertices),
        candidates=[(
            endpoints=map(vertex -> graph.vertices[vertex], edge),
            link=_interlink_link_identity(graph.edge_properties[index]),
        ) for (index, edge) in enumerate(graph.edges)],
        forbidden_pairs=sort!(collect(model.forbidden_pairs)),
        _scheduling_policy_entry(args.scheduling_policy_model.objective)...,
        terminal_counts=[spacecraft.n_terminal for spacecraft in args.dynamics_model.spacecraft],
    )
    # Terminal numbers past `n_terminal` name mounted terminal hardware, so a run
    # with any records each spacecraft's terminals, in order; a run without keeps
    # the configuration it had before hardware terminals.
    hardware = [_mounted_terminal_identities(spacecraft) for spacecraft in args.dynamics_model.spacecraft]
    all(isempty, hardware) || (configuration = merge(configuration, (hardware_terminals=hardware,)))
    # Likewise, a run with ground stations records each station's site, mask and
    # antennas, in order.
    isempty(model.ground_stations) ||
        (configuration = merge(configuration, (ground_stations=map(_ground_station_identity, model.ground_stations),)))
    mounts = _interlink_terminal_mounts(args)
    isempty(mounts) && return configuration
    return merge(configuration, (terminal_mounts=mounts,))
end

# Geometry-aware candidates also depend on the terminal's mounting link. Keep
# this separate from the hardware identity so center-to-center links retain their
# historical checkpoint configuration. Only endpoints that read geometry matter.
function _interlink_terminal_mounts(args)
    graph = args.interlink_model.graph
    models = SimulationModel.InterLinkModels
    endpoints = Set{Tuple{Int, Int}}()
    for (edge, link) in zip(graph.edges, graph.edge_properties)
        (models._link_needs_attitude(link) || models.link_applies_torque(link)) || continue
        for vertex in edge
            graph.vertices[vertex][1] > 0 && push!(endpoints, graph.vertices[vertex])
        end
    end
    mounts = Any[]
    for (position, number) in sort!(collect(endpoints))
        spacecraft = args.dynamics_model.spacecraft[position]
        for entry in SimulationModel.mounted_link_terminals(spacecraft)
            entry.number == number || continue
            link = entry.link
            index = findfirst(candidate -> candidate === link, spacecraft.links)
            # Root attitude belongs to integrated spacecraft state, not the
            # link-to-body transform. A hinge's angle is restored separately;
            # comparing its current q would reject an unchanged moving mount.
            rotation = if link.root
                (kind=:body_frame,)
            elseif isempty(link.hinges)
                for control in args.control_model.control_effectors
                    panel = control isa SimulationModel.AerobrakingEnergyDepletionControlModel ? control.aoa_effector : control
                    if panel isa SimulationModel.SolarPanelAngleOfAttackControlModel && index in panel.controlled_panel_links
                        throw(ArgumentError("Checkpointing a geometry-aware interlink on controlled panel link $(index) " *
                            "of spacecraft $(position) requires a PanelHinge to restore its mounting attitude."))
                    end
                end
                (kind=:fixed, q=_interlink_kind_parameter(Tuple(link.q)))
            else
                (kind=:hinged, hinges=[(axis=_interlink_kind_parameter(Tuple(h.axis)),
                    min_angle_rad=_interlink_kind_parameter(h.min_angle_rad),
                    max_angle_rad=_interlink_kind_parameter(h.max_angle_rad),
                    max_rate_rad_s=_interlink_kind_parameter(h.max_rate_rad_s)) for h in link.hinges])
            end
            push!(mounts, (endpoint=(position, number), link=index, root=link.root,
                r=_interlink_kind_parameter(Tuple(link.r)), rotation=rotation))
        end
    end
    return mounts
end

# A ground station's identity: its site and mask, each float with its bits, and
# each antenna as a terminal's identity. Its name is a label, not part of it.
_ground_station_identity(station) = (
    site=map(_interlink_kind_parameter, (lat_deg=station.lat_deg, lon_deg=station.lon_deg, alt_m=station.alt_m,
        min_elevation_deg=station.min_elevation_deg)),
    antennas=Any[_interlink_terminal_identity(antenna) for antenna in station.antennas])

# The scheduling policy's entry: `policy`, `target_idx` and `target_idxs`, the keys
# checkpoints have always been written with. A policy that 0.5 chose by a symbol
# (`SemimajorAxisRatePolicy`, `EccentricityRatePolicy`, `PowerDeliveryPolicy`,
# `DataThroughputPolicy`) is recorded as 0.5 recorded it, with its target index and
# target set, so those checkpoints resume and the entry is unchanged for the same
# policy. Any other policy is recorded by its kind and every one of its fields
# (`_scheduling_policy_identity`), with no target keys.
function _scheduling_policy_entry(policy)
    recorded = SimulationModel.SchedulingPolicies._symbol_policy_fields(policy)
    recorded === nothing || return (policy=recorded[1], target_idx=recorded[2], target_idxs=recorded[3])
    return (policy=_scheduling_policy_identity(policy), target_idx=nothing, target_idxs=nothing)
end

# A policy's identity: its kind, named as link kinds are, and every one of its
# fields, recursively, by value, so that an equal policy rebuilt in a new process
# has the same identity:
# - a float with its bits, a bits composite field by field, as terminal fields
#   are; another number, a string or a symbol by its value; a type by its name;
# - a field that is a policy, or a weighted mix's term, by its own identity, a
#   tuple element by element, an array by its type, its size and its elements;
# - a function by its type and, for a closure, every value it captures (its
#   fields), so that two closures made by one expression with different captures
#   differ. A named function is its type alone: a change to its code is not part
#   of the identity;
# - any other value with fields (a struct, mutable or not, a named tuple, a
#   `Ref`) by its type and its fields.
# What it cannot record by value is refused: a closure that captures a variable
# reassigned after the closure was made (a `Core.Box`); a dictionary or a set; a
# channel, stream, task or lock; a pointer, or a value holding one (a `Regex`),
# whose address differs from one process to the next; a value that is not a
# struct and not one of the above; a field or element that is not assigned; and a
# payload that contains itself.
function _scheduling_policy_identity(policy, path=Base.IdSet{Any}())
    kind = typeof(policy)
    names = fieldnames(kind)
    return (kind=_interlink_kind_identity(kind),
        fields=NamedTuple{names}(map(name -> _policy_field_identity(getfield(policy, name), path), names)))
end
_policy_field_identity(value::SimulationModel.AbstractSchedulingPolicy, path) = _scheduling_policy_identity(value, path)
_policy_field_identity(value::SimulationModel.PolicyTerm, path) = _scheduling_policy_identity(value, path)
_policy_field_identity(value::Tuple, path) = map(element -> _policy_field_identity(element, path), value)
_policy_field_identity(value::Union{Number, AbstractString, Symbol, Type, Module}, path) = _interlink_kind_parameter(value)
_policy_field_identity(value::AbstractArray, path) = _on_path(path, value) do
    all(index -> isassigned(value, index), eachindex(value)) || _unrecordable_field(value, :element)
    (type=_interlink_kind_identity(typeof(value)), size=size(value),
        elements=[_policy_field_identity(element, path) for element in vec(value)])
end
function _policy_field_identity(value::Function, path)
    kind = typeof(value)
    names = fieldnames(kind)
    isempty(names) && return (callable=_interlink_kind_identity(kind),)
    return _on_path(path, value) do
        (callable=_interlink_kind_identity(kind),
            captures=NamedTuple{names}(map(name -> _policy_field_identity(getfield(value, name), path), names)))
    end
end
_policy_field_identity(::Core.Box, path) = throw(ArgumentError(
    "A scheduling policy holds a closure that captures a variable reassigned after the closure was made, " *
    "which a checkpoint cannot record. Capture the value instead (for example in a `let` block)."))
_policy_field_identity(value::Union{AbstractDict, AbstractSet, AbstractChannel, IO, Task, Base.AbstractLock}, path) =
    _unrecordable_field(value)
_policy_field_identity(value::Ptr, path) = _unrecordable_field(value, :pointer)
function _policy_field_identity(value, path)
    kind = typeof(value)
    _holds_pointer(kind) && _unrecordable_field(value, :pointer)
    isbits(value) && return _interlink_kind_parameter(value)
    isstructtype(kind) || _unrecordable_field(value)
    names = fieldnames(kind)
    for name in names
        isdefined(value, name) || _unrecordable_field(value, :field, name)
    end
    encode() = (type=_interlink_kind_identity(kind),
        fields=NamedTuple{names}(map(name -> _policy_field_identity(getfield(value, name), path), names)))
    return ismutable(value) ? _on_path(encode, path, value) : encode()
end

# Whether a value of `kind` holds a pointer in a field, or in a bits field's own
# fields. A field of an abstract type that holds a pointer is found when its
# value is encoded.
_holds_pointer(kind::DataType) =
    any(field -> field <: Ptr || (isbitstype(field) && _holds_pointer(field)), fieldtypes(kind))

# The refusal of a value the identity cannot record: `reason` is `:value` (its
# kind), `:element` or `:field` (`name`, not assigned) or `:pointer`.
function _unrecordable_field(value, reason::Symbol=:value, name=nothing)
    held = "a value of type `$(repr(typeof(value)))`" * (
        reason === :element ? " with an unassigned element" :
        reason === :field ? " whose field `$(name)` is not assigned" :
        reason === :pointer ? ", which holds a pointer whose address differs from one process to the next" : "")
    advice = reason === :pointer ? "Hold the values it is built from (for a `Regex`, its pattern string)" :
        reason === :value ? "Hold numbers, strings, symbols, functions, arrays, tuples, named tuples and structs of these" :
        "Assign it before the run"
    throw(ArgumentError("A scheduling policy holds $(held), which a checkpoint cannot record by value. " *
        "$(advice), or run without checkpoints."))
end

# Encode `value` with `encode` while `value` is on the path from the policy,
# refusing a payload that contains itself.
function _on_path(encode, path, value)
    value in path && throw(ArgumentError(
        "A scheduling policy's fields contain themselves, which a checkpoint cannot record."))
    push!(path, value)
    try
        return encode()
    finally
        delete!(path, value)
    end
end

# A snapshot's policy entry in the current form: a symbol with its targets (all
# checkpoints written before 0.6, and those of the policies they had) is rebuilt as
# the typed policy it stands for and recorded as that policy is, so a checkpoint
# resumes under the typed policy that schedules the same. A symbol that stands for no
# policy is left as written, and the comparison refuses it.
function _canonical_policy_configuration(configuration)
    hasproperty(configuration, :policy) && configuration.policy isa Symbol || return configuration
    policy = SimulationModel.SchedulingPolicies._symbol_scheduling_policy(configuration.policy,
        get(configuration, :target_idx, 1), get(configuration, :target_idxs, nothing))
    policy === nothing && return configuration
    return merge(configuration, _scheduling_policy_entry(policy))
end

_mounted_terminal_identities(spacecraft::SimulationModel.SpacecraftModel) =
    Any[_interlink_terminal_identity(entry.terminal) for entry in SimulationModel.mounted_link_terminals(spacecraft)]
_mounted_terminal_identities(spacecraft) = Any[]

# A terminal's identity is its kind, named as link kinds are, and every one of
# its fields, encoded as type parameters are: a float with its bits, a composite
# bits value (an `SVector`) field by field. Its fields set what its links do (a
# receiver's aperture, efficiency and limit, a transmitter's efficiency and
# limits), and the pointing models to come will read the rest, so a rebuilt
# terminal with any other value refuses the checkpoint.
function _interlink_terminal_identity(terminal)
    kind = typeof(terminal)
    names = fieldnames(kind)
    return (kind=_interlink_kind_identity(kind),
        fields=NamedTuple{names}(map(name -> _interlink_kind_parameter(getfield(terminal, name)), names)))
end

function _interlink_checkpoint_runtime_state(args, is_active)
    args.interlink_model === nothing && return nothing
    count = length(args.dynamics_model.spacecraft)
    is_active isa AbstractVector{Bool} && length(is_active) == count ||
        throw(ArgumentError("Writing an interlink checkpoint requires one activity flag per spacecraft."))
    return (interlinks=(
        version=2,
        configuration=_interlink_checkpoint_configuration(args),
        is_active=Vector{Bool}(is_active),
        # Scheduling now uses the previous matching for receiver-outflow credit.
        # Keep it independently of the optional/recoverable history sidecar.
        active_links=sort!([key for (key, edge) in args.interlink_model.graph.edge_index_lookup
            if args.interlink_model.graph.active[edge]]),
    ),)
end

# A snapshot's configuration in the current form. Version 1 (before link kinds)
# held one laser link type for the run and P, B and range per candidate. Both
# versions' policy entries are rebuilt as typed policies.
function _interlink_snapshot_configuration(snapshot)
    configuration = snapshot.configuration
    snapshot.version === 2 && return _canonical_policy_configuration(configuration)
    hasproperty(configuration, :link_type) && configuration.link_type === :laser || return nothing
    candidates = [(endpoints=candidate.endpoints,
        link=_interlink_link_identity(SimulationModel.LaserMomentumLink(candidate.P, candidate.B, candidate.range)))
        for candidate in configuration.candidates]
    rest = (name => getproperty(configuration, name) for name in keys(configuration)
        if !(name in (:link_type, :vertices, :candidates)))
    return _canonical_policy_configuration((; vertices=configuration.vertices, candidates=candidates, rest...))
end

function _restore_interlink_checkpoint_state!(args, checkpoint, is_active)
    runtime = get(checkpoint, :runtime_state, nothing)
    snapshot = runtime !== nothing && hasproperty(runtime, :interlinks) ? runtime.interlinks : nothing
    if args.interlink_model === nothing
        snapshot === nothing || throw(ArgumentError(
            "An interlink checkpoint cannot resume with interlink_model=nothing."))
        return nothing
    end
    snapshot !== nothing && hasproperty(snapshot, :version) && snapshot.version in (1, 2) ||
        throw(ArgumentError("Interlink checkpoint runtime metadata version 1 or 2 is required; " *
            "older checkpoints without this metadata cannot safely resume an interlinked run."))
    hasproperty(snapshot, :configuration) &&
        isequal(_interlink_snapshot_configuration(snapshot), _interlink_checkpoint_configuration(args)) ||
        throw(ArgumentError("Interlink configuration differs from the checkpoint. Resume with the same " *
            "candidate topology, link kinds and parameters, forbidden pairs, scheduling policy, terminal counts, " *
            "terminal mounting geometry and ground stations."))
    count = length(args.dynamics_model.spacecraft)
    hasproperty(snapshot, :is_active) && snapshot.is_active isa AbstractVector{Bool} &&
        length(snapshot.is_active) == count && is_active isa AbstractVector{Bool} &&
        length(is_active) == count ||
        throw(ArgumentError("Interlink checkpoint activity flags must match the spacecraft count."))
    hasproperty(checkpoint.u, :sc) && length(checkpoint.u.sc) == count ||
        throw(ArgumentError("Interlink checkpoint state must match the spacecraft count."))
    stations = SimulationModel.InterLinkModels.has_station_counters(args.interlink_model) ?
        length(args.interlink_model.ground_stations) : 0
    (hasproperty(checkpoint.u, :gs) ? length(checkpoint.u.gs) : 0) == stations ||
        throw(ArgumentError("Interlink checkpoint state must hold one data counter pair per ground station " *
            "($(stations) in this run)."))
    graph = args.interlink_model.graph
    restored_active = falses(length(graph.active))
    if hasproperty(snapshot, :active_links)
        selected = snapshot.active_links
        selected isa AbstractVector && allunique(selected) &&
            all(key -> haskey(graph.edge_index_lookup, key), selected) ||
            throw(ArgumentError("Checkpoint active interlinks must be distinct candidates in the restored topology."))
        allunique(collect(Iterators.flatten(selected))) || throw(ArgumentError(
            "Checkpoint active interlinks must not share terminals."))
        for key in selected
            restored_active[graph.edge_index_lookup[key]] = true
        end
    elseif SimulationModel.SchedulingPolicies._uses_backlog_outflow(args.scheduling_policy_model.objective)
        @warn "Checkpoint has no previous active interlink matching; receiver-outflow credit starts at zero."
    end
    # Availability is refreshed at the checkpoint's time before scoring; the
    # old matching never grants credit for a link that is now unavailable.
    copyto!(is_active, snapshot.is_active)
    copyto!(graph.active, restored_active)
    return nothing
end

# The interlink states' names. A run's state calls them `interlink_dv` and
# `interlink_delta_sma`; a checkpoint keeps the names they had before 0.6.0,
# `laser_dv` and `laser_delta_sma`, so its integrated state is written as earlier
# versions wrote it and those versions can resume it. A resumed run reads either
# pair.
const _INTERLINK_STATE_NAMES = (:interlink_dv => :laser_dv, :interlink_delta_sma => :laser_delta_sma)

# Which pair of names a state gives the interlink states: `:current`, `:legacy`
# (the names before 0.6.0) or `:none`; `:mixed` for a state that has parts of
# both pairs, or whose spacecraft differ.
function _interlink_state_naming(u)
    namings = map(eachindex(u.sc)) do i
        names = propertynames(u.sc[i])
        current = count(in(names), first.(_INTERLINK_STATE_NAMES))
        legacy = count(in(names), last.(_INTERLINK_STATE_NAMES))
        current == 2 && legacy == 0 ? :current : legacy == 2 && current == 0 ? :legacy :
            current == 0 && legacy == 0 ? :none : :mixed
    end
    isempty(namings) && return :none
    return allequal(namings) ? first(namings) : :mixed
end

# A checkpointed state's interlink states fit the run: both names of one pair in a
# run with interlinks, neither pair in a run without.
function _check_interlink_checkpoint_state(args, checkpoint)
    naming = _interlink_state_naming(checkpoint.u)
    (args.interlink_model === nothing ? naming === :none : naming in (:current, :legacy)) ||
        throw(ArgumentError("Checkpoint interlink state does not match the configured interlink_model; resume with the same interlink setting."))
    return nothing
end

# The state a checkpoint stores: the run's state, with the interlink states under
# their checkpoint names in a run with interlinks.
_checkpoint_state(args, u) = args.interlink_model === nothing ? u : _renamed_spacecraft_state(u, _INTERLINK_STATE_NAMES)

# The state a resumed run starts from: the checkpoint's state, with interlink states
# written under their checkpoint names given back the run's names.
_resumed_checkpoint_state(args, u) = args.interlink_model !== nothing && _interlink_state_naming(u) === :legacy ?
    _renamed_spacecraft_state(u, map(reverse, _INTERLINK_STATE_NAMES)) : u

# A copy of the state `u` with the spacecraft components renamed by `names` (old =>
# new pairs). Every other component, every shape and every value stay (the ground
# stations' counters, `gs`, as they are); the copy's values are checked bit for bit
# against `u`'s.
function _renamed_spacecraft_state(u::ComponentVector, names)
    spacecraft = [_renamed_components(u.sc[i], names) for i in eachindex(u.sc)]
    state = hasproperty(u, :gs) ?
        ComponentVector(sc=spacecraft, gs=[_renamed_components(u.gs[s], ()) for s in eachindex(u.gs)]) :
        ComponentVector(sc=spacecraft)
    old, new = ComponentArrays.getdata(u), ComponentArrays.getdata(state)
    length(old) == length(new) && all(index -> old[index] === new[index], eachindex(old, new)) ||
        error("Renaming the interlink states changed the state's values; this is a SpaceAGORA bug.")
    return state
end
# A state given as named tuples, as unit tests build one.
_renamed_spacecraft_state(u::NamedTuple, names) = merge(u, (sc=[_renamed_components(spacecraft, names) for spacecraft in u.sc],))

function _renamed_components(spacecraft, names)
    renamed(name) = (index = findfirst(pair -> first(pair) === name, names); index === nothing ? name : last(names[index]))
    components = propertynames(spacecraft)
    return NamedTuple{map(renamed, components)}(map(name -> _state_component_copy(getproperty(spacecraft, name)), components))
end
_state_component_copy(value::Number) = value
_state_component_copy(value::AbstractArray) = copy(value)

# The run state the callbacks keep on the spacecraft records: the batteries'
# flags, the faults' activity and held levels, the sensors' streams, latest
# measurements and gyro biases, the hinges' and gimbals' angles, the firing
# state of the mounted thrusters with firing dynamics (the ramp factor `κ` and
# the average thrust `thrust` their last control tick scheduled), and, in a run
# with data stores, the stores' flags. A plain-data snapshot written with every
# checkpoint and restored, after checking it belongs to the same configuration,
# before a resumed run's callbacks initialize.
_record_faults(sc) = sc.root.faults.faults
_record_hinges(sc) = [(link, h) for link in sc.links for h in link.hinges]
_record_gimbals(sc) = [(link, g) for link in sc.links for g in link.gimbals]
_record_firing(sc) = [(link, t) for link in sc.links for t in link.thrusters if t.firing_dynamics]

# Runs whose resume needs the snapshot (their resume was refused before it
# existed, or its firing state is missing from it), and runs with articulated
# actuators whose angles it also carries.
_needs_callback_state(sc) = SimulationModel.battery(sc) !== nothing || !isempty(_record_faults(sc)) ||
    !isempty(SimulationModel.mounted_sensors(sc)) || !isempty(_record_firing(sc)) || SimulationModel.has_data_store(sc)

# A run with data stores records each spacecraft's store, `nothing` for one
# without; a run without records none, so its snapshot is as before data stores.
_records_data_stores(args) = any(SimulationModel.has_data_store, args.dynamics_model.spacecraft)
_has_callback_state(sc) = _needs_callback_state(sc) || !isempty(_record_hinges(sc)) || !isempty(_record_gimbals(sc))

# Version 2 added the thrusters' firing state, and the link of a faulted
# actuator. A version 1 snapshot (SpaceAGORA 0.4 and earlier) still resumes a run
# without firing dynamics; a run with them refuses it, since it cannot restore
# their firing state.
const CALLBACK_CHECKPOINT_STATE_VERSION = 2

# Each record with where it sits, so a record moved to another device of the same
# kind refuses the snapshot: a fault or a thruster with firing dynamics with its
# actuator's place (`_actuator_place`), a sensor, hinge or gimbal with its link's
# position, a gimbal also with its thruster's on that link. A thruster also with
# the parameters of its firing dynamics, which its saved state follows.
_link_position(sc, link) = findfirst(l -> l === link, sc.links)

# Where an actuator is mounted: its position in the mounting order (as
# `resolve_faults!` finds it), its link's position, and its position among the
# actuators on that link, in mounting order. The mounting order alone stays the
# same when an actuator moves to another link with nothing mounted between.
function _actuator_place(sc, actuator)
    entries = SimulationModel.mounted_actuators(sc)
    k = findfirst(e -> e.actuator === actuator, entries)
    k === nothing && return nothing
    link = entries[k].link
    return (k, _link_position(sc, link), count(e -> e.link === link, view(entries, 1:k)))
end

# The configuration a snapshot of `version` was written with. Version 1 had no
# thrusters, so a run with firing dynamics has none to compare (it is refused),
# and placed a faulted actuator by its position in the mounting order alone.
function _callback_checkpoint_configuration(args, version=CALLBACK_CHECKPOINT_STATE_VERSION)
    stores = _records_data_stores(args)
    return [begin
        configuration = (
            battery=(b = SimulationModel.battery(sc); b === nothing ? nothing : (b.capacity_wh, b.recovery_soc)),
            faults=[(nameof(typeof(f.actuator)), (place = _actuator_place(sc, f.actuator); version === 1 ? first(place) : place),
                     f.kind, f.start_s, f.stop_s, f.value) for f in _record_faults(sc)],
            sensors=[(nameof(typeof(e.sensor)), _link_position(sc, e.link), e.sensor.period_s, e.sensor.seed)
                     for e in SimulationModel.mounted_sensors(sc)],
            hinges=[_link_position(sc, link) for (link, _) in _record_hinges(sc)],
            gimbals=[(_link_position(sc, link), findfirst(t -> t === g.thruster, link.thrusters)) for (link, g) in _record_gimbals(sc)],
        )
        configuration = version === 1 ? configuration : merge(configuration, (thrusters=[(_actuator_place(sc, t), t.max_thrust,
            t.cutoff_frequency, t.min_firing_time, t.level_on, t.level_off) for (_, t) in _record_firing(sc)],))
        stores ? merge(configuration, (data_store=(s = SimulationModel.data_store(sc);
            s === nothing ? nothing : (s.capacity_bits, s.initial_bits)),)) : configuration
    end for sc in args.dynamics_model.spacecraft]
end

function _callback_snapshot_configuration(args, version)
    version === 1 && any(sc -> !isempty(_record_firing(sc)), args.dynamics_model.spacecraft) && throw(ArgumentError(
        "This run has thrusters with firing dynamics, whose firing state this checkpoint does not hold " *
        "(it was written by SpaceAGORA 0.4 or earlier). Run from the start."))
    return _callback_checkpoint_configuration(args, version)
end

function _callback_checkpoint_state(args)
    spacecraft = args.dynamics_model.spacecraft
    any(_has_callback_state, spacecraft) || return nothing
    return (
        version=CALLBACK_CHECKPOINT_STATE_VERSION,
        configuration=_callback_checkpoint_configuration(args),
        spacecraft=[merge((
            battery=(b = SimulationModel.battery(sc); b === nothing ? nothing : (full=b.full, empty=b.empty, shed_tier=b.shed_tier)),
            faults=[(active=f.active, held_level=f.held_level) for f in _record_faults(sc)],
            sensors=[(rng=(s = e.sensor; s.rng === nothing ? nothing : copy(s.rng)), latest=deepcopy(e.sensor.latest),
                      latest_time_s=e.sensor.latest_time_s, bias=hasfield(typeof(e.sensor), :bias) ? e.sensor.bias : nothing)
                     for e in SimulationModel.mounted_sensors(sc)],
            hinges=[(angle_rad=h.angle_rad, last_update_s=h.last_update_s) for (_, h) in _record_hinges(sc)],
            gimbals=[(angle_rad=g.angle_rad, last_update_s=g.last_update_s, previous_angle_rad=g.previous_angle_rad)
                     for (_, g) in _record_gimbals(sc)],
            thrusters=[(κ=t.κ, thrust=t.thrust) for (_, t) in _record_firing(sc)],
        ), _data_store_checkpoint_state(args, sc)) for sc in spacecraft],
    )
end

_data_store_checkpoint_state(args, sc) = _records_data_stores(args) ?
    (data_store=(s = SimulationModel.data_store(sc); s === nothing ? nothing : (full=s.full, empty=s.empty)),) : (;)

function _restore_callback_checkpoint_state!(args, checkpoint, p)
    spacecraft = args.dynamics_model.spacecraft
    runtime = get(checkpoint, :runtime_state, nothing)
    snapshot = runtime !== nothing && hasproperty(runtime, :callbacks) ? runtime.callbacks : nothing
    # A snapshot is always checked, so a run without the devices it holds
    # refuses it rather than dropping them.
    if snapshot === nothing
        any(_has_callback_state, spacecraft) || return nothing
        any(_needs_callback_state, spacecraft) && throw(ArgumentError(
            "This run has batteries, actuator faults, sensors, thrusters with firing dynamics or data stores, whose run state " *
            "this checkpoint does not hold (it was written by an earlier SpaceAGORA or by another configuration). " *
            "Run from the start."))
        @warn "The checkpoint holds no hinge or gimbal angles; the resumed run starts them from their current values."
        return nothing
    end
    version = hasproperty(snapshot, :version) ? snapshot.version : nothing
    version === 1 || version === CALLBACK_CHECKPOINT_STATE_VERSION ||
        throw(ArgumentError("Callback checkpoint state version 1 or $(CALLBACK_CHECKPOINT_STATE_VERSION) is required to resume this run."))
    hasproperty(snapshot, :configuration) && isequal(snapshot.configuration, _callback_snapshot_configuration(args, version)) ||
        throw(ArgumentError("The batteries, actuator faults, sensors, hinges, gimbals, thrusters with firing dynamics " *
            "or data stores differ from the checkpoint's. Resume with the same configuration, or run from the start."))
    length(snapshot.spacecraft) == length(spacecraft) ||
        throw(ArgumentError("The callback checkpoint state must match the spacecraft count."))
    # Checked above: from here on the restore cannot fail part-way.
    for (sc, saved) in zip(spacecraft, snapshot.spacecraft)
        b = SimulationModel.battery(sc)
        if b !== nothing
            b.full, b.empty, b.shed_tier = saved.battery.full, saved.battery.empty, saved.battery.shed_tier
        end
        store = SimulationModel.data_store(sc)
        if store !== nothing
            store.full, store.empty = saved.data_store.full, saved.data_store.empty
        end
        for (f, s) in zip(_record_faults(sc), saved.faults)
            f.active, f.held_level = s.active, s.held_level
        end
        for (e, s) in zip(SimulationModel.mounted_sensors(sc), saved.sensors)
            sensor = e.sensor
            sensor.rng = s.rng === nothing ? nothing : copy(s.rng)
            sensor.latest = deepcopy(s.latest)
            sensor.latest_time_s = s.latest_time_s
            s.bias === nothing || (sensor.bias = s.bias)
        end
        for ((link, h), s) in zip(_record_hinges(sc), saved.hinges)
            h.angle_rad, h.last_update_s = s.angle_rad, s.last_update_s
            SimulationModel.ActuatorMounting._set_hinge_attitude!(link, h, h.angle_rad)
        end
        for ((link, g), s) in zip(_record_gimbals(sc), saved.gimbals)
            g.angle_rad, g.last_update_s, g.previous_angle_rad = s.angle_rad, s.last_update_s, s.previous_angle_rad
            SimulationModel.ActuatorMounting._update_gimbal_rotation!(link, g.thruster)
        end
        # Version 1 holds no thrusters: the run has none with firing dynamics (checked).
        version === 1 && continue
        for ((_, t), s) in zip(_record_firing(sc), saved.thrusters)
            t.κ, t.thrust = s.κ, s.thrust
        end
    end
    p.shared_buffers.callback_resume[] = true
    return nothing
end

# The run state the models keep themselves (`checkpoint_state`): the guidance,
# navigation and control models in the configuration's order, then each
# spacecraft's sensors in mounting order. Written next to the callbacks' snapshot
# when any of them keeps some, with each model's kind at its position; restored,
# after checking those match the run's, once the run-start initializers ran.
function _checkpoint_models(args)
    return (
        navigation=collect(Any, args.navigation_model.navigation_effectors),
        guidance=collect(Any, args.guidance_model.guidance_effectors),
        control=collect(Any, args.control_model.control_effectors),
        sensors=[Any[e.sensor for e in SimulationModel.mounted_sensors(sc)] for sc in args.dynamics_model.spacecraft],
    )
end

_each_model(groups) = Iterators.flatten((groups.navigation, groups.guidance, groups.control, Iterators.flatten(groups.sensors)))
_model_kinds(group) = [nameof(typeof(m)) for m in group]
_model_checkpoint_configuration(models) = (navigation=_model_kinds(models.navigation),
    guidance=_model_kinds(models.guidance), control=_model_kinds(models.control),
    sensors=[_model_kinds(group) for group in models.sensors])

_missing_restore(model) = ArgumentError("$(nameof(typeof(model))) returns run state from checkpoint_state but has no " *
    "restore_checkpoint_state! method for it; implement SpaceAGORA.restore_checkpoint_state!(::$(nameof(typeof(model))), data), " *
    "or return nothing from checkpoint_state.")

# A checkpointed run's models that keep state can put it back.
function _check_model_checkpoint_hooks(args)
    for model in _each_model(_checkpoint_models(args))
        data = SimulationLifecycle.checkpoint_state(model)
        data === nothing || applicable(SimulationLifecycle.restore_checkpoint_state!, model, data) ||
            throw(_missing_restore(model))
    end
    return nothing
end

function _model_checkpoint_state(args)
    models = _checkpoint_models(args)
    saved = map(group -> Any[SimulationLifecycle.checkpoint_state(m) for m in group],
        (navigation=models.navigation, guidance=models.guidance, control=models.control))
    sensors = [Any[SimulationLifecycle.checkpoint_state(m) for m in group] for group in models.sensors]
    all(data -> data === nothing, _each_model((saved..., sensors=sensors))) && return nothing
    return (version=1, configuration=_model_checkpoint_configuration(models), saved..., sensors=sensors)
end

# Checked when the checkpoint is read, before anything is restored; `nothing`
# when the checkpoint holds no model state.
function _check_model_checkpoint_state(args, checkpoint)
    runtime = get(checkpoint, :runtime_state, nothing)
    snapshot = runtime !== nothing && hasproperty(runtime, :models) ? runtime.models : nothing
    snapshot === nothing && return nothing
    hasproperty(snapshot, :version) && snapshot.version === 1 ||
        throw(ArgumentError("Model checkpoint state version 1 is required to resume this run."))
    models = _checkpoint_models(args)
    hasproperty(snapshot, :configuration) && isequal(snapshot.configuration, _model_checkpoint_configuration(models)) ||
        throw(ArgumentError("The guidance, navigation or control models or the sensors differ from the checkpoint's. " *
            "Resume with the same configuration, or run from the start."))
    for (model, data) in zip(_each_model(models), _each_model(snapshot))
        data === nothing || applicable(SimulationLifecycle.restore_checkpoint_state!, model, data) ||
            throw(_missing_restore(model))
    end
    return snapshot
end

# After the run-start initializers of a resumed run. Before anything is
# restored, a model that keeps run state now and has none in the file is
# refused: no snapshot at all, or `nothing` saved at its position (another
# model wrote the snapshot, and this one kept no state then).
function _restore_model_checkpoint_state!(args, snapshot)
    models = collect(_each_model(_checkpoint_models(args)))
    saved = snapshot === nothing ? fill(nothing, length(models)) : collect(_each_model(snapshot))
    kept = unique(nameof(typeof(m)) for (m, data) in zip(models, saved)
                  if data === nothing && SimulationLifecycle.checkpoint_state(m) !== nothing)
    isempty(kept) || throw(ArgumentError("The run state of $(join(kept, ", ")) (checkpoint_state) is not in this checkpoint " *
        "(it was written by an earlier SpaceAGORA or by another configuration, or the model kept none then). Run from the start."))
    for (model, data) in zip(models, saved)
        data === nothing || SimulationLifecycle.restore_checkpoint_state!(model, data)
    end
    return nothing
end

# The spacecraft each guidance, navigation and control model acts on, as its
# `spacecraft_scope` reads once the run-start initializers ran (the RPO planner
# and controller set their positions there, from the spacecraft ids). A
# checkpointed or resumed run takes this baseline after its initializers and
# writes it with every checkpoint (`gnc_scopes`). A resumed run compares its own
# baseline with the checkpoint's before any model's run state is restored, so
# saved state never goes back into a model that acts on other spacecraft, such as
# a momentum manager's held command into a manager of another spacecraft, or a
# bound copy declared on another spacecraft. Each group's scopes are compared in
# configuration order (the configuration's models, then each spacecraft's bound
# copies, as `_flatten_spacecraft_gnc` orders them). The models' types are
# compared with their run state (`models`), as before. A run without these
# models writes the entry too, with empty groups, so a checkpoint without it was
# written by SpaceAGORA 0.5.1 or earlier; a run with any of them refuses such a
# checkpoint when it is read (`_gnc_scope_checkpoint_entry`).
#
# A scope is plain data built afresh, so a selection changed in place later
# cannot change a baseline already taken:
# - `AllSpacecraft()`: `(scope=:all,)`;
# - `BoundSpacecraft`: `(scope=:bound, positions)`, the run positions its ticks
#   visit: those it names within the run, ascending and each once. A position
#   outside the run, which an initializer can set after the pre-run check, is
#   never ticked and not recorded, so checkpointing does not decide whether such
#   a run starts;
# - `SpacecraftSelection` (not yet taken by these models in a run): the same form,
#   with the run positions `resolve_spacecraft_scope` resolves it to;
# - `CoupledSpacecraft`: `(scope=:coupled, primary, partners)`, as declared.
# Before the run, `check_spacecraft_scopes` refuses any other scope.
_gnc_scope_identity(::GNCRuntime.AllSpacecraft, spacecraft, owner) = (scope=:all,)
_gnc_scope_identity(scope::GNCRuntime.BoundSpacecraft, spacecraft, owner) =
    (scope=:bound, positions=sort!(unique!(Int[i for i in scope.indices if 1 <= i <= length(spacecraft)])))
_gnc_scope_identity(scope::SimulationModel.SpacecraftSelection, spacecraft, owner) =
    (scope=:bound, positions=SimulationModel.resolve_spacecraft_scope(scope, spacecraft; owner=owner))
_gnc_scope_identity(scope::GNCRuntime.CoupledSpacecraft, spacecraft, owner) =
    (scope=:coupled, primary=scope.primary, partners=collect(Int, scope.partners))

const _GNC_SCOPE_GROUPS = (navigation="navigation", guidance="guidance", control="control")

_gnc_groups(args) = (navigation=args.navigation_model.navigation_effectors,
    guidance=args.guidance_model.guidance_effectors, control=args.control_model.control_effectors)

_gnc_scope_owner(group::Symbol, j::Int, model) = "$(_GNC_SCOPE_GROUPS[group]) model $j ($(nameof(typeof(model))))"

function _gnc_scope_baseline(args)
    spacecraft = args.dynamics_model.spacecraft
    groups = _gnc_groups(args)
    scopes(group) = Any[_gnc_scope_identity(GNCRuntime.spacecraft_scope(m), spacecraft, _gnc_scope_owner(group, j, m))
        for (j, m) in enumerate(getproperty(groups, group))]
    return (version=1, navigation=scopes(:navigation), guidance=scopes(:guidance), control=scopes(:control))
end

# The checkpoint's entry, checked for its form when the checkpoint is read, before
# anything is restored. A checkpoint without one cannot show which spacecraft its
# models acted on, so a run with any guidance, navigation or control model, of
# any scope, refuses it; a run without them resumes from it (`nothing`), having
# nothing to compare.
function _gnc_scope_checkpoint_entry(args, checkpoint)
    runtime = get(checkpoint, :runtime_state, nothing)
    entry = runtime !== nothing && hasproperty(runtime, :gnc_scopes) ? runtime.gnc_scopes : nothing
    if entry === nothing
        models = [_gnc_scope_owner(group, j, m) for (group, group_models) in pairs(_gnc_groups(args))
            for (j, m) in enumerate(group_models)]
        isempty(models) && return nothing
        named = length(models) <= 3 ? join(models, ", ") : join(models[1:3], ", ") * " and $(length(models) - 3) more"
        throw(ArgumentError(
            "The checkpoint records no spacecraft scopes for guidance, navigation and control models: it was written " *
            "by SpaceAGORA 0.5.1 or earlier, so this run's $named cannot be checked against it. Resume it with " *
            "SpaceAGORA 0.5.1, or run from the start."))
    end
    hasproperty(entry, :version) && entry.version === 1 && all(group -> hasproperty(entry, group) &&
        getproperty(entry, group) isa AbstractVector, keys(_GNC_SCOPE_GROUPS)) ||
        throw(ArgumentError("GNC scope checkpoint state version 1 is required to resume this run."))
    for group in keys(_GNC_SCOPE_GROUPS), (j, scope) in enumerate(getproperty(entry, group))
        _gnc_scope_form(scope) || throw(ArgumentError(
            "The checkpoint's GNC scope entry is damaged: it gives $(_GNC_SCOPE_GROUPS[group]) model $j the scope " *
            "$(repr(scope)), which is not a version 1 scope. Run from the start."))
    end
    return entry
end

# Whether `scope` has one of the forms `_gnc_scope_identity` writes.
_gnc_scope_form(scope) = scope isa NamedTuple && hasproperty(scope, :scope) && (scope.scope === :all ||
    scope.scope === :bound && hasproperty(scope, :positions) && scope.positions isa AbstractVector{Int} ||
    scope.scope === :coupled && hasproperty(scope, :primary) && scope.primary isa Int &&
        hasproperty(scope, :partners) && scope.partners isa AbstractVector{Int})

function _describe_gnc_scope(scope)
    scope.scope === :all && return "every spacecraft"
    scope.scope === :bound && return isempty(scope.positions) ? "no spacecraft" :
        "the spacecraft at run position$(length(scope.positions) == 1 ? "" : "s") $(join(scope.positions, ", "))"
    scope.scope === :coupled && return "the spacecraft at run position $(scope.primary), coupled with run " *
        "position$(length(scope.partners) == 1 ? "" : "s") $(join(scope.partners, ", "))"
    return repr(scope)
end

# `entry` is `nothing` only for a run without guidance, navigation or control
# models resumed from a checkpoint written before the entry, which has nothing
# to compare (`_gnc_scope_checkpoint_entry` refuses one with models).
function _check_gnc_scope_baseline(args, entry, baseline)
    entry === nothing && return nothing
    for (group, models) in pairs(_gnc_groups(args))
        saved, current = getproperty(entry, group), getproperty(baseline, group)
        isequal(saved, current) && continue
        name = _GNC_SCOPE_GROUPS[group]
        length(saved) == length(current) || throw(ArgumentError(
            "The GNC models' spacecraft scopes differ from the checkpoint's: it has $(length(saved)) $name " *
            "model$(length(saved) == 1 ? "" : "s") and this run $(length(current)). Resume with the same models, " *
            "or run from the start."))
        j = findfirst(k -> !isequal(saved[k], current[k]), eachindex(current))
        owner = _gnc_scope_owner(group, j, models[j])
        throw(ArgumentError(
            "The GNC models' spacecraft scopes differ from the checkpoint's: $owner acts on " *
            "$(_describe_gnc_scope(current[j])) here and on $(_describe_gnc_scope(saved[j])) in the " *
            "checkpoint. Resume with the scopes the checkpoint was written with, or run from the start."))
    end
    return nothing
end

# Every checkpoint a run writes records, for each of the run's periodic
# navigation, guidance and control clocks (in the order of their callbacks), whether
# it had not yet taken the tick due at the checkpoint's time, up to the clock's
# rounding tolerance (`_same_run_time`, on both sides of that time). That is what
# each clock did, not a rule about the time: at a segment boundary the ticks due
# there ran before the checkpoint, and at the run's end they did not. Where a solve
# stopped early, at an event's `terminate!` (impact, touchdown, orbit or entry
# count) or a discrete callback's, a tick exactly then ran; one an ulp or so off ran
# only if its callback came after the termination (the solver applies the step's
# discrete callbacks after its events, in order), and none ran when the stop came
# just before the run's end, whose ticks the run leaves.
function _periodic_tick_checkpoint_state(p, t::Float64)
    pending = NamedTuple{(:clock, :period_s), Tuple{Int, Float64}}[]
    for (position, clock) in pairs(p.shared_buffers.periodic_clocks)
        next_tick = clock.origin[] + clock.index[] * clock.period_s
        SimulationModel.SimulationCallbacks._same_run_time(next_tick, t) &&
            push!(pending, (clock=position, period_s=clock.period_s))
    end
    return (version=2, time_s=t, pending=pending)
end

# The ticks a run resumed from `checkpoint` takes at its start, at `at_s` (the
# checkpoint's time; NaN when there are none), up to `_same_run_time`: those of the
# clocks the checkpoint lists, as (position, period) pairs, matched to the resumed
# run's clocks, or, when it lists none, every clock's tick from `from_s` on.
# Checkpoints written before the list (SpaceAGORA 0.4.0 and earlier) are read as
# that version left them: it compared tick times exactly, so it took a tick just
# before its checkpoint's time at that tick's own time and left one just after it
# (3 * 0.1 after a boundary at 0.3) to the next segment, where the solver could not
# reach it.
# - A version 1 entry (0.4.0, at the run's end): the ticks at the end and just after
#   it are pending.
# - No entry (at a segment boundary or an early stop): the ticks at that time ran
#   before the checkpoint, and those just after it are pending.
function _pending_periodic_ticks(checkpoint)
    t = Float64(checkpoint.t)
    runtime = get(checkpoint, :runtime_state, nothing)
    entry = runtime !== nothing && hasproperty(runtime, :periodic_ticks) ? runtime.periodic_ticks : nothing
    entry === nothing && return (at_s=t, from_s=nextfloat(t), clocks=nothing)
    version = hasproperty(entry, :version) ? entry.version : nothing
    if version === 2 && hasproperty(entry, :time_s) && hasproperty(entry, :pending) &&
            entry.time_s == t && entry.pending isa AbstractVector &&
            all(c -> hasproperty(c, :clock) && hasproperty(c, :period_s) && c.clock isa Int &&
                c.period_s isa Float64, entry.pending)
        return (at_s=t, from_s=-Inf, clocks=Tuple{Int, Float64}[(c.clock, c.period_s) for c in entry.pending])
    elseif version === 1 && hasproperty(entry, :pending_s) && entry.pending_s == t
        return (at_s=t, from_s=t, clocks=nothing)
    end
    throw(ArgumentError("Periodic tick checkpoint state version 1 or 2, at the checkpoint's time, is required to resume this run."))
end

# The results rows kept beside the checkpoint, so that a resumed run's results
# hold the whole run and not only the rows after its resume. A run that records
# rows (it writes or returns results, or runs the gravity backbone) appends the
# rows saved since the last checkpoint at each checkpoint
# (`IOSerialization._append_checkpoint_rows!`). The stream is only ever appended
# to: a record can drop rows from the end of the ones before it, but earlier
# records stay as they are, so an earlier checkpoint still finds its own rows.
# The checkpoint's `results_rows` entry records how much of that stream is its own
# (`bytes`, `rows`), the stream's `token`, which ties the two together, and when
# the rows start (`start_s`): 0 for a run started from its initial conditions, the
# resume time for a resumed run that could not include the rows before it.
#
# A run that streams its results (`_ResultsStream`) keeps in the stream only the
# rows still in memory: the last saved row at each checkpoint. Its entry is
# version 2, with the results streams' lengths beside.
mutable struct _ResultsStream
    spacecraft::String  # the Arrow IPC streams, beside the checkpoint
    run::String
    bytes::Vector{Int}  # the spacecraft and run streams, then their CSVs
    released::Int  # rows of the run that have left memory, written or thinned away
    written::Int  # saved times written
    first_s::Float64  # the first written time, NaN before one
    fresh::Bool  # the next append starts new files
    schemas::Vector{Any}  # the spacecraft and run tables' column names and types, once written
end

mutable struct _ResultsRowsLog
    path::String
    header::NamedTuple
    token::String
    recording::Bool
    start_s::Float64
    logged::Int  # rows of the run's accumulators already in the stream
    bytes::Int
    fresh::Bool  # the next append starts a new stream
    drop::Int  # rows the next record drops from the end of the stream's rows
    retake::Bool  # the restored last row is saved again by the first segment
    stream::Union{Nothing, _ResultsStream}  # the results streams of a run that streams them
    keep_earlier::Bool  # preserve the previous run only when its checkpoint is replaced
    preserve_paths::NTuple{3, String}  # rows and both long-result streams
end

# Whether a run streams its results: `results_layout=:long`, results written with
# the bundle, and checkpoints (the caller's condition for a rows log at all).
function _new_results_stream(args, recording::Bool)::Union{Nothing, _ResultsStream}
    settings = args.simulation_settings
    recording && settings.results_layout === :long && settings.results && _typed_save_bundle_enabled() ||
        return nothing
    paths = _checkpoint_paths(args)
    return _ResultsStream(paths.stream_spacecraft, paths.stream_run, zeros(Int, 4), 0, 0, NaN, true, Any[nothing, nothing])
end

# What the rows were saved with: a resumed run can include them only if it saves
# the same fields for the same number of spacecraft.
_results_rows_header(save_fields, spacecraft_count::Int) = (kind=:spaceagora_results_rows, version=1,
    spacecraft_count=spacecraft_count, fields=[(field.name, field.per_satellite, field.column_prefix) for field in save_fields])

# A token for a new stream, from the clock and the process: it ties the stream to
# the checkpoints of the run that started it, and draws on no random stream.
_results_rows_token(path::String) = string(hash((time_ns(), getpid(), path)); base=16)

# The checkpoint's `results_rows` entry, or `nothing` for a checkpoint written
# before SpaceAGORA kept the rows.
function _results_rows_checkpoint_state(checkpoint)
    runtime = get(checkpoint, :runtime_state, nothing)
    return runtime !== nothing && hasproperty(runtime, :results_rows) ? runtime.results_rows : nothing
end

# The rows up to the checkpoint at `t`, or `nothing` and the reason they cannot be used.
function _checkpoint_results_rows(entry, path::String, header, t::Float64)
    entry === nothing &&
        return nothing, "the checkpoint holds no record of them (it was written before SpaceAGORA kept them)"
    hasproperty(entry, :version) && entry.version in (1, 2) ||
        return nothing, "the checkpoint's record of them has an unknown version"
    entry.recording || return nothing, "the run that wrote the checkpoint did not save results rows"
    rows = SimulationModel.IOSerialization._read_checkpoint_rows(path, entry.bytes)
    rows === nothing && return nothing, "$(path) is missing, shorter than the checkpoint records, or unreadable"
    hasproperty(entry, :token) && hasproperty(rows.header, :token) && rows.header.token == entry.token ||
        return nothing, "$(path) was written by another run than the checkpoint"
    Base.structdiff(rows.header, NamedTuple{(:token,)}) == header ||
        return nothing, "they were saved with other save fields or another number of spacecraft"
    length(rows.times) == entry.rows && !isempty(rows.times) && rows.times[end] == t ||
        return nothing, "they do not reach the checkpoint's time"
    return rows, nothing
end

# Start the run's rows: a resumed run that records rows first takes back the rows
# up to its checkpoint into `times_acc` and `data_acc`, and continues their stream
# after cutting away anything written after the checkpoint. When it cannot, it
# warns: its rows then start at the resume, and its results are written beside the
# earlier ones (`_results_destination`). A resume that solves nothing and returns
# no table needs no rows and reads none.
#
# `retake_last_row`: the checkpoint's run ended before the ticks due at its time,
# so its last row, at the checkpoint's time, does not hold them. The resumed run
# takes those ticks first and saves that row again, as a run never interrupted
# saves it after them. The restored row stays until that row arrives
# (`_retake_restored_row!`), so a solve that fails first keeps it.
function _start_results_rows!(args, recording::Bool, save_fields, resumed::Bool, t_start::Float64, entry,
                              times_acc::Vector{Float64}, data_acc::Vector; solves::Bool=true,
                              returns_table::Bool=true, retake_last_row::Bool=false)::_ResultsRowsLog
    paths = _checkpoint_paths(args)
    path = paths.rows
    header = _results_rows_header(save_fields, length(args.dynamics_model.spacecraft))
    log = _ResultsRowsLog(path, header, _results_rows_token(path), recording, resumed ? t_start : 0.0, 0, 0, true,
        0, false, _new_results_stream(args, recording), resumed && solves && !recording,
        (path, paths.stream_spacecraft, paths.stream_run))
    # A resumed run that saves no rows would remove the earlier rows file at its
    # first checkpoint; it keeps it instead.
    # The checkpoint's run writes its results where its rows start: its results
    # directory, or `resumed_at_<start>s` for a run that was itself a resume
    # without the earlier rows. A resume that solves nothing takes the rows back
    # even without returning a table when that run was killed while writing its
    # streamed files: it finishes them (`_save_streamed_results!`).
    streamed = log.stream !== nothing && entry !== nothing && hasproperty(entry, :version) && entry.version === 2
    own_args = streamed ? _results_destination(args, (start_s=entry.start_s,)) : args
    unfinished = !solves && streamed && _streamed_results_unfinished(own_args)
    (resumed && recording && (solves || returns_table || unfinished)) || return log
    rows, reason = _checkpoint_results_rows(entry, path, header, t_start)
    if rows !== nothing && entry.version === 2
        # A resume that solves nothing cuts the streams back only to finish them.
        reason = _resume_results_stream!(log.stream, entry, own_args; cut=solves || unfinished)
        reason === nothing || (rows = nothing)
    end
    if rows === nothing
        # Until a replacement checkpoint is ready, these files still belong to
        # the old checkpoint. An interrupted first segment must leave it usable.
        log.keep_earlier = solves
        consequence = solves ?
            "They start at the resume, and any results files are written to " *
            "$(_results_destination(args, log).simulation_settings.results_directory), so the earlier ones are kept." :
            "This resume solves nothing, so it returns an empty table and writes no files."
        @warn "The resumed run's results cannot include the rows before its checkpoint at t = $(t_start) s: " *
              "$(reason). $(consequence)"
        return log
    end
    append!(times_acc, rows.times)
    append!(data_acc, rows.data)
    SimulationModel.IOSerialization._truncate_checkpoint_rows!(path, entry.bytes)
    if log.stream !== nothing && entry.version === 2
        stream = log.stream
        stream.bytes = collect(Int, entry.stream.bytes)
        stream.released = entry.stream.released
        stream.written = entry.stream.written
        stream.first_s = entry.stream.first_s
        stream.schemas = collect(Any, entry.stream.schemas)
        stream.fresh = false
    end
    log.token = rows.header.token
    log.start_s = entry.start_s
    log.logged = length(times_acc)
    log.bytes = entry.bytes
    log.fresh = false
    log.keep_earlier = false
    log.retake = retake_last_row
    return log
end

# A resumed run that does not continue the rows file beside its checkpoint (it
# cannot include the rows there, or it saves none) starts a new one, or removes
# it, at its first checkpoint. The earlier file can hold the only copy of the
# earlier run's rows (a run that was killed writes no results files), so it is
# renamed instead: `<rows file>.before_resume_at_<t>s.bin`, numbered when that
# name is taken. Returns the new name, or `nothing` without a file.
function _keep_earlier_rows!(path::String, t_start::Float64)
    isfile(path) || return nothing
    stem, extension = splitext(path)
    stem *= ".before_resume_at_$(t_start)s"
    kept = stem * extension
    n = 1
    while ispath(kept)
        n += 1
        kept = stem * "_$(n)" * extension
    end
    mv(path, kept)
    return kept
end

# Move the old run's files aside immediately before the first replacement
# checkpoint writes its streams or rows, never at resume startup. The stream
# paths matter even when the resumed run itself uses a different output layout.
function _preserve_earlier_results!(log::_ResultsRowsLog)
    log.keep_earlier || return nothing
    kept = [_keep_earlier_rows!(path, log.start_s) for path in log.preserve_paths]
    log.keep_earlier = false
    if log.recording
        kept[1] === nothing || @warn "The rows file beside the checkpoint is kept as $(kept[1])."
        streams = filter(!isnothing, kept[2:3])
        isempty(streams) || @warn "The results streams beside it are kept as $(join(streams, " and "))."
    end
    return nothing
end

# Before a segment's rows join the run's: when the restored last row is to be
# saved again (`retake_last_row`), drop it as the segment's own row at that time
# arrives, and drop it from the stream with the next record.
function _retake_restored_row!(log, times_acc::Vector{Float64}, data_acc::Vector, saved_times)
    (log === nothing || !log.retake) && return nothing
    log.retake = false
    !isempty(saved_times) && !isempty(times_acc) && first(saved_times) == times_acc[end] || return nothing
    pop!(times_acc)
    pop!(data_acc)
    log.logged -= 1
    log.drop += 1
    return nothing
end

# The results streams a version-2 entry recorded, taken back for a resumed run
# that streams too, after checking that each reaches the checkpoint's length and
# is that run's. `args` are those the checkpoint's run wrote its results with, so
# its CSVs are the ones checked and cut. A resume that solves cuts each back to
# that length (`cut`); one that solves nothing leaves the files as the run left
# them, unless it finishes them. Returns why they cannot be taken back, or `nothing`.
function _resume_results_stream!(stream::Union{Nothing, _ResultsStream}, entry, args; cut::Bool=true)
    stream === nothing && return "they were written to a long-layout run's streamed results files, and this run " *
        "does not stream its results (it writes another layout, no results, or no results bundle)"
    hasproperty(entry, :stream) || return "the checkpoint's record of the streamed results is incomplete"
    bytes = entry.stream.bytes
    prefix = _results_bundle_prefix(args)
    paths = (stream.spacecraft, stream.run, prefix * "_spacecraft.csv", prefix * "_run.csv")
    hasproperty(entry.stream, :schemas) || return "the checkpoint's record of the streamed results is incomplete"
    for (path, length_bytes) in zip(paths, bytes)
        length_bytes == 0 || (isfile(path) && filesize(path) >= length_bytes) ||
            return "$(path) is missing or shorter than the checkpoint records"
    end
    # The tokens before anything is cut: a stream another run wrote is left as it is.
    for (path, length_bytes) in zip(paths[1:2], bytes[1:2])
        length_bytes == 0 || _long_stream_token(path) == entry.token ||
            return "$(path) was written by another run than the checkpoint"
    end
    for (path, length_bytes) in zip(paths, bytes)
        cut || break
        length_bytes == 0 ? rm(path; force=true) :
            SimulationModel.IOSerialization._truncate_checkpoint_rows!(path, length_bytes)
    end
    return nothing
end

# Whether the run's saved row `index` is written under the results thinning,
# which counts rows over the whole run as `_thin_results_segment` does: every
# row by default, only the run's last row with SPACEAGORA_RESULTS_FINAL_ONLY,
# every k-th and the last with SPACEAGORA_RESULTS_THIN_STRIDE=k.
function _results_row_written(index::Int, last::Bool)::Bool
    _results_final_only() && return last
    stride = _results_thin_stride()
    return last || stride <= 1 || (index - 1) % stride == 0
end

# A run that streams its results writes the rows it no longer keeps: at a
# checkpoint every row but the last, which the checkpoint keeps beside it (a
# resumed run may save that row again, see `_retake_restored_row!`); at the end
# of the run, or before the results of a failed solve, every row (`final`). The
# rows leave memory, and the rows log then holds only the rows still kept.
# `args` are those the run writes its results with.
function _stream_results_rows!(log, times_acc::Vector{Float64}, data_acc::Vector, save_fields, args; final::Bool=false)
    (log === nothing || log.stream === nothing) && return nothing
    stream = log.stream
    release = final ? length(times_acc) : length(times_acc) - 1
    keep = Int[i for i in 1:max(release, 0) if _results_row_written(stream.released + i, final && i == release)]
    if !isempty(keep) || (final && stream.written == 0)
        tables = _build_long_results_tables(times_acc[keep], data_acc[keep], save_fields, args)
        _append_results_stream!(stream, tables, args, log.token)
        isempty(keep) || isnan(stream.first_s) && (stream.first_s = times_acc[first(keep)])
        stream.written += length(keep)
        final || _write_long_manifest!(args, Dict{String, Any}(), stream.written,
            isnan(stream.first_s) ? nothing : stream.first_s, nothing;
            extra=("streaming" => true, "time_reached_s" => times_acc[release]))
    end
    release > 0 || return nothing
    stream.released += release
    deleteat!(times_acc, 1:release)
    deleteat!(data_acc, 1:release)
    log.drop += log.logged
    log.logged = 0
    return nothing
end

# Append long tables to the run's streams and CSVs; the first append starts new
# files in place of any an earlier run left.
function _append_results_stream!(stream::_ResultsStream, tables, args, token::String)
    prefix = _results_bundle_prefix(args)
    csv(part) = args.simulation_settings.save_csv ? prefix * "_" * part * ".csv" : nothing
    if stream.fresh
        SimulationModel.IOOutputs._invalidate_mission_history!(args)
        foreach(path -> rm(path; force=true), (stream.spacecraft, stream.run, prefix * "_spacecraft.csv", prefix * "_run.csv"))
        stream.fresh = false
    end
    spacecraft = _conformed_to_stream!(stream, 1, tables.spacecraft, "spacecraft")
    run = tables.run === nothing ? nothing : _conformed_to_stream!(stream, 2, tables.run, "run")
    stream.bytes[1], stream.bytes[3] = _append_long_stream!(stream.spacecraft, csv("spacecraft"), spacecraft, token)
    if run !== nothing
        stream.bytes[2], stream.bytes[4] = _append_long_stream!(stream.run, csv("run"), run, token)
    end
    return nothing
end

# Rows to append to a stream must have its columns, in its order and of its
# types, since an Arrow stream has one schema. The first rows set it. A column
# whose stream type also holds a later batch's (`Union{Missing, Float64}` for
# `Float64`) is converted; any other change is refused, naming the column.
function _conformed_to_stream!(stream::_ResultsStream, part::Int, df::DataFrame, table::String)::DataFrame
    schema = stream.schemas[part]
    if schema === nothing
        stream.schemas[part] = Pair{String, Any}[name => eltype(df[!, name]) for name in names(df)]
        return df
    end
    advice = "A long-layout run that streams its results needs every column to keep its names and type from " *
        "one checkpoint to the next. Make the save field return values of one type (NaN rather than missing, " *
        "say), or run without checkpoints or in the wide layout."
    names(df) == first.(schema) || throw(ArgumentError(
        "The streamed $(table) table's columns were $(join(first.(schema), ", ")), and are now " *
        "$(join(names(df), ", ")). $(advice)"))
    for (name, T) in schema
        S = eltype(df[!, name])
        S === T && continue
        S <: T || throw(ArgumentError(
            "The streamed $(table) table's column \"$(name)\" was written as $(T), and is now $(S). $(advice)"))
        df[!, name] = convert(Vector{T}, df[!, name])
    end
    return df
end

# Whether a run that streamed its results did not finish writing them: its
# manifest still says it is streaming, or its spacecraft Feather file is missing.
function _streamed_results_unfinished(args)::Bool
    prefix = _results_bundle_prefix(args)
    isfile(prefix * "_spacecraft.feather") || return true
    manifest_path = prefix * ".manifest.toml"
    isfile(manifest_path) || return true
    return get(TOML.parsefile(manifest_path), "streaming", false) === true
end

# At a checkpoint: append the rows saved since the last one, and return the
# checkpoint's entry. A run that records no rows removes a stream an earlier run
# left, which is not its own.
function _append_results_rows!(log::_ResultsRowsLog, times_acc::Vector{Float64}, data_acc::Vector)
    _preserve_earlier_results!(log)
    if log.recording
        new = (log.logged + 1):length(times_acc)
        log.bytes = SimulationModel.IOSerialization._append_checkpoint_rows!(log.path, view(times_acc, new),
            view(data_acc, new); drop=log.drop, header=log.fresh ? merge(log.header, (token=log.token,)) : nothing)
        log.logged = length(times_acc)
        log.drop = 0
    elseif log.fresh
        rm(log.path; force=true)
    end
    log.fresh = false
    entry = (version=1, recording=log.recording, bytes=log.bytes, rows=log.logged, start_s=log.start_s,
        token=log.token)
    stream = log.stream
    stream === nothing && return entry
    return merge(entry, (version=2, stream=(bytes=Tuple(stream.bytes), released=stream.released,
        schemas=copy(stream.schemas),
        written=stream.written, first_s=stream.first_s)))
end

# The interlink history kept beside the checkpoint, so that a resumed run's link
# schedule, which is built from the history, holds the whole run as its results
# do. It follows the results rows' stream (`_append_checkpoint_rows!`): a header,
# then records appended at each checkpoint, each the samples since the last one;
# a record drops the last sample before it when the scheduler replaced it (a
# segment, or a resumed run, starts with a sample at the time the last one ended).
# The checkpoint's `interlink_history` entry records how much of the stream is its
# own (`bytes`, `samples`), its SHA-256 digest, and the stream's `token`. The live
# digest context hashes only new bytes and is never serialized. A run without interlinks
# keeps none and writes no entry.
mutable struct _InterlinkHistoryLog
    path::String
    token::String
    recording::Bool
    logged::Int            # the history's samples already in the stream
    last::Any              # the last of them, to see whether the scheduler replaced it
    bytes::Int
    fresh::Bool            # the next append starts a new stream
    start_s::Float64
    keep_earlier::Bool
    hash_context::SHA256_CTX
end

_InterlinkHistoryLog(path, token, recording, logged, last, bytes, fresh, start_s, keep_earlier) =
    _InterlinkHistoryLog(path, token, recording, logged, last, bytes, fresh, start_s, keep_earlier, SHA256_CTX())

const _INTERLINK_HISTORY_HEADER = (kind=:spaceagora_interlink_history, version=2)

_interlink_recording(model) = (record_available=SimulationModel.InterLinkModels.records_available(model),
    record_terms=model.record_terms)

function _interlink_history_checkpoint_state(checkpoint)
    runtime = get(checkpoint, :runtime_state, nothing)
    return runtime !== nothing && hasproperty(runtime, :interlink_history) ? runtime.interlink_history : nothing
end

# The history up to the checkpoint at `t`, or `nothing` and the reason it cannot be used.
function _checkpoint_interlink_history(entry, path::String, t::Float64, model; hash_context=SHA256_CTX())
    entry === nothing &&
        return nothing, "the checkpoint holds no record of it (it was written before SpaceAGORA kept it)"
    hasproperty(entry, :version) && entry.version === 2 ||
        return nothing, "the checkpoint's record of it has an unknown version"
    all(key -> hasproperty(entry, key), (:bytes, :token, :samples, :recording)) ||
        return nothing, "the checkpoint's history record is incomplete"
    entry.recording == _interlink_recording(model) || throw(ArgumentError(
        "Cannot resume interlink history with different record_available or record_terms settings. " *
        "Restore the checkpoint's recording settings or start a new run."))
    entry.bytes isa Integer && entry.bytes >= 0 ||
        return nothing, "the checkpoint's history byte length is invalid"
    if hasproperty(entry, :sha256)
        entry.sha256 isa String && occursin(r"^[0-9a-f]{64}$", entry.sha256) ||
            return nothing, "the checkpoint's history checksum is malformed"
    end
    try
        SimulationModel.IOSerialization._hash_checkpoint_rows!(hash_context, path, 0, entry.bytes)
    catch err
        err isa InterruptException && rethrow()
        return nothing, "$(path) is missing, shorter than the checkpoint records, or unreadable"
    end
    if hasproperty(entry, :sha256)
        entry.sha256 == bytes2hex(SHA.digest!(deepcopy(hash_context))) ||
            return nothing, "$(path) does not match the checkpoint's history checksum"
    end
    stream = SimulationModel.IOSerialization._read_checkpoint_rows(path, entry.bytes)
    stream === nothing && return nothing, "$(path) is missing, shorter than the checkpoint records, or unreadable"
    hasproperty(stream.header, :token) && stream.header.token == entry.token ||
        return nothing, "$(path) was written by another run than the checkpoint"
    Base.structdiff(stream.header, NamedTuple{(:token,)}) ==
        merge(_INTERLINK_HISTORY_HEADER, (recording=entry.recording,)) ||
        return nothing, "it was kept in another format"
    length(stream.data) == entry.samples && !isempty(stream.times) && stream.times[end] == t ||
        return nothing, "it does not reach the checkpoint's time"
    return stream, nothing
end

# Start the run's interlink history: a resumed run takes back the history up to
# its checkpoint and continues its stream after cutting away anything written
# after the checkpoint. When it cannot, it warns, and its history, and so its
# link schedule, starts at the resume (`history_start_s` says so).
function _start_interlink_history!(args, resumed::Bool, t_start::Float64, entry; solves::Bool=true)::_InterlinkHistoryLog
    path = _checkpoint_paths(args).interlinks
    model = args.interlink_model
    log = _InterlinkHistoryLog(path, _results_rows_token(path), model !== nothing, 0, nothing, 0, true, t_start, resumed && solves)
    (resumed && model !== nothing) || return log
    stream, reason = _checkpoint_interlink_history(entry, path, t_start, model; hash_context=log.hash_context)
    if stream === nothing
        consequence = solves ? "It covers the resumed interval only, from history_start_s = $(t_start) s." :
            "This resume solves nothing, so its returned schedule is empty and earlier schedule files are preserved."
        @warn "The resumed run's link schedule cannot include the windows before its checkpoint at t = $(t_start) s: " *
              "$(reason). $(consequence)"
        return log
    end
    if !hasproperty(entry, :sha256)
        @warn "This checkpoint's interlink history has no checksum. Its readable history is retained, " *
              "but earlier byte corruption cannot be detected. New checkpoints will record a checksum."
    end
    append!(model.history, stream.data)
    solves && SimulationModel.IOSerialization._truncate_checkpoint_rows!(path, entry.bytes)
    log.keep_earlier = false
    log.token = stream.header.token
    log.logged = length(model.history)
    log.last = isempty(model.history) ? nothing : model.history[end]
    log.bytes = entry.bytes
    log.fresh = false
    return log
end

# At a checkpoint: append the history's samples since the last record, and return
# the checkpoint's entry; `nothing` for a run without interlinks, which removes a
# stream an earlier run left.
function _append_interlink_history!(log::_InterlinkHistoryLog, model)
    if log.keep_earlier
        _keep_earlier_rows!(log.path, log.start_s)
        log.keep_earlier = false
    end
    if !log.recording
        log.fresh && rm(log.path; force=true)
        log.fresh = false
        return nothing
    end
    history = model.history
    # A sample the scheduler replaced since it was logged is logged again.
    drop = log.logged > 0 && log.logged <= length(history) && history[log.logged] !== log.last ? 1 : 0
    new = (log.logged - drop + 1):length(history)
    previous_bytes = log.fresh ? 0 : log.bytes
    next_context = log.fresh ? SHA256_CTX() : deepcopy(log.hash_context)
    next_bytes = SimulationModel.IOSerialization._append_checkpoint_rows!(log.path,
        Float64[sample.time for sample in view(history, new)], collect(view(history, new)); drop,
        header=log.fresh ? merge(_INTERLINK_HISTORY_HEADER, (token=log.token, recording=_interlink_recording(model))) : nothing)
    SimulationModel.IOSerialization._hash_checkpoint_rows!(next_context, log.path, previous_bytes, next_bytes)
    log.hash_context = next_context
    log.bytes = next_bytes
    log.logged = length(history)
    log.last = isempty(history) ? nothing : history[end]
    log.fresh = false
    return (version=2, samples=log.logged, bytes=log.bytes, token=log.token, recording=_interlink_recording(model),
        sha256=bytes2hex(SHA.digest!(deepcopy(log.hash_context))))
end

# Where the run writes its results: its results directory when its rows start at
# the start of the run, otherwise (a resumed run that could not include the rows
# before its resume) `resumed_at_<start>s` inside it, so the earlier results are
# not replaced. Every other setting is carried by field name.
function _results_destination(args, log)
    (log === nothing || log.start_s == 0.0) && return args
    settings = args.simulation_settings
    names = fieldnames(typeof(settings))
    fields = NamedTuple{names}(map(name -> getfield(settings, name), names))
    directory = joinpath(settings.results_directory, "resumed_at_$(log.start_s)s")
    return SimulationModel.SimConfig._with_configuration(args;
        simulation_settings=SimulationModel.SimulationSettings(; merge(fields, (results_directory=directory,))...))
end

# A checkpointed run's segment from `t` ends one interval later, or at the run's
# end when that comes first or is the same time to the periodic tick clock
# (`_same_run_time`): accumulated decimal intervals stop short of the end by an
# ulp (ten 0.1 s intervals reach 0.9999999999999999), and a last segment that
# short is a step the solver can reject and then abort as unstable.
function _checkpoint_segment_end(t::Float64, interval::Float64, mission_end::Float64)
    t_next = t + interval
    t_next >= mission_end && return mission_end
    return SimulationModel.SimulationCallbacks._same_run_time(t_next, mission_end) ? mission_end : t_next
end

@inline function _write_checkpoint!(args, t::Float64, u_state, solver_mode::String; is_active=nothing, callback_state=nothing,
                                     model_state=nothing, periodic_ticks=nothing, gnc_scopes=nothing,
                                     results_rows=nothing, interlink_history=nothing, maneuver_buffers=nothing, orbits=nothing, entries=nothing)
    runtime = _interlink_checkpoint_runtime_state(args, is_active)
    if callback_state !== nothing
        runtime = runtime === nothing ? (callbacks=callback_state,) : merge(runtime, (callbacks=callback_state,))
    end
    if orbits !== nothing
        runtime = runtime === nothing ? (orbits=orbits,) : merge(runtime, (orbits=orbits,))
    end
    if entries !== nothing
        runtime = runtime === nothing ? (entries=entries,) : merge(runtime, (entries=entries,))
    end
    if maneuver_buffers !== nothing
        runtime = runtime === nothing ? (maneuver_buffers=maneuver_buffers,) :
            merge(runtime, (maneuver_buffers=maneuver_buffers,))
    end
    if is_active !== nothing
        activity = _activity_checkpoint_state(is_active)
        runtime = runtime === nothing ? (activity=activity,) : merge(runtime, (activity=activity,))
    end
    if model_state !== nothing
        runtime = runtime === nothing ? (models=model_state,) : merge(runtime, (models=model_state,))
    end
    if periodic_ticks !== nothing
        runtime = runtime === nothing ? (periodic_ticks=periodic_ticks,) : merge(runtime, (periodic_ticks=periodic_ticks,))
    end
    if gnc_scopes !== nothing
        runtime = runtime === nothing ? (gnc_scopes=gnc_scopes,) : merge(runtime, (gnc_scopes=gnc_scopes,))
    end
    effector_scopes = _effector_scope_checkpoint_state(args)
    if effector_scopes !== nothing
        runtime = runtime === nothing ? (effector_scopes=effector_scopes,) : merge(runtime, (effector_scopes=effector_scopes,))
    end
    if results_rows !== nothing
        runtime = runtime === nothing ? (results_rows=results_rows,) : merge(runtime, (results_rows=results_rows,))
    end
    if interlink_history !== nothing
        runtime = runtime === nothing ? (interlink_history=interlink_history,) :
            merge(runtime, (interlink_history=interlink_history,))
    end
    return SimulationModel.IOSerialization._write_checkpoint!(
        args,
        t,
        _checkpoint_state(args, u_state),
        CHECKPOINT_SCHEMA_VERSION;
        solver_mode=solver_mode,
        runtime_state=runtime
    )
end

@inline _load_checkpoint(args) = SimulationModel.IOSerialization._load_checkpoint(args)
@inline _clear_checkpoint!(args) = SimulationModel.IOSerialization._clear_checkpoint!(args)
