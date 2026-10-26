# Dynamic effector scopes (`scoped`, `SpacecraftSelection`), resolved once when a
# run starts: after the per-spacecraft GNC is flattened and before any check that
# reads the effectors by type. The run's configuration then holds the bare
# effectors, so routing, thread whitelists, workspaces and type lookups see the
# effectors themselves, and an `EffectorScopes` table of the run positions each
# one acts on. A run with no scoped effector keeps its configuration and has no
# table, so it runs the code it ran before scopes existed.

import ..GNCRuntime

const _EffectorScopes = SimulationModel.SpacecraftModels.EffectorScopes

# Whether effector `eff_idx` of the run acts on spacecraft `sat_idx` (always, without a table).
@inline _effector_applies(scopes, eff_idx::Int, sat_idx::Int)::Bool =
    SimulationModel.SpacecraftModels._effector_applies(scopes, eff_idx, sat_idx)

# The run's scope table, or nothing (no scoped effector, or a parameter set without a run).
@inline _effector_scopes(::Nothing) = nothing
@inline _effector_scopes(p) = hasproperty(p, :args) ? _effector_scopes_of(p.args) : nothing
@inline _effector_scopes_of(args) = hasproperty(args, :dynamics_model) ? _scopes_of_model(args.dynamics_model) : nothing
@inline _scopes_of_model(model::SimulationModel.DynamicsModel) = model.effector_scopes
@inline _scopes_of_model(::Any) = nothing

# Whether some atmosphere-consuming effector acts on spacecraft `sat_idx`.
@inline _atmosphere_wanted(::Nothing, sat_idx::Int)::Bool = true
@inline function _atmosphere_wanted(scopes::_EffectorScopes, sat_idx::Int)::Bool
    mask = scopes.atmosphere
    return mask === nothing || @inbounds mask[sat_idx]
end

# The active flags of a flat batch kernel for effector `eff_idx`: the run's flags
# themselves when the effector acts on every spacecraft, otherwise the flags of
# the active spacecraft in its scope.
struct _ScopedActiveFlags{A}
    active::A
    mask::BitVector
end
@inline Base.getindex(f::_ScopedActiveFlags, i::Int) = @inbounds(f.active[i]) && @inbounds(f.mask[i])

@inline _scoped_active_flags(active, ::Nothing, eff_idx::Int) = active
@inline function _scoped_active_flags(active, scopes::_EffectorScopes, eff_idx::Int)
    @inbounds mask = scopes.masks[eff_idx]
    return mask === nothing ? active : _ScopedActiveFlags(active, mask)
end

@inline _bare_effector(e::SimulationModel.ScopedEffector) = e.effector
@inline _bare_effector(e) = e
@inline _declared_effector_scope(e::SimulationModel.ScopedEffector) = e.scope
@inline _declared_effector_scope(e) = GNCRuntime.spacecraft_scope(e)

_effector_owner(j::Int, e) = "Dynamic effector $j ($(nameof(typeof(e))))"

# Built-ins that already pick their spacecraft by an index field or keep state per
# spacecraft; scoping them is refused until they declare their own scope.
_index_bound_effector(e) =
    e isa SimulationModel.DynamicEffectors.RobotArmReactionEffector ||
    e isa SimulationModel.DynamicEffectors.ReactionWheelMomentumModel ||
    e isa SimulationModel.DynamicEffectors.PlumeSurfaceInteractionModel

# The gravity models EDG's prediction reads from the run's whole effector tuple.
_edg_predicted_gravity(e) =
    e isa SimulationModel.InverseSquaredJ2GravityModel || e isa SimulationModel.NBodyGravityModel

# Guidance whose prediction applies the run's gravity models to every spacecraft
# (EDG's). The guidance types say so where they are defined, through the private
# GNCRuntime predicate, so this file imports no module that loads after the engine.
_edg_guidance(g)::Bool = GNCRuntime._prediction_uses_run_gravity(g) === true

"""
    _resolve_effector_scopes(args) -> SimulationConfiguration

Resolve the scopes of `args`' dynamic effectors (a `scoped` entry's, or an
effector type's own `spacecraft_scope`) against its spacecraft, refusing scopes
the run cannot honour, and return the configuration the run uses: `args` itself
when no effector is scoped, otherwise one whose `DynamicsModel` holds the bare
effectors and an `EffectorScopes` table. A scope that names every spacecraft is
the same as none. Idempotent: a configuration that already holds a table is
checked and returned.
"""
function _resolve_effector_scopes(args)
    model = args.dynamics_model
    entries = model.dynamic_effectors
    spacecraft = model.spacecraft
    n = length(spacecraft)
    if model.effector_scopes !== nothing
        table = model.effector_scopes
        any(e -> e isa SimulationModel.ScopedEffector, entries) && throw(ArgumentError(
            "The dynamics model holds both scoped entries and resolved effector scopes; build it from the " *
            "scoped entries alone."))
        (table.n_spacecraft == n && length(table.masks) == length(entries)) || throw(ArgumentError(
            "The dynamics model's effector scopes were resolved for $(table.n_spacecraft) spacecraft and " *
            "$(length(table.masks)) effectors, but it has $n spacecraft and $(length(entries)) effectors. " *
            "Build it again from the scoped entries."))
        _check_scoped_run(args, entries, table)
        return args
    end
    declared = Any[_declared_effector_scope(e) for e in entries]
    if !any(e -> e isa SimulationModel.ScopedEffector, entries) && all(s -> s isa GNCRuntime.AllSpacecraft, declared)
        return args
    end
    bare = map(_bare_effector, entries)
    positions = Vector{Union{Nothing, Vector{Int}}}(undef, length(entries))
    masks = Vector{Union{Nothing, BitVector}}(undef, length(entries))
    for (j, entry) in enumerate(entries)
        effector = bare[j]
        owner = _effector_owner(j, effector)
        scope = declared[j]
        if entry isa SimulationModel.ScopedEffector
            effector isa SimulationModel.ScopedEffector && throw(ArgumentError(
                "$owner is scoped twice; give scoped the effector itself."))
            GNCRuntime.spacecraft_scope(effector) isa GNCRuntime.AllSpacecraft || throw(ArgumentError(
                "$owner already declares its scope through spacecraft_scope; the run applies that scope, so give " *
                "the effector without scoped."))
        end
        scope isa GNCRuntime.SpacecraftScope || throw(ArgumentError(
            "$owner has the scope $(repr(scope)); a scope is AllSpacecraft(), BoundSpacecraft(positions) or " *
            "SpacecraftSelection(...)."))
        scope isa GNCRuntime.CoupledSpacecraft && throw(ArgumentError(
            "$owner has a CoupledSpacecraft scope: a force model acts on one spacecraft at a time, and forces " *
            "between spacecraft are interlinks (link_force). Give the spacecraft it acts on as BoundSpacecraft " *
            "or SpacecraftSelection."))
        if scope isa GNCRuntime.AllSpacecraft
            positions[j] = nothing
            masks[j] = nothing
            continue
        end
        selected = SimulationModel.resolve_spacecraft_scope(scope, spacecraft; owner=owner)
        isempty(selected) && throw(ArgumentError(
            "$owner selects no spacecraft of the run, so it would never act. Remove it, or widen its scope."))
        if length(selected) == n
            positions[j] = nothing
            masks[j] = nothing
        else
            mask = falses(n)
            mask[selected] .= true
            positions[j] = selected
            masks[j] = mask
        end
    end
    all(isnothing, masks) && return SimulationModel.SimConfig._with_configuration(args;
        dynamics_model=SimulationModel.DynamicsModel(spacecraft, bare))
    table = _EffectorScopes(positions, masks, declared, n, _atmosphere_mask(bare, masks, n))
    _check_scoped_run(args, bare, table)
    return SimulationModel.SimConfig._with_configuration(args;
        dynamics_model=SimulationModel.DynamicsModel(spacecraft, bare, table))
end

# The spacecraft some atmosphere-consuming effector acts on, or nothing when one acts on all.
function _atmosphere_mask(effectors::Tuple, masks, n::Int)
    wanted = falses(n)
    for (j, effector) in enumerate(effectors)
        SimulationModel.environment_requirements(effector).atmosphere || continue
        masks[j] === nothing && return nothing
        wanted .|= masks[j]
    end
    return wanted
end

# Refusals that depend on the whole run, for a run with a scope table.
function _check_scoped_run(args, effectors::Tuple, table)
    spacecraft = args.dynamics_model.spacecraft
    n = length(spacecraft)
    scoped_positions = [j for j in eachindex(table.masks) if table.masks[j] !== nothing]
    for j in scoped_positions
        effector = effectors[j]
        owner = _effector_owner(j, effector)
        _index_bound_effector(effector) && throw(ArgumentError(
            "$owner picks its spacecraft itself (by its spacecraft index or its per-spacecraft state), so it " *
            "cannot be scoped. Give it without scoped."))
    end
    gravity = [any(j -> SimulationModel.is_gravitational(effectors[j]) && _effector_applies(table, j, i),
                   eachindex(table.masks)) for i in 1:n]
    if any(gravity) && !all(gravity)
        without = findall(!, gravity)
        throw(ArgumentError(
            "The run's scoped gravity leaves the spacecraft at run positions $(join(without, ", ")) without any " *
            "gravitational effector while the others have one. Scope another gravity model to them, or widen " *
            "the scope."))
    end
    owned = _external_owned_indices(args)
    for j in scoped_positions, i in table.positions[j]
        i in owned && throw(ArgumentError(
            "$(_effector_owner(j, effectors[j])) names the spacecraft at run position $i in its scope, but that " *
            "spacecraft is propagated externally and dynamic effectors never act on it."))
    end
    for (i, sc) in enumerate(spacecraft)
        articulated = SimulationModel.articulated_has_moving_joints(sc)
        (articulated || !isempty(sc.attachments)) && throw(ArgumentError(
            "Scoped dynamic effectors are not supported in a run with " *
            "$(articulated ? "an articulated spacecraft (moving joints)" : "compliant attachments") " *
            "(spacecraft at run position $i): that path sums the effectors per body. Give the effectors " *
            "without scoped."))
    end
    _any_robot_arm_effector(args) && throw(ArgumentError(
        "Scoped dynamic effectors are not supported in a run with a robot arm, whose link gravity sums every " *
        "effector. Give the effectors without scoped."))
    if any(_edg_guidance, args.guidance_model.guidance_effectors)
        for j in scoped_positions
            _edg_predicted_gravity(effectors[j]) && throw(ArgumentError(
                "$(_effector_owner(j, effectors[j])) is scoped in a run with EDG guidance, whose prediction " *
                "applies the run's gravity models to every spacecraft. Give it without scoped."))
        end
    end
    return nothing
end

# The checkpoint entry of a run with scoped effectors, `effector_scopes`: each
# effector's kind (as the interlink entry names kinds) and the run positions it
# acts on (`nothing` for every spacecraft). A run without scoped effectors writes
# no entry, so its checkpoint files keep the runtime keys they had.
function _effector_scope_checkpoint_state(args)
    scopes = _effector_scopes_of(args)
    scopes === nothing && return nothing
    return (version=1, effectors=[(kind=_interlink_kind_identity(typeof(effector)),
        positions=scopes.positions[j] === nothing ? nothing : copy(scopes.positions[j]))
        for (j, effector) in enumerate(args.dynamics_model.dynamic_effectors)])
end

# Before a resume restores anything: the checkpoint and the run must agree on the
# effectors' scopes, both having none, or both the same kinds acting on the same
# positions.
function _check_effector_scope_checkpoint_state(args, checkpoint)
    runtime = get(checkpoint, :runtime_state, nothing)
    entry = runtime !== nothing && hasproperty(runtime, :effector_scopes) ? runtime.effector_scopes : nothing
    current = _effector_scope_checkpoint_state(args)
    entry === nothing && current === nothing && return nothing
    entry === nothing && throw(ArgumentError(
        "This run has scoped dynamic effectors, but the checkpoint records no effector scopes: it was written by a " *
        "run without them, or by a SpaceAGORA that did not record them. Start the run again rather than resume it."))
    current === nothing && throw(ArgumentError(
        "The checkpoint was written by a run with scoped dynamic effectors, but this run has none; resume it with " *
        "the same scoped effectors."))
    hasproperty(entry, :version) && entry.version === 1 && hasproperty(entry, :effectors) || throw(ArgumentError(
        "Effector scope checkpoint state version 1 is required to resume this run."))
    if !isequal(entry.effectors, current.effectors)
        length(entry.effectors) == length(current.effectors) || throw(ArgumentError(
            "The dynamic effectors' spacecraft scopes differ from the checkpoint's: the checkpoint has " *
            "$(length(entry.effectors)) dynamic effectors and this run $(length(current.effectors))."))
        j = findfirst(k -> !isequal(entry.effectors[k], current.effectors[k]), eachindex(current.effectors))
        describe(positions) = positions === nothing ? "every spacecraft" : "run positions $(join(positions, ", "))"
        throw(ArgumentError(
            "The dynamic effectors' spacecraft scopes differ from the checkpoint's: effector $j " *
            (isequal(entry.effectors[j].kind, current.effectors[j].kind) ?
                "acts on $(describe(current.effectors[j].positions)) here and on " *
                "$(describe(entry.effectors[j].positions)) in the checkpoint." :
                "is a $(current.effectors[j].kind.name) here and a $(entry.effectors[j].kind.name) in the checkpoint.")))
    end
    return nothing
end
