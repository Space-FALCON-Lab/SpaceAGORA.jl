# Spacecraft scope and inactive spacecraft, as the GNC tick loops apply them
# (GNCRuntime declarations). A periodic tick reads the model's scope once and
# visits only the run positions it names (`_tick_positions`), in ascending order,
# so a tick costs the size of the scope rather than the number of spacecraft. A
# control model that consumes a link has its scope read once more, first, to
# decide whether its producers prepare it (`_ticks_any`).
#
# When the run starts, once the run-start initializers ran, the engine takes each
# model's scope (the run's GNC scopes, `GNCRunScopes`). From then on the scope is
# fixed: each read compares it with the one the run took (`_check_run_scope`), and
# a SpacecraftSelection is visited at the positions the run resolved it to.
using ...GNCRuntime: SpacecraftScope, AllSpacecraft, BoundSpacecraft, CoupledSpacecraft,
    spacecraft_scope, skip_inactive, release_command!
using ..SpacecraftModels: SpacecraftSelection
using ..ConfigTypes: GNCScopeSlot, GNCRunScopes

# Whether a model with this scope ticks for run position `spacecraft_index`.
@inline _ticks_for(::AllSpacecraft, spacecraft_index::Int) = true
@inline _ticks_for(scope::BoundSpacecraft, spacecraft_index::Int) = spacecraft_index in scope.indices
@inline _ticks_for(scope::CoupledSpacecraft, spacecraft_index::Int) = spacecraft_index == scope.primary

# The run positions a scope names (none for every spacecraft).
_scope_indices(::AllSpacecraft) = ()
_scope_indices(scope::BoundSpacecraft) = scope.indices
_scope_indices(scope::CoupledSpacecraft) = (scope.primary, scope.partners...)

# The run positions, among `n` spacecraft, that a model with this scope ticks for:
# the positions `_ticks_for` accepts, in ascending order and each once. That is
# every position, the `BoundSpacecraft` indices within `1:n`, or the
# `CoupledSpacecraft` primary. Indices already strictly increasing and within the
# run are returned as they are; any others are sorted, deduplicated and filtered
# into a new vector.
@inline _tick_positions(::AllSpacecraft, n::Int) = Base.OneTo(n)
@inline function _tick_positions(scope::BoundSpacecraft, n::Int)
    _increasing_within(scope.indices, n) && return scope.indices
    return sort!(unique!(Int[i for i in scope.indices if 1 <= i <= n]))
end
@inline _tick_positions(scope::CoupledSpacecraft, n::Int) = 1 <= scope.primary <= n ? (scope.primary,) : ()

# A SpacecraftSelection names spacecraft only once the run resolved it, when it
# starts; outside a run (direct calls) it names none that a tick could visit.
const _UNRESOLVED_SELECTION = "A GNC model's SpacecraftSelection is resolved when the run starts; " *
    "only a running simulation can say which spacecraft it names."
_tick_positions(::SpacecraftSelection, n::Int) = throw(ArgumentError(_UNRESOLVED_SELECTION))
_ticks_for(::SpacecraftSelection, spacecraft_index::Int) = throw(ArgumentError(_UNRESOLVED_SELECTION))
_scope_indices(::SpacecraftSelection) = ()

@inline function _increasing_within(indices::Tuple, n::Int)
    previous = 0
    for i in indices
        previous < i <= n || return false
        previous = i
    end
    return true
end

# One spacecraft's step of a GNC tick loop, for a spacecraft in the model's
# scope: the tick, or `release_command!` in its place when the model declares
# `skip_inactive` and the spacecraft is inactive.
@inline function _gnc_tick_in_scope!(tick!::F, model, u, p, t, sat_idx::Int) where {F}
    if skip_inactive(model) && !p.is_active[sat_idx]
        release_command!(model, t, sat_idx)
    else
        tick!(model, u, p, t, sat_idx)
    end
    return nothing
end

# One spacecraft's step for any spacecraft, as the callbacks that visit one
# spacecraft at a time take it. A model whose scope leaves spacecraft `sat_idx`
# out is not called, and the step returns `false`. Otherwise it returns `true`
# after the step above.
@inline function _gnc_tick!(tick!::F, model, u, p, t, sat_idx::Int) where {F}
    _ticks_for(spacecraft_scope(model), sat_idx) || return false
    _gnc_tick_in_scope!(tick!, model, u, p, t, sat_idx)
    return true
end

# A periodic tick of `model` for every spacecraft in its scope, in ascending order.
# `slot` is the model's entry in the run's GNC scopes (`_run_scope_slot`), or
# `nothing` outside a run.
@inline function _gnc_tick_scope!(tick!::F, model, u, p, t, num_sats::Int, slot=nothing) where {F}
    scope = spacecraft_scope(model)
    _check_run_scope(slot, scope, num_sats)
    @inbounds for sat_idx in _run_tick_positions(slot, scope, num_sats)
        _gnc_tick_in_scope!(tick!, model, u, p, t, sat_idx)
    end
    return nothing
end

# The run's GNC scopes, set when the run starts: `nothing` before, and for
# parameters that are not a run's.
@inline _run_scopes(p) = hasproperty(p, :shared_buffers) && hasproperty(p.shared_buffers, :gnc_scopes) ?
    p.shared_buffers.gnc_scopes[] : nothing

# The run's entry for model `j` of `group` (`:navigation`, `:guidance` or
# `:control`, in configuration order), or `nothing` before the run starts.
@inline function _run_scope_slot(p, group::Symbol, j::Int)
    scopes = _run_scopes(p)
    scopes === nothing && return nothing
    return getfield(scopes, group)[j]
end

# The models of `group` (`:guidance` or `:control`, whose tuple is `models`) that
# act on spacecraft `sat_idx` in this run, as their positions in the tuple, in
# configuration order: from the run's member lists, or every position when the
# run keeps none (every model acts on every spacecraft, or the run has not
# started). The per-spacecraft readers and passes call a model's hooks for those.
@inline function _scope_members(p, group::Symbol, models, sat_idx::Int)
    scopes = _run_scopes(p)
    members = scopes === nothing ? nothing :
        group === :control ? scopes.control_members : scopes.guidance_members
    return members === nothing ? eachindex(models) : members[sat_idx]
end

# Whether model `j` of `group` acts on spacecraft `sat_idx`: as the run took its
# scope when it started, or as it reads now outside a run.
@inline function _acts_on(p, group::Symbol, j::Int, model, sat_idx::Int)
    slot = _run_scope_slot(p, group, j)
    slot === nothing && return _ticks_for(spacecraft_scope(model), sat_idx)
    return insorted(sat_idx, slot.positions)
end

# The positions a tick visits: the scope's, or a SpacecraftSelection's as the
# run resolved it when it started.
@inline _run_tick_positions(slot, scope, n::Int) = _tick_positions(scope, n)
@inline _run_tick_positions(slot::GNCScopeSlot, ::SpacecraftSelection, n::Int) = slot.positions

# A model's scope is fixed once the run-start initializers ran: a scope read later
# that differs from the one the run took then stops the run. The same kind of
# scope naming the same spacecraft is the same scope, so a BoundSpacecraft whose
# indices come in another order passes. A SpacecraftSelection is compared with
# the run's own copy of it, so a selection changed in place does not pass.
@inline _check_run_scope(::Nothing, scope, n::Int) = nothing
@inline function _check_run_scope(slot::GNCScopeSlot, scope, n::Int)
    _same_run_scope(slot, slot.declared, scope, n) && return nothing
    throw(ArgumentError(_scope_change_message(slot, scope)))
end
@inline _same_run_scope(slot, ::AllSpacecraft, ::AllSpacecraft, n::Int) = true
@inline _same_run_scope(slot, a::BoundSpacecraft, b::BoundSpacecraft, n::Int) =
    a.indices == b.indices || _same_positions(slot.positions, _tick_positions(b, n))
@inline _same_run_scope(slot, a::CoupledSpacecraft, b::CoupledSpacecraft, n::Int) = a == b
@inline _same_run_scope(slot, a::SpacecraftSelection, b::SpacecraftSelection, n::Int) = a == b
@inline _same_run_scope(slot, a, b, n::Int) = false

function _same_positions(positions::Vector{Int}, other)
    length(positions) == length(other) || return false
    for (a, b) in zip(positions, other)
        a == b || return false
    end
    return true
end

# A scope as a user writes it, without its type parameters (`BoundSpacecraft((1,))`,
# not `BoundSpacecraft{1}((1,))`), for messages.
_scope_text(::AllSpacecraft) = "AllSpacecraft()"
_scope_text(scope::BoundSpacecraft) = "BoundSpacecraft($(scope.indices))"
_scope_text(scope::CoupledSpacecraft) = "CoupledSpacecraft($(scope.primary), $(scope.partners))"
function _scope_text(scope::SpacecraftSelection)
    criteria = ["$(name)=$(repr(getfield(scope, name)))" for name in (:positions, :ids, :planes, :slots)
        if getfield(scope, name) !== nothing]
    if scope.where !== nothing
        name = string(scope.where isa Function ? nameof(scope.where) : nameof(typeof(scope.where)))
        push!(criteria, startswith(name, "#") ? "where=<function>" : "where=$name")
    end
    return "SpacecraftSelection($(join(criteria, ", ")))"
end
_scope_text(scope) = repr(scope)

_scope_change_message(slot, scope) = _scope_change_message(slot, slot.declared, scope)
_scope_change_message(slot, declared, scope) =
    "$(slot.owner)'s spacecraft_scope changed during the run: it was $(_scope_text(declared)) when the run " *
    "started and is $(_scope_text(scope)) now. A guidance, navigation or control model's scope is fixed once its " *
    "run-start initializers ran; set it before then (in the model, in bind_spacecraft or in an initializer)."
# A selection with the same criteria and another `where` object, as a scope method
# that builds its selection, and a new closure, at each call returns: the advice is
# to return one object.
function _scope_change_message(slot, declared::SpacecraftSelection, scope::SpacecraftSelection)
    same_criteria = (declared.positions, declared.ids, declared.planes, declared.slots) ==
        (scope.positions, scope.ids, scope.planes, scope.slots)
    (same_criteria && declared.where !== scope.where) ||
        return invoke(_scope_change_message, Tuple{Any, Any, Any}, slot, declared, scope)
    return "$(slot.owner)'s spacecraft_scope changed during the run: it returns a SpacecraftSelection with " *
        "the criteria it had when the run started but another where function. A selection's where is " *
        "compared by identity (===), so a closure built anew at each call is a different scope. Keep the " *
        "selection in the model and return that object at every call, or give where a function that is " *
        "the same object at every call."
end

# Whether `model`'s command for spacecraft `sat_idx` is released at this time, as
# `_gnc_tick!` releases it: the spacecraft is in the model's scope and inactive, and
# the model declares `skip_inactive`. A path that writes the command outside the
# tick (the burn scheduler's scheduling and burn edges) releases it again after its
# bookkeeping, so nothing undoes a release.
@inline _released(model, p, sat_idx::Int) =
    _ticks_for(spacecraft_scope(model), sat_idx) && skip_inactive(model) && !p.is_active[sat_idx]

# Whether a tick of `model` at this time runs its tick for some spacecraft, rather
# than skipping or releasing every one: the condition for preparing it first
# (`before_reference_control!` for its linked producers). It reads the scope
# itself, before the producers prepare; the tick reads it again afterwards (see
# `_periodic_control_tick!`). Order and repeats do not matter to the answer, so a
# `BoundSpacecraft`'s indices are checked as given, without the sorted copy
# `_tick_positions` builds for unsorted ones.
function _ticks_any(model, p, num_sats::Int, slot=nothing)
    scope, skip = spacecraft_scope(model), skip_inactive(model)
    _check_run_scope(slot, scope, num_sats)
    for sat_idx in _any_order_positions(slot, scope, num_sats)
        1 <= sat_idx <= num_sats || continue
        !(skip && !p.is_active[sat_idx]) && return true
    end
    return false
end
@inline _any_order_positions(slot, scope::BoundSpacecraft, n::Int) = scope.indices
@inline _any_order_positions(slot, scope, n::Int) = _run_tick_positions(slot, scope, n)

"""
    check_spacecraft_scopes(args) -> nothing

Before a run: refuse it with an `ArgumentError` when a guidance, navigation or
control model's `spacecraft_scope` is not a `SpacecraftScope`, its
`skip_inactive` is not a `Bool`, or its scope names a run position outside `1:n`
for the run's `n` spacecraft, where the model would never act. (SpaceAGORA 0.3
warned about such a position instead.) A `SpacecraftSelection` is checked when
the run resolves it, as the run starts.
"""
function check_spacecraft_scopes(args)
    n = length(args.dynamics_model.spacecraft)
    models = (args.guidance_model.guidance_effectors..., args.navigation_model.navigation_effectors...,
              args.control_model.control_effectors...)
    for model in models
        T = typeof(model)
        scope = spacecraft_scope(model)
        scope isa Union{AllSpacecraft, BoundSpacecraft, CoupledSpacecraft, SpacecraftSelection} || throw(ArgumentError(
            "spacecraft_scope($(nameof(T))) must return AllSpacecraft(), BoundSpacecraft(indices), " *
            "CoupledSpacecraft(primary, partners) or a SpacecraftSelection, not $(repr(scope))."))
        inactive = skip_inactive(model)
        inactive isa Bool || throw(ArgumentError(
            "skip_inactive($(nameof(T))) must return true or false, not $(repr(inactive))."))
        outside = [i for i in _scope_indices(scope) if !(1 <= i <= n)]
        isempty(outside) || throw(ArgumentError(
            "$(nameof(T)) names spacecraft $(join(outside, ", ")) in its spacecraft_scope, outside the run's " *
            "positions 1 to $(n), so the model would never act there. Set its spacecraft indices to run positions."))
    end
    return nothing
end
