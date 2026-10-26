# The run's GNC scopes (`GNCRunScopes`): the spacecraft each guidance, navigation
# and control model acts on, taken once when the run starts (`_start_run_models!`),
# after the run-start initializers (the RPO pair sets its positions there, from the
# spacecraft ids) and, on a resumed run, as its checkpoint's `gnc_scopes` entry is
# compared with them (`_check_gnc_scope_baseline`). From then on:
# - every scope read compares the scope with the one taken here and stops the run
#   on a change (SimulationCallbacks._check_run_scope): the periodic ticks at each
#   tick, and the event-driven control models once per accepted step;
# - a SpacecraftSelection acts on the positions resolved here, its `where`
#   called once per spacecraft;
# - the right-hand side and the per-spacecraft readers call a control or guidance
#   model's hooks only for the spacecraft its scope names (the member lists), when
#   some model of the group acts on fewer than every spacecraft.
# A checkpoint records the same scopes as its `gnc_scopes` v1 entry
# (`_gnc_scope_record`).

# Test switch: `false` builds no member lists, so every model's hooks are called
# for every spacecraft, as before SpaceAGORA 0.7.0. The scopes are still taken and
# checked.
const _GNC_MEMBER_LISTS = Ref(true)

function _gnc_run_scopes(args)
    spacecraft = args.dynamics_model.spacecraft
    n = length(spacecraft)
    groups = _gnc_groups(args)
    slots(group) = SimulationModel.GNCScopeSlot[
        _gnc_scope_slot(GNCRuntime.spacecraft_scope(m), spacecraft, _gnc_scope_owner(group, j, m))
        for (j, m) in enumerate(getproperty(groups, group))]
    navigation, guidance, control = slots(:navigation), slots(:guidance), slots(:control)
    return SimulationModel.GNCRunScopes(navigation, guidance, control, _gnc_members(guidance, n),
        _gnc_members(control, n))
end

# One model's slot: its scope as declared now, and the positions it acts on. A
# BoundSpacecraft or CoupledSpacecraft acts where its ticks do (the indices within
# the run, the primary); a SpacecraftSelection where it resolves to, refused if it
# cannot be resolved; and the run keeps its own copy of the selection's criteria.
_gnc_scope_slot(scope::GNCRuntime.AllSpacecraft, spacecraft, owner) =
    SimulationModel.GNCScopeSlot(owner, scope, collect(1:length(spacecraft)))
_gnc_scope_slot(scope::Union{GNCRuntime.BoundSpacecraft, GNCRuntime.CoupledSpacecraft}, spacecraft, owner) =
    SimulationModel.GNCScopeSlot(owner, scope,
        collect(Int, SimulationModel.SimulationCallbacks._tick_positions(scope, length(spacecraft))))
function _gnc_scope_slot(scope::SimulationModel.SpacecraftSelection, spacecraft, owner)
    positions = SimulationModel.resolve_spacecraft_scope(scope, spacecraft; owner=owner)
    return SimulationModel.GNCScopeSlot(owner, _owned_selection(scope), positions)
end
_gnc_scope_slot(scope, spacecraft, owner) = throw(ArgumentError(
    "$owner returns $(repr(scope)) from spacecraft_scope once its run-start initializers ran; a guidance, " *
    "navigation or control model's scope is AllSpacecraft(), BoundSpacecraft(indices), " *
    "CoupledSpacecraft(primary, partners) or a SpacecraftSelection."))

_owned_selection(s) = SimulationModel.SpacecraftSelection(_owned_copy(s.positions), _owned_copy(s.ids),
    _owned_copy(s.planes), _owned_copy(s.slots), s.where)
_owned_copy(::Nothing) = nothing
_owned_copy(values::AbstractVector) = copy(values)

# For each spacecraft, the positions in the group of the models whose scope names
# it, in configuration order; `nothing` when every model acts on every spacecraft.
function _gnc_members(slots, n::Int)
    (_GNC_MEMBER_LISTS[] && any(slot -> length(slot.positions) != n, slots)) || return nothing
    members = [Int[] for _ in 1:n]
    for (j, slot) in enumerate(slots), i in slot.positions
        push!(members[i], j)
    end
    return members
end

# The checkpoint's `gnc_scopes` v1 entry (as `_gnc_scope_baseline` writes it) from
# the run's scopes. A SpacecraftSelection is recorded at the positions the run
# resolved it to, so its `where` is not called again.
function _gnc_scope_record(args, scopes)
    spacecraft = args.dynamics_model.spacecraft
    entry(slot) = slot.declared isa SimulationModel.SpacecraftSelection ?
        (scope=:bound, positions=copy(slot.positions)) : _gnc_scope_identity(slot.declared, spacecraft, slot.owner)
    return (version=1, navigation=Any[entry(s) for s in scopes.navigation],
        guidance=Any[entry(s) for s in scopes.guidance], control=Any[entry(s) for s in scopes.control])
end
