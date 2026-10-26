# Spacecraft selections: the spacecraft of a run chosen by run position, id,
# constellation plane or slot, or a predicate, resolved to run positions once per
# run. A selection is a `SpacecraftScope` (GNCRuntime), the vocabulary GNC models
# use for the spacecraft they act on; dynamic effectors take one through `scoped`.
# Included by `SpacecraftModels` (model.jl), after the constellation builder;
# the scope types come from GNCRuntime through constellation_configuration.jl.

"""
    SpacecraftSelection(; positions=nothing, ids=nothing, planes=nothing, slots=nothing, where=nothing)

The spacecraft of a run that meet every criterion given; give at least one.

- `positions`: run positions, from 1, in the order of `DynamicsModel.spacecraft`.
- `ids`: spacecraft `id`s.
- `planes`: constellation planes, matched against `spacecraft.constellation_slot.plane`
  (`build_constellation` sets the slots).
- `slots`: constellation slots, as `ConstellationSlot`s or `(plane, slot)` pairs.
- `where`: a function `(position, spacecraft) -> Bool`, called once per spacecraft
  when the run starts.

Each criterion takes one value or a collection (`planes=1`, `ids=101:108`). A run
resolves the selection to run positions once, when it starts, and refuses a
selection that names an id no spacecraft or two spacecraft hold, a plane or slot
no spacecraft holds, or a position outside the run. Use it as the scope of
`scoped`, or return it from a GNC model's `spacecraft_scope`.

Two selections are equal when their criteria are equal and their `where` is the
same function object (`===`). A GNC model's scope must not change during the run,
so such a model keeps its selection in a field and returns that object: a selection
built anew at each call, with a new closure for `where`, is a different scope.
"""
struct SpacecraftSelection{W} <: SpacecraftScope
    positions::Union{Nothing, Vector{Int}}
    ids::Union{Nothing, Vector{Int}}
    planes::Union{Nothing, Vector{Int}}
    slots::Union{Nothing, Vector{ConstellationSlot}}
    where::W
end

function SpacecraftSelection(; positions=nothing, ids=nothing, planes=nothing, slots=nothing, where=nothing)
    all(isnothing, (positions, ids, planes, slots, where)) && throw(ArgumentError(
        "A SpacecraftSelection needs at least one of positions, ids, planes, slots or where."))
    positions = _selection_integers(positions, "positions"; minimum=1)
    ids = _selection_integers(ids, "ids")
    planes = _selection_integers(planes, "planes"; minimum=1)
    slots = slots === nothing ? nothing : ConstellationSlot[_selection_slot(s) for s in _selection_slot_values(slots)]
    where === nothing || hasmethod(where, Tuple{Int, SpacecraftModel}) ||
        throw(ArgumentError("A SpacecraftSelection's where is a function (position, spacecraft) -> Bool."))
    return SpacecraftSelection(positions, ids, planes, slots, where)
end

Base.:(==)(a::SpacecraftSelection, b::SpacecraftSelection) =
    (a.positions, a.ids, a.planes, a.slots) == (b.positions, b.ids, b.planes, b.slots) && a.where === b.where

_selection_values(value) = value isa Union{AbstractVector, Tuple, AbstractRange, AbstractSet} ? value : (value,)
# Slots: one ConstellationSlot or one (plane, slot) pair is a single slot.
_selection_slot_values(value) = _selection_values(value)
_selection_slot_values(value::ConstellationSlot) = (value,)
_selection_slot_values(value::Tuple{Integer, Integer}) = (value,)

function _selection_integers(value, name::String; minimum::Union{Nothing, Int}=nothing)
    value === nothing && return nothing
    values = Int[]
    for v in _selection_values(value)
        v isa Integer || throw(ArgumentError("A SpacecraftSelection's $name are integers; got $(repr(v))."))
        minimum === nothing || v >= minimum || throw(ArgumentError(
            "A SpacecraftSelection's $name are numbered from $minimum; got $v."))
        push!(values, Int(v))
    end
    return values
end

_selection_slot(slot::ConstellationSlot) = slot
_selection_slot(slot::Tuple{Integer, Integer}) = ConstellationSlot(slot...)
_selection_slot(slot) = throw(ArgumentError(
    "A SpacecraftSelection's slots are ConstellationSlot(plane, slot)s or (plane, slot) pairs; got $(repr(slot))."))

"""
    resolve_spacecraft_scope(scope, spacecraft; owner="the model") -> Vector{Int}

The run positions, sorted and without repeats, that `scope` names among the run's
`spacecraft` (a vector in run order): every position for `AllSpacecraft()`, the
`BoundSpacecraft` indices, or the spacecraft a `SpacecraftSelection` selects.
Throws an `ArgumentError` naming `owner` for a position outside the run, an id no
spacecraft or two spacecraft hold, a plane or slot no spacecraft holds, or a `where`
that throws or does not return a `Bool`. A `CoupledSpacecraft` is not resolved here.
"""
resolve_spacecraft_scope(::AllSpacecraft, spacecraft::AbstractVector; owner::AbstractString="the model") =
    collect(1:length(spacecraft))

function resolve_spacecraft_scope(scope::BoundSpacecraft, spacecraft::AbstractVector; owner::AbstractString="the model")
    return _positions_in_run(scope.indices, length(spacecraft), owner)
end

resolve_spacecraft_scope(scope::CoupledSpacecraft, spacecraft::AbstractVector; owner::AbstractString="the model") =
    throw(ArgumentError("$owner has a CoupledSpacecraft scope, which names a primary and the partners it reads, " *
                        "not a set of spacecraft to resolve."))

function resolve_spacecraft_scope(scope::SpacecraftSelection, spacecraft::AbstractVector; owner::AbstractString="the model")
    n = length(spacecraft)
    selected = trues(n)
    if scope.positions !== nothing
        keep = falses(n)
        keep[_positions_in_run(scope.positions, n, owner)] .= true
        selected .&= keep
    end
    if scope.ids !== nothing
        keep = falses(n)
        for id in unique(scope.ids)
            holders = [i for (i, sc) in enumerate(spacecraft) if sc.id == id]
            isempty(holders) && throw(ArgumentError(
                "$owner selects spacecraft id $id, which no spacecraft of the run has."))
            length(holders) == 1 || throw(ArgumentError(
                "$owner selects spacecraft id $id, which the spacecraft at run positions $(join(holders, ", ")) " *
                "share. Select them by position, or give each spacecraft its own id."))
            keep[only(holders)] = true
        end
        selected .&= keep
    end
    if scope.planes !== nothing || scope.slots !== nothing
        laid_out = [sc.constellation_slot for sc in spacecraft]
        all(isnothing, laid_out) && throw(ArgumentError(
            "$owner selects constellation planes or slots, but no spacecraft of the run has a constellation_slot. " *
            "Lay the constellation out with build_constellation, or select by position or id."))
        if scope.planes !== nothing
            keep = falses(n)
            for plane in unique(scope.planes)
                holders = [i for (i, slot) in enumerate(laid_out) if slot !== nothing && slot.plane == plane]
                isempty(holders) && throw(ArgumentError(
                    "$owner selects constellation plane $plane, which no spacecraft of the run is in."))
                keep[holders] .= true
            end
            selected .&= keep
        end
        if scope.slots !== nothing
            keep = falses(n)
            for wanted in unique(scope.slots)
                holders = [i for (i, slot) in enumerate(laid_out) if slot == wanted]
                isempty(holders) && throw(ArgumentError(
                    "$owner selects constellation slot (plane $(wanted.plane), slot $(wanted.slot)), which no " *
                    "spacecraft of the run is in."))
                keep[holders] .= true
            end
            selected .&= keep
        end
    end
    if scope.where !== nothing
        for (i, sc) in enumerate(spacecraft)
            answer = try
                scope.where(i, sc)
            catch err
                throw(ArgumentError("$owner selects spacecraft with a where function that threw for the " *
                                    "spacecraft at run position $i: $(sprint(showerror, err))"))
            end
            answer isa Bool || throw(ArgumentError(
                "$owner selects spacecraft with a where function that returned $(repr(answer)) for the spacecraft " *
                "at run position $i; it must return true or false."))
            answer || (selected[i] = false)
        end
    end
    return findall(selected)
end

function _positions_in_run(positions, n::Int, owner::AbstractString)::Vector{Int}
    outside = [i for i in positions if !(1 <= i <= n)]
    isempty(outside) || throw(ArgumentError(
        "$owner names spacecraft $(join(outside, ", ")) in its scope, outside the run's positions 1 to $n, so it " *
        "would never act there. Use run positions, from 1 in the order of DynamicsModel.spacecraft."))
    return sort!(unique(collect(Int, positions)))
end
