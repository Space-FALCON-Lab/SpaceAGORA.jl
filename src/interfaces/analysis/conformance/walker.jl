# Walkers: the concrete types a package defines, and the state reachable from an
# object.

"""
    package_types(T::Type, root::Module=SpaceAGORA) -> Vector{Any}

Every concrete type that subtypes `T` and is defined in `root` or one of its
submodules (types are found by name, so each appears once). Core's conformance
tests use it to run a check on every built-in model; a package can use it the same
way on its own module.
"""
function package_types(T::Type, root::Module=parentmodule(@__MODULE__))
    found = Any[]
    seen = Set{Module}()
    function walk(m::Module)
        m in seen && return
        push!(seen, m)
        for n in names(m; all=true)
            isdefined(m, n) || continue
            v = getfield(m, n)
            if v isa Module
                parentmodule(v) === m && v !== m && walk(v)
            elseif v isa Union{DataType, UnionAll} && parentmodule(v) === m
                body = Base.unwrap_unionall(v)
                body isa DataType && !isabstracttype(body) && body <: T && push!(found, v)
            end
        end
    end
    walk(root)
    return unique(found)
end

# Whether `a` and `b` hold the same state: equal bits for bits values, equal values
# for big numbers (a `BigInt`'s limbs live behind a pointer that differs between
# copies, as in a `MersenneTwister`), equal elements for arrays, the same entries for
# dictionaries (keys paired one to one and compared like values, `_same_entries`), and
# the same field values for every other object, followed recursively. An
# `ImmutableDict` is a chain of entries, newest first, in which a newer entry shadows an
# older one with an equal key, so its entries are compared in that order
# (`_same_chain`): the same entries in another order can look up other values. Mutable
# objects are matched one to one in both directions, so aliasing is part of the
# state: two objects of `a` cannot both correspond to one object of `b`, nor the
# reverse. Tasks, locks and IO streams are not model state and are not compared.
#
# Dictionary entries that look alike can be paired more than one way, and an alias met
# later in the state can rule out the pairing first chosen. The comparison therefore
# searches: when it fails, it starts again from the beginning with the next pairing at
# the last entry that had another candidate (`_Pairings`), and the states are the same
# when one complete pass succeeds. A pass costs one comparison, and states without
# look-alike entries need one pass. After `_MAX_REPAIRINGS` further passes it gives up,
# warns, and reports the states as different.
function _same_state(a, b)
    pairings = _Pairings()
    for _ in 0:_MAX_REPAIRINGS
        _same_state(a, b, IdDict{Any, Any}(), IdDict{Any, Any}(), pairings) && return true
        _next_pairings!(pairings) || return false
    end
    @warn "The state comparison stopped after $(_MAX_REPAIRINGS) alternative pairings of dictionary " *
        "entries that look alike, and reports the states as different. Dictionaries with many keys and " *
        "values of equal content need more; keys of numbers, symbols or strings need none."
    return false
end

function _same_state(a, b, seen::IdDict{Any, Any}, back::IdDict{Any, Any}, pairings)
    typeof(a) === typeof(b) || return false
    T = typeof(a)
    isbitstype(T) && return a === b
    (a isa Module || a isa Function && !_has_captures(a) || a isa Type || a isa Symbol) && return a === b
    a isa AbstractString && return a == b
    a isa Union{BigInt, BigFloat} && return isequal(a, b)
    (a isa Task || a isa Base.AbstractLock || a isa IO || a isa Threads.Condition) && return true
    if ismutable(a)
        haskey(seen, a) && return seen[a] === b
        haskey(back, b) && return false
        seen[a] = b
        back[b] = a
    end
    a isa Base.ImmutableDict && return _same_chain(a, b, seen, back, pairings)
    a isa AbstractDict && return _same_entries(a, b, seen, back, pairings)
    if a isa Array
        size(a) == size(b) || return false
        for i in eachindex(a)
            isassigned(a, i) == isassigned(b, i) || return false
            isassigned(a, i) && !_same_state(a[i], b[i], seen, back, pairings) && return false
        end
        return true
    end
    for i in 1:fieldcount(T)
        isdefined(a, i) == isdefined(b, i) || return false
        isdefined(a, i) && !_same_state(getfield(a, i), getfield(b, i), seen, back, pairings) && return false
    end
    return true
end

_has_captures(f::Function) = fieldcount(typeof(f)) > 0

# How many further passes a comparison may make before it gives up (see `_same_state`).
const _MAX_REPAIRINGS = 1000

# The search of one comparison. A choice point is an entry with two or more candidates,
# numbered in the order a pass meets it, which is the same in every pass until the
# point whose choice changed. `forced[p]` is the candidate the next pass starts from at
# point `p`; `taken[p]` and `more[p]` are the candidate this pass took there and
# whether a later candidate could also match. Inside a trial (`trial`), choices are
# guesses, neither forced nor recorded, and `guessed` notes that one was made.
mutable struct _Pairings
    forced::Vector{Int}
    taken::Vector{Int}
    more::Vector{Bool}
    trial::Bool
    guessed::Bool
end
_Pairings() = _Pairings(Int[], Int[], Bool[], false, false)

# Prepares the next pass after a failed one: the last point with another candidate
# moves on to it, the points before it keep their choices, and later points start
# afresh. False when no point has another candidate.
function _next_pairings!(pairings::_Pairings)
    p = findlast(pairings.more)
    p === nothing && return false
    pairings.forced = push!(pairings.taken[1:(p - 1)], pairings.taken[p] + 1)
    empty!(pairings.taken)
    empty!(pairings.more)
    return true
end

# Keys whose identity carries no state (bits values, symbols, strings). In a `Dict` or
# an `IdDict` such a key is found by lookup and then compared with the stored key the
# lookup gives (`_looked_up_key`); other dictionaries pair all their entries.
_plain_key(k) = isbitstype(typeof(k)) || k isa Symbol || k isa String

# The stored key of `b` that `k` looks up to, or a `_NoKey()` when there is none. A
# `Dict` finds keys by `isequal`, so it can hold another value (`1.0` for `1`), and
# `getkey` gives it; an `IdDict` finds keys by identity, so the stored key is `k`
# itself. Other dictionaries do not say which key they hold, so they give none: their
# lookup would hide a stored `1.0` behind a query `1`. Their entries are paired as they
# are iterated, each key held once (an `ImmutableDict`, which can hold a key more than
# once, is compared as its chain instead, `_same_chain`).
struct _NoKey end
_looked_up_key(b::Dict, k) = getkey(b, k, _NoKey())
_looked_up_key(b::IdDict, k) = haskey(b, k) ? k : _NoKey()
_looked_up_key(b::AbstractDict, k) = _NoKey()

# A hash of the state of `x` that does not depend on identity: states `_same_state`
# finds the same hash alike. It follows at most `nodes[]` values, in a fixed order, and
# hashes a dictionary by its type and length only, since the order of its entries can
# differ between equal dictionaries.
function _shape(x, h::UInt, nodes::Base.RefValue{Int})
    T = typeof(x)
    h = hash(objectid(T), h)
    nodes[] > 0 || return h
    nodes[] -= 1
    (isbitstype(T) || x isa Module || x isa Function && !_has_captures(x) || x isa Type || x isa Symbol) &&
        return hash(objectid(x), h)
    x isa Union{AbstractString, BigInt, BigFloat} && return hash(x, h)
    (x isa Task || x isa Base.AbstractLock || x isa IO || x isa Threads.Condition) && return h
    x isa AbstractDict && return hash(length(x), h)
    if x isa Array
        h = hash(size(x), h)
        for i in eachindex(x)
            nodes[] > 0 || break
            h = isassigned(x, i) ? _shape(x[i], h, nodes) : hash(:unassigned, h)
        end
        return h
    end
    for i in 1:fieldcount(T)
        nodes[] > 0 || break
        h = isdefined(x, i) ? _shape(getfield(x, i), h, nodes) : hash(:undefined, h)
    end
    return h
end

_entry_shape(k, v) = _shape(v, _shape(k, zero(UInt), Ref(32)), Ref(32))

# Whether two dictionaries hold the same entries. Keys are paired one to one and
# compared like any other state, through the same alias maps as the values, so a key
# that is shared with another part of the state must stay shared, and a copied
# identity dictionary (an `IdDict`, whose keys are found by identity) still matches.
# Plain keys of a `Dict` or an `IdDict` are paired by lookup. Every other entry of `a`
# is paired with a stored entry of `b` of the same shape (`_entry_shape`): the one the
# alias maps already give its key, or else one chosen by `_paired_entry`.
function _same_entries(a::AbstractDict, b::AbstractDict, seen::IdDict{Any, Any}, back::IdDict{Any, Any},
        pairings::_Pairings)
    length(a) == length(b) || return false
    lookup = a isa Union{Dict, IdDict}
    others = Pair[]
    for (k, v) in a
        if lookup && _plain_key(k)
            _same_state(k, _looked_up_key(b, k), seen, back, pairings) || return false
            _same_state(v, b[k], seen, back, pairings) || return false
        else
            push!(others, k => v)
        end
    end
    entries = Pair[k => v for (k, v) in b if !(lookup && _plain_key(k))]
    length(entries) == length(others) || return false
    groups = Dict{UInt, Vector{Int}}()
    for (j, (k, v)) in enumerate(entries)
        push!(get!(Vector{Int}, groups, _entry_shape(k, v)), j)
    end
    for (k, v) in others
        group = get(groups, _entry_shape(k, v), Int[])
        if ismutable(k) && haskey(seen, k)
            # The key is already matched elsewhere in the state: only its partner will do.
            i = findfirst(j -> first(entries[j]) === seen[k], group)
            i === nothing && return false
            _same_state(v, last(entries[group[i]]), seen, back, pairings) || return false
        else
            i = _paired_entry(k, v, b, entries, group, seen, back, pairings)
            i === nothing && return false
        end
        deleteat!(group, i)
    end
    return true
end

# Whether two `ImmutableDict`s hold the same entries in the same order, newest first.
# The chains are compared entry by entry rather than field by field: the empty
# dictionary at the end of each leaves its value field unset.
function _same_chain(a::Base.ImmutableDict, b::Base.ImmutableDict, seen::IdDict{Any, Any}, back::IdDict{Any, Any},
        pairings::_Pairings)
    length(a) == length(b) || return false
    for ((ka, va), (kb, vb)) in zip(a, b)
        _same_state(ka, kb, seen, back, pairings) && _same_state(va, vb, seen, back, pairings) || return false
    end
    return true
end

# Pairs the entry `k => v` of `a` with one of the unpaired entries `group` of `b` of
# the same shape, the one `k` looks up to first, and returns its place in `group`, or
# `nothing`. With one candidate, it is compared. With more, this is a choice point:
# candidates are tried on copies of the alias maps (`_try_entry`), from the one the
# search forces, and the first that can match is taken. A sure match commits its
# copies; a candidate whose trial had to guess is compared again on the maps, with
# its own choices recorded.
function _paired_entry(k, v, b::AbstractDict, entries::Vector{Pair}, group::Vector{Int},
        seen::IdDict{Any, Any}, back::IdDict{Any, Any}, pairings::_Pairings)
    isempty(group) && return nothing
    order = collect(eachindex(group))
    looked_up = _looked_up_key(b, k)
    i = findfirst(j -> first(entries[j]) === looked_up, group)
    i === nothing || (deleteat!(order, i); pushfirst!(order, i))
    compare(i) = _same_state(k, first(entries[group[i]]), seen, back, pairings) &&
        _same_state(v, last(entries[group[i]]), seen, back, pairings)
    length(order) == 1 && return compare(only(order)) ? only(order) : nothing
    if pairings.trial
        pairings.guessed = true
        return compare(first(order)) ? first(order) : nothing
    end
    point = length(pairings.taken) + 1
    start = point <= length(pairings.forced) ? pairings.forced[point] : 1
    trial(t) = _try_entry(k, v, entries[group[order[t]]], seen, back, pairings)
    t = start
    outcome = nothing
    while t <= length(order)
        outcome = trial(t)
        first(outcome) === :mismatch || break
        t += 1
    end
    if t > length(order)
        push!(pairings.taken, t)
        push!(pairings.more, false)
        return nothing
    end
    push!(pairings.taken, t)
    push!(pairings.more, any(u -> first(trial(u)) !== :mismatch, (t + 1):length(order)))
    status, trial_seen, trial_back = outcome
    if status === :match
        merge!(empty!(seen), trial_seen)
        merge!(empty!(back), trial_back)
    else
        compare(order[t]) || return nothing
    end
    return order[t]
end

# Compares `k => v` with `entry` on copies of the alias maps, guessing at any choice
# inside: `:match` or `:mismatch` when no guess was needed, `:unsure` otherwise, with
# the copies.
function _try_entry(k, v, entry::Pair, seen::IdDict{Any, Any}, back::IdDict{Any, Any}, pairings::_Pairings)
    trial_seen, trial_back = copy(seen), copy(back)
    pairings.trial, pairings.guessed = true, false
    same = _same_state(k, first(entry), trial_seen, trial_back, pairings) &&
        _same_state(v, last(entry), trial_seen, trial_back, pairings)
    status = pairings.guessed ? :unsure : same ? :match : :mismatch
    pairings.trial, pairings.guessed = false, false
    return (status, trial_seen, trial_back)
end

# Every mutable object reachable from `x` (arrays, dictionaries, mutable structs),
# by identity, skipping what `_same_state` does not compare.
function _reachable_mutables(x, found::IdDict{Any, Nothing}=IdDict{Any, Nothing}())
    T = typeof(x)
    isbitstype(T) && return found
    (x isa Module || x isa Type || x isa Symbol || x isa AbstractString || x isa Task ||
        x isa Base.AbstractLock || x isa IO || x isa Threads.Condition) && return found
    if ismutable(x)
        haskey(found, x) && return found
        found[x] = nothing
    end
    if x isa AbstractDict
        for (k, v) in x
            _reachable_mutables(k, found)
            _reachable_mutables(v, found)
        end
    elseif x isa Array
        for i in eachindex(x)
            isassigned(x, i) && _reachable_mutables(x[i], found)
        end
    else
        for i in 1:fieldcount(T)
            isdefined(x, i) && _reachable_mutables(getfield(x, i), found)
        end
    end
    return found
end

# Mutable objects reachable from both `a` and `b`.
function _shared_mutables(a, b)
    from_a = _reachable_mutables(a)
    return [object for object in keys(_reachable_mutables(b)) if haskey(from_a, object)]
end
