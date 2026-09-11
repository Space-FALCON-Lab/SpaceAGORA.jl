module Constellations

using ..SpacecraftModels: SpacecraftModel

export constellation_struct, build_constellation
export activate_link!, deactivate_link!, reset_active_links!

# ── Constellation struct ───────────────────────────────────────────────────────

mutable struct constellation_struct
    spacecraft::Vector{SpacecraftModel}
    possible_links::Vector{Tuple{Int, Int}}  # all pairs that geometry could ever allow
    active_links::Vector{Tuple{Int, Int}}    # subset currently firing; mutated at runtime
end

# ── Constructors ───────────────────────────────────────────────────────────────

# Builds a constellation from caller-supplied spacecraft — one SpacecraftModel per satellite,
# each already carrying its own orbit (in initial_condition), mass, and bus properties.
# Inputs: spacecraft_model_list, a list/array/tuple of SpacecraftModel (any order).
# Output: constellation_struct; ids are (re)assigned 1..N by position, and possible_links is every
# directed pair between distinct satellites. For a non-default topology, construct
# constellation_struct(spacecraft, possible_links, active_links) directly.
function build_constellation(spacecraft_model_list)
    spacecraft = collect(spacecraft_model_list)
    for (id, sc) in enumerate(spacecraft)
        sc.id = id
    end
    possible = _all_directed_pairs(spacecraft)
    return constellation_struct(spacecraft, possible, copy(possible))
end

# ── Link update helpers ────────────────────────────────────────────────────────

# Adds (emitter_id, receiver_id) to active_links if it exists in possible_links; errors otherwise.
# Inputs: constellation_struct (mutated), emitter/receiver satellite IDs (Int).
# Output: nothing; mutates c.active_links in place.
function activate_link!(c::constellation_struct, emitter_id::Int, receiver_id::Int)
    pair = (emitter_id, receiver_id)                                                          # build the link tuple
    pair in c.possible_links || throw(ArgumentError("($emitter_id, $receiver_id) is not a possible link."))  # reject invalid link
    pair in c.active_links   || push!(c.active_links, pair)                                    # add only if not already active
end

# Removes (emitter_id, receiver_id) from active_links; silently a no-op if the pair is already absent.
# Inputs: constellation_struct (mutated), emitter/receiver satellite IDs (Int).
# Output: nothing; mutates c.active_links in place.
function deactivate_link!(c::constellation_struct, emitter_id::Int, receiver_id::Int)
    filter!(p -> p != (emitter_id, receiver_id), c.active_links)  # drop the pair if present
end

# Restores active_links to the full set of possible_links, undoing any runtime activations/deactivations.
# Input: constellation_struct (mutated).
# Output: nothing; active_links is a fresh copy of possible_links after the call.
function reset_active_links!(c::constellation_struct)
    resize!(c.active_links, length(c.possible_links))  # match possible_links length
    copyto!(c.active_links, c.possible_links)          # copy all pairs back in
end

# ── Private helpers ────────────────────────────────────────────────────────────

# Computes every directed (emitter_id, receiver_id) pair between distinct satellites — the topology
# used for constellations built from a flat list of elements.
function _all_directed_pairs(spacecraft::Vector{SpacecraftModel})
    ids = [sc.id for sc in spacecraft]
    return [(i, j) for i in ids for j in ids if i != j]
end

end # module Constellations
