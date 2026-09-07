module Constellations

using ..SpacecraftModels: SpacecraftModel

export constellation_struct, build_constellation
export activate_link!, deactivate_link!, reset_active_links!

# ── Constellation struct ───────────────────────────────────────────────────────

mutable struct constellation_struct
    spacecraft::Vector{SpacecraftModel}
    possible_links::Vector{Tuple{Int, Int}}          # all pairs that geometry could ever allow
    active_links::Vector{Tuple{Int, Int}}            # subset currently firing; mutated at runtime
    previous_in_range_links::Vector{Tuple{Int, Int}} # in-range subset of possible_links as of the last scheduler tick
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
    return constellation_struct(spacecraft, possible, copy(possible), Tuple{Int, Int}[])
end

# ── Private helpers ────────────────────────────────────────────────────────────

# Computes every directed (emitter_id, receiver_id) pair between distinct satellites — the topology
# used for constellations built from a flat list of elements.
function _all_directed_pairs(spacecraft::Vector{SpacecraftModel})
    ids = [sc.id for sc in spacecraft]
    return [(i, j) for i in ids for j in ids if i != j]
end

end # module Constellations
