module Constellations

using ..SpacecraftModels: SpacecraftModel

export constellation_struct, build_constellation
export activate_link!, deactivate_link!, reset_active_links!

# ── Constellation struct ───────────────────────────────────────────────────────

mutable struct constellation_struct
    spacecraft::Vector{SpacecraftModel}
    possible_links::Vector{Tuple{Int, Int}} # links currently in range; mutated by the scheduler each tick
    active_links::Vector{Tuple{Int, Int}}   # subset currently firing; mutated at runtime
end

# ── Constructors ───────────────────────────────────────────────────────────────

# Builds a constellation from caller-supplied spacecraft — one SpacecraftModel per satellite,
# each already carrying its own orbit (in initial_condition), mass, and bus properties.
# Inputs: spacecraft_model_list, a list/array/tuple of SpacecraftModel (any order).
# Output: constellation_struct; ids are (re)assigned 1..N by position. possible_links and active_links
# both start empty and are filled in by the scheduler once the simulation runs.
function build_constellation(spacecraft_model_list)
    spacecraft = collect(spacecraft_model_list)
    for (id, sc) in enumerate(spacecraft)
        sc.id = id
    end
    return constellation_struct(spacecraft, Tuple{Int, Int}[], Tuple{Int, Int}[])
end

end # module Constellations
