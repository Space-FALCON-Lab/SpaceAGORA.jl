# Full-feature test fixtures deliberately opt into the local companion.
# Independent installation tests do not use this helper.
let companion = normpath(joinpath(@__DIR__, "..", "..", "packages", "SpaceAGORAHYPR"))
    companion in LOAD_PATH || push!(LOAD_PATH, companion)
end
using SpaceAGORAHYPR
