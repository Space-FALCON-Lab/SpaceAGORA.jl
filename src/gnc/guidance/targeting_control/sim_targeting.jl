# Compatibility wrapper: canonical path forwarding to legacy implementation.
include(joinpath(let
    p = @__DIR__
    while basename(p) != "src"
        nextp = dirname(p)
        nextp == p && error("Could not locate src root from $(@__DIR__)")
        p = nextp
    end
    p
end, "control/targeting_control/sim_targeting.jl"))
