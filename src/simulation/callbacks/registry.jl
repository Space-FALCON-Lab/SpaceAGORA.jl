include(joinpath(@__DIR__, "save_fields.jl"))
include(joinpath(@__DIR__, "density_callbacks.jl"))
include(joinpath(@__DIR__, "thermal_callbacks.jl"))
include(joinpath(@__DIR__, "control_callbacks.jl"))
include(joinpath(@__DIR__, "navigation_guidance_callbacks.jl"))
include(joinpath(@__DIR__, "event_callbacks.jl"))
include(joinpath(@__DIR__, "gram_track_cache.jl"))

# Compatibility window source of truth.
include(joinpath(let
    p = @__DIR__
    while basename(p) != "src"
        nextp = dirname(p)
        nextp == p && error("Could not locate src root from $(@__DIR__)")
        p = nextp
    end
    p
end, "simulation_model/callbacks.jl"))
