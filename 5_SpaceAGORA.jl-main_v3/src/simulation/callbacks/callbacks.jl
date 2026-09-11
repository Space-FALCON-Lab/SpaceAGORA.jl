# Canonical aggregator: no behavior ownership.
module SimulationCallbacks
include(joinpath(@__DIR__, "registry.jl"))
include(joinpath(@__DIR__, "save_fields.jl"))
include(joinpath(@__DIR__, "density_callbacks.jl"))
include(joinpath(@__DIR__, "gram_track_cache.jl"))
include(joinpath(@__DIR__, "thermal_callbacks.jl"))
include(joinpath(@__DIR__, "event_callbacks.jl"))
include(joinpath(@__DIR__, "navigation_guidance_callbacks.jl"))
include(joinpath(@__DIR__, "control_callbacks.jl"))
end # module SimulationCallbacks
