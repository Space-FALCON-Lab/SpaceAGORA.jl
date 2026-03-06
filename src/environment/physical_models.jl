module EnvironmentModels
    using ..AbstractTypes: AbstractPlanet, AbstractDensityModel
    using ..Analysis
    using ..Kinematics
    using ..Effectors
    using ..SimConfig: InitialTime
    using Reexport

    export NoAtmosphereModel, ExponentialAtmosphereModel, GRAMAtmosphereModel, GRAMAtmosphereModelSurrogate
    export getDensity, getDensityBatch!, precompute_gram_static_grids!, clear_gram_static_grid_cache!

    include(joinpath(@__DIR__, "..", "environment", "atmosphere", "density_models.jl"))
end # module EnvironmentModels
