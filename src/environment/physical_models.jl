module EnvironmentModels
    using ..AbstractTypes: AbstractPlanet, AbstractDensityModel
    using ..Kinematics
    using ..SimConfig: InitialTime
    using Reexport

    export NoAtmosphereModel, ExponentialAtmosphereModel, PiecewiseExponentialAtmosphereModel, TabulatedFlightAtmosphereModel
    export NRLMSISE00AtmosphereModel, init_nrlmsise_space_indices!
    export GRAMAtmosphereModel, GRAMAtmosphereModelSurrogate, ConstantDensityModel
    export getDensity, getDensityBatch!, precompute_gram_static_grids!, clear_gram_static_grid_cache!

    include(joinpath(@__DIR__, "..", "environment", "atmosphere", "density_models.jl"))
end # module EnvironmentModels
