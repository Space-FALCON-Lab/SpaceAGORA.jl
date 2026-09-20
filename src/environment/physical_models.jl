module EnvironmentModels
    using ..AbstractTypes: AbstractPlanet, AbstractDensityModel
    using ..Kinematics
    using ..SimConfig: InitialTime
    using Reexport

    export NoAtmosphereModel, ExponentialAtmosphereModel, PiecewiseExponentialAtmosphereModel, TabulatedFlightAtmosphereModel, TimeTabulatedAtmosphereModel
    export NRLMSISE00AtmosphereModel, init_nrlmsise_space_indices!
    export GRAMAtmosphereModel, GRAMAtmosphereModelSurrogate, GRAMGridAtmosphereModel, ConstantDensityModel
    export with_density_model_epoch
    export getDensity, getDensityBatch!, precompute_gram_static_grids!, clear_gram_static_grid_cache!

    include(joinpath(@__DIR__, "..", "environment", "atmosphere", "density_models.jl"))
    include(joinpath(@__DIR__, "atmosphere", "density_epoch.jl"))
    include(joinpath(@__DIR__, "atmosphere", "gram_grid_atmosphere_model.jl"))
    include(joinpath(@__DIR__, "atmosphere", "surrogate_presets.jl"))
end # module EnvironmentModels
