module EnvironmentModels
    using ..AbstractTypes: AbstractPlanet, AbstractDensityModel, AbstractThermalModel
    using ..Analysis
    using ..Kinematics
    using ..Effectors
    using ..SimConfig: InitialTime
    using Reexport

    export NoAtmosphereModel, ExponentialAtmosphereModel, GRAMAtmosphereModel, GRAMAtmosphereModelSurrogate
    export getDensity, getDensityBatch!, precompute_gram_static_grids!, clear_gram_static_grid_cache!
    
    export MaxwellianHeat, getHeatRate

    include(joinpath(@__DIR__, "..", "environment", "atmosphere", "density_models.jl"))
    # @reexport DensityModels
    include(joinpath(@__DIR__, "..", "environment", "thermal", "thermal_models.jl"))
end # module PhysicalModels
