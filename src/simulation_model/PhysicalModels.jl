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

    include("../physical_models/Density_models.jl")
    # @reexport DensityModels
    include("../physical_models/Thermal_models.jl")
end # module PhysicalModels
