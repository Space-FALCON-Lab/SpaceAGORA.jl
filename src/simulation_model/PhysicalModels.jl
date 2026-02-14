module EnvironmentModels
    using ..AbstractTypes: AbstractPlanet, AbstractDensityModel, AbstractThermalModel
    using ..Analysis
    using ..Kinematics
    using ..Effectors
    using ..SimConfig: InitialTime

    export NoAtmosphereModel, ExponentialAtmosphereModel, GRAMAtmosphereModel
    export getDensity
    
    export MaxwellianHeat, getHeatRate

    include("../physical_models/Density_models.jl")
    include("../physical_models/Thermal_models.jl")
end # module PhysicalModels