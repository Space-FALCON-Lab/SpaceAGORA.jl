module VehicleThermalModels
    using ..AbstractTypes: AbstractThermalModel, AbstractPlanet

    export MaxwellianHeat, getHeatRate

    include(joinpath(@__DIR__, "thermal_models.jl"))
end
