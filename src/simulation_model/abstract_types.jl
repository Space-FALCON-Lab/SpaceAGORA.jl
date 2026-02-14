module AbstractTypes

export AbstractForceTorqueModel, AbstractPlanet, AbstractDensityModel, AbstractThermalModel

# The "contract" for all force/torque models
abstract type AbstractForceTorqueModel end
abstract type AbstractPlanet end
abstract type AbstractDensityModel end
abstract type AbstractThermalModel end

end # module AbstractTypes