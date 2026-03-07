module AbstractTypes

export AbstractForceTorqueModel, AbstractPlanet, AbstractDensityModel, AbstractThermalModel, AbstractThrusterModel, AbstractControlEffectorModel, AbstractEphemeridesModel

# The "contract" for all force/torque models
abstract type AbstractForceTorqueModel end
abstract type AbstractControlEffectorModel end
abstract type AbstractPlanet end
abstract type AbstractDensityModel end
abstract type AbstractThermalModel end
abstract type AbstractEphemeridesModel end
abstract type AbstractThrusterModel end
abstract type AbstractGuidanceModel end
end # module AbstractTypes
