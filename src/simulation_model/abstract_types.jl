module AbstractTypes

export AbstractForceTorqueModel, Planet

# The "contract" for all force/torque models
abstract type AbstractForceTorqueModel end
abstract type Planet end

end # module AbstractTypes