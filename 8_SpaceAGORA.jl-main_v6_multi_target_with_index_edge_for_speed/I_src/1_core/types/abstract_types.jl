module AbstractTypes

export AbstractForceTorqueModel, AbstractPlanet, AbstractDensityModel, AbstractThermalModel, AbstractThrusterModel, AbstractControlEffectorModel, AbstractEphemeridesModel, AbstractGuidanceModel

"""
    AbstractForceTorqueModel

Stable extension interface for custom force/torque model implementations that
contribute translational and/or rotational dynamics terms.
"""
abstract type AbstractForceTorqueModel end

"""
    AbstractControlEffectorModel

Stable extension interface for control/execution effectors that produce runtime
commands or actuation-side control behavior.
"""
abstract type AbstractControlEffectorModel end

"""
    AbstractPlanet

Stable extension interface for planet/body definitions used by environment,
ephemerides, and dynamics models.
"""
abstract type AbstractPlanet end

"""
    AbstractDensityModel

Stable extension interface for atmosphere density models, including no-atmosphere,
analytic, GRAM-backed, and user-defined models.
"""
abstract type AbstractDensityModel end

"""
    AbstractThermalModel

Stable extension interface for vehicle thermal/heat-rate models.
"""
abstract type AbstractThermalModel end

"""
    AbstractEphemeridesModel

Stable extension interface for ephemerides/frame backends such as SPICE-backed or
simple analytic models.
"""
abstract type AbstractEphemeridesModel end

"""
    AbstractThrusterModel

Stable extension interface for thruster hardware models used by actuator and
control layers.
"""
abstract type AbstractThrusterModel end

"""
    AbstractGuidanceModel

Stable extension interface for guidance strategy models that produce guidance-side
commands or trajectory directives.
"""
abstract type AbstractGuidanceModel end
end # module AbstractTypes
