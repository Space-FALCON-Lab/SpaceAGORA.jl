module AbstractTypes

export AbstractTerrainModel, AbstractForceTorqueModel, AbstractPlanet, AbstractDensityModel, AbstractThermalModel, AbstractThrusterModel, AbstractControlEffectorModel, AbstractEphemeridesModel, AbstractGuidanceModel
export AbstractNavigationModel
export AbstractActuator, AbstractSensor

"""
    AbstractForceTorqueModel

Stable extension interface for custom force/torque model implementations that
contribute translational and/or rotational dynamics terms.
"""
abstract type AbstractForceTorqueModel end

"""
    AbstractControlEffectorModel

Stable extension interface for control/execution effectors that produce runtime
commands or actuation-side control behavior. A subtype implements
`update_control!(model, u, p, t, spacecraft_index)` and, unless its `actuation_mode`
is `:mounted`, `control_wrench(model, x::StateSample, p, t)`.
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
    AbstractActuator

Stable extension interface for a piece of hardware mounted on a spacecraft link
that acts on the vehicle when commanded: a thruster, and the kinds added after
it. A kind implements `actuator_wrench` and `actuator_limits`; the other
actuator methods have defaults. Mount it with `add_actuator!`.
"""
abstract type AbstractActuator end

"""
    AbstractSensor

Stable extension interface for a measuring instrument mounted on a spacecraft
link: a rate gyro, a star tracker, a magnetometer, and kinds users add. A kind
implements `sensor_measure` and `sensor_dimension`; mount it with
`add_sensor!`. Sensors observe the state; they do not act on the dynamics.
"""
abstract type AbstractSensor end

"""
    AbstractGuidanceModel

Stable extension interface for guidance strategy models that produce guidance-side
commands or trajectory directives. A subtype implements
`update_guidance!(model, u, p, t, spacecraft_index)`.
"""
abstract type AbstractGuidanceModel end

"""
    AbstractNavigationModel

Stable extension interface for navigation models: sensor and state-estimator
effectors listed in a `NavigationModel`. A subtype implements
`update_navigation!(model, u, p, t, spacecraft_index)`.
"""
abstract type AbstractNavigationModel end

"""
    AbstractTerrainModel

Extension interface for terrain height queries in metres above an explicitly
chosen reference sphere. Implement `terrain_height(model, lat_deg, lon_deg)`
with planetocentric latitude and east-positive longitude in degrees.
"""
abstract type AbstractTerrainModel end
end # module AbstractTypes
