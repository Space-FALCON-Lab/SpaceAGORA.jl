"""
    Sensors

The sensor contract: a sensor kind is a mutable record subtyping
`AbstractSensor` with methods of `sensor_measure` and `sensor_dimension`;
mounted on a link with `add_sensor!`, it is sampled at its period from the
accepted state, and navigation effectors read its latest measurement with
`sensor_measurement`.
"""
module Sensors

using StaticArrays
using ..AbstractTypes: AbstractSensor
using ..Actuators: ActuatorContext, ActuatorRequirements

export SensorTruth, sensor_measure, sensor_dimension, sensor_requirements, sensor_period_s, initialize_sensor!
export sensor_needs_specific_force

"""
    SensorTruth

What a sensor measures from, at a sample: the inertial position and velocity
(m, m/s), the attitude quaternion (inertial to body, scalar last), the body
rate (rad/s), the mass (kg), the time (s into the run), and, each NaN unless
the sensor asks for it, the local magnetic field in the body frame (T), the
direction to the Sun in the body frame (a unit vector) with the lit fraction
of the Sun's disk (0 in umbra, 1 in full sunlight) and the planet's
equatorial radius (m), the specific force (the non-gravitational
acceleration, m/s², in the body frame) with the body's angular acceleration
(rad/s², body frame), and the position, velocity (planet-fixed frame) and
altitude (m) over the planet. Without attitude propagation the attitude is
the identity and the rate and angular acceleration zero.

`SensorTruth(pos_ii, vel_ii, q_ib, ω_body, mass_kg, t_s, magnetic_field_body)`
leaves everything after the field NaN; the ten-argument form (through
`specific_force_body`) leaves the rest NaN.
"""
struct SensorTruth
    pos_ii::SVector{3, Float64}
    vel_ii::SVector{3, Float64}
    q_ib::SVector{4, Float64}
    ω_body::SVector{3, Float64}
    mass_kg::Float64
    t_s::Float64
    magnetic_field_body::SVector{3, Float64}
    sun_direction_body::SVector{3, Float64}
    sun_lit_fraction::Float64
    specific_force_body::SVector{3, Float64}
    angular_acceleration_body::SVector{3, Float64}
    pos_pp::SVector{3, Float64}
    vel_pp::SVector{3, Float64}
    altitude_m::Float64
    planet_radius_m::Float64
end

const _NAN3_TRUTH = SVector{3, Float64}(NaN, NaN, NaN)
SensorTruth(pos_ii, vel_ii, q_ib, ω_body, mass_kg, t_s, magnetic_field_body) =
    SensorTruth(pos_ii, vel_ii, q_ib, ω_body, mass_kg, t_s, magnetic_field_body, _NAN3_TRUTH, NaN, _NAN3_TRUTH)
SensorTruth(pos_ii, vel_ii, q_ib, ω_body, mass_kg, t_s, magnetic_field_body, sun_direction_body, sun_lit_fraction, specific_force_body) =
    SensorTruth(pos_ii, vel_ii, q_ib, ω_body, mass_kg, t_s, magnetic_field_body, sun_direction_body, sun_lit_fraction,
                specific_force_body, _NAN3_TRUTH, _NAN3_TRUTH, _NAN3_TRUTH, NaN, NaN)

"""
    sensor_measure(sensor, truth::SensorTruth, ctx::ActuatorContext, rng) -> measurement

Required: one measurement (an `SVector` or a number) of `truth` by `sensor`,
mounted where `ctx` says (`ctx.link_to_body` maps the sensor's frame, its
link's, to the body frame). Draw any noise from `rng`, the sensor's own
random stream. Called once per sample, so a sensor may advance its own state
(a bias random walk) here.
"""
function sensor_measure end

"""
    initialize_sensor!(sensor, truth::SensorTruth, ctx::ActuatorContext, rng) -> nothing

Optional: reset a sensor's own run state at the start of a run. It is called
once per run, after the engine has reset the fields it owns (the random stream,
the latest measurement and, for a gyro, the bias) and just before the first
measurement. It receives that measurement's `truth` and `ctx` and the fresh
`rng`. It is not called again for later checkpoint segments. A sensor that keeps
more state than those fields, such as a filter, resets it here. Default: nothing.
"""
initialize_sensor!(sensor, truth, ctx, rng) = nothing  # untyped, so no sensor's method is ambiguous with it

"""
    sensor_dimension(sensor) -> Int

Required: the number of components of the sensor's measurement, which sets
its result columns.
"""
function sensor_dimension end

"""
    sensor_requirements(sensor) -> ActuatorRequirements

Optional: what the sensor needs sampled. Sensors may ask for the magnetic
field (`magnetic_field=true`), the Sun (`solar=true`: its direction, the lit
fraction and the planet's radius) and the planet frame (`planet_frame=true`:
the planet-fixed position, velocity and altitude); the atmosphere is refused
when it is mounted.
Default: nothing.
"""
@inline sensor_requirements(::Any) = ActuatorRequirements()

"""
    sensor_needs_specific_force(sensor) -> Bool

Optional: whether the sensor measures the specific force, which the engine
then computes at its samples: the vehicle's acceleration less the
acceleration of the run's gravitational effectors (see
`is_gravitational`). Default: false.
"""
@inline sensor_needs_specific_force(::Any) = false

"""
    sensor_period_s(sensor) -> Float64

Optional: the sampling period, s. Default: the sensor's `period_s` field.

Besides its parameters, a sensor kind holds its run state in four fields the
engine fills: `seed` (an integer the user sets, which selects the random
stream), `rng` (the stream), `latest` (the latest measurement: a vector of
`sensor_dimension` components, or a number for a one-component sensor; NaN
before the first sample) and `latest_time_s` (its time, s into the run). The
period must be finite and positive.
"""
@inline sensor_period_s(sensor) = Float64(sensor.period_s)

end # module Sensors
