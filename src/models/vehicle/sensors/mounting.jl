"""
    SensorMounting

Mounting sensors on spacecraft links and reading their measurements.
"""
module SensorMounting

using SHA
using Random: Xoshiro
using ..AbstractTypes: AbstractSensor
using ..SpacecraftModels: Link, SpacecraftModel
using ..Sensors: sensor_dimension, sensor_requirements, sensor_period_s

export add_sensor!, mounted_sensors, sensor_measurement, sensor_measurement_time, reset_sensors!

const _RUN_STATE_FIELDS = (:period_s, :seed, :rng, :latest, :latest_time_s)

"""
    add_sensor!(model::SpacecraftModel, link::Link, sensor) -> sensor

Mount `sensor` on one of `model`'s links; the link's frame is the sensor's
frame (a sensor on a hinged panel turns with it). Sensors are numbered in
mounting order, links in order, for [`sensor_measurement`](@ref) and the
result columns. Returns the sensor.
"""
function add_sensor!(model::SpacecraftModel, link::Link, sensor::AbstractSensor)
    any(l -> l === link, model.links) || throw(ArgumentError(
        "Sensors mount on one of the spacecraft's own links; this link is not in its link list."))
    absent = [f for f in _RUN_STATE_FIELDS if !hasfield(typeof(sensor), f)]
    isempty(absent) || throw(ArgumentError(
        "A sensor kind needs the fields $(join(_RUN_STATE_FIELDS, ", ")); $(typeof(sensor)) lacks $(join(absent, ", "))."))
    ismutabletype(typeof(sensor)) || throw(ArgumentError("A sensor kind must be a mutable struct: the engine fills its run state."))
    req = sensor_requirements(sensor)
    req.atmosphere && throw(ArgumentError(
        "Sensors may ask for the magnetic field, the Sun and the planet frame only; this sensor asks for the atmosphere."))
    sensor_dimension(sensor) >= 1 || throw(ArgumentError("A sensor's sensor_dimension must be positive."))
    _check_period(sensor)
    any(s -> s === sensor, (s for l in model.links for s in l.sensors)) &&
        throw(ArgumentError("This sensor is already mounted on the spacecraft."))
    push!(link.sensors, sensor)
    return sensor
end

"""
    mounted_sensors(model) -> Vector{NamedTuple{(:link, :sensor)}}

Every sensor mounted on `model`, with its link, in mounting order.
"""
function mounted_sensors(model::SpacecraftModel)
    out = NamedTuple{(:link, :sensor), Tuple{Link, AbstractSensor}}[]
    for link in model.links, s in link.sensors
        push!(out, (link=link, sensor=s))
    end
    return out
end

"""
    sensor_measurement(sensor) -> measurement
    sensor_measurement(model::SpacecraftModel, k::Int) -> measurement

The latest measurement of `sensor`, or of the `k`-th sensor mounted on
`model` (NaN before its first sample). A navigation effector reads it in
`update_navigation!`, from the spacecraft of the run
(`p.args.dynamics_model.spacecraft[i]`).
"""
sensor_measurement(sensor::AbstractSensor) = sensor.latest
sensor_measurement(model::SpacecraftModel, k::Int) = sensor_measurement(mounted_sensors(model)[k].sensor)

"""
    sensor_measurement_time(sensor) -> Float64
    sensor_measurement_time(model::SpacecraftModel, k::Int) -> Float64

The time of the latest measurement, s into the run (NaN before the first
sample).
"""
sensor_measurement_time(sensor::AbstractSensor)::Float64 = sensor.latest_time_s
sensor_measurement_time(model::SpacecraftModel, k::Int)::Float64 = sensor_measurement_time(mounted_sensors(model)[k].sensor)

# A sensor's sampling period must be finite and positive.
function _check_period(sensor::AbstractSensor)
    period = sensor_period_s(sensor)
    (period isa Real && isfinite(period) && period > 0) || throw(ArgumentError(
        "A sensor's sampling period must be finite and positive; $(typeof(sensor)) has $(repr(period))."))
    return Float64(period)
end

# The random stream of a sensor: seeded from the SHA-256 digest of the
# sensor's seed, the spacecraft's ID and the sensor's mounting position, so
# runs with the same seeds repeat and each sensor has its own stream.
function _sensor_rng(seed::Integer, spacecraft_id::Integer, position::Integer)
    words = reinterpret(UInt64, sha256("SpaceAGORA/sensor/v1/$(seed)/$(spacecraft_id)/$(position)"))
    return Xoshiro(words[1], words[2], words[3], words[4])
end

"""
    reset_sensors!(model) -> model

At the start of a run: each sensor's random stream from its seed, its
measurement unset, a gyro's bias at its initial value.
"""
function reset_sensors!(model::SpacecraftModel)
    for (k, entry) in enumerate(mounted_sensors(model))
        s = entry.sensor
        s.rng = _sensor_rng(s.seed, model.id, k)
        s.latest = map(_ -> NaN, s.latest)
        s.latest_time_s = NaN
        hasfield(typeof(s), :initial_bias) && hasfield(typeof(s), :bias) && (s.bias = s.initial_bias)
    end
    return model
end

end # module SensorMounting
