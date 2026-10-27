# Sensors

Measuring instruments mounted on spacecraft links.

| File | Holds |
|---|---|
| `sensors.jl` | `Sensors`, the contract: `SensorTruth`, `sensor_measure`, `sensor_dimension`, `sensor_requirements`, `sensor_period_s`, `initialize_sensor!` and `sensor_needs_specific_force` |
| `sensor_kinds.jl` | the built-in kinds, owned by `Components`: `RateGyro`, `StarTracker`, `Magnetometer`, `GNSSReceiver`, `SunSensor` and `Accelerometer` |
| `mounting.jl` | `SensorMounting`: `add_sensor!`, `mounted_sensors`, `sensor_measurement`, `sensor_measurement_time`, `reset_sensors!` and the random stream of each sensor |

The engine samples each sensor from the accepted state at the start of the run
and then once per sensor period, in
`src/orchestration/simulation/callbacks/sensor_callbacks.jl`, which also writes
the sensor columns. Sensors do not act on the dynamics; navigation models read
their latest measurements.
