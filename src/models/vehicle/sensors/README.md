# Sensors

Measuring instruments mounted on spacecraft links.

| File | Holds |
|---|---|
| `sensors.jl` | `Sensors`: the contract (`SensorTruth`, `sensor_measure`, `sensor_dimension`, `sensor_requirements`, `sensor_period_s`, `sensor_needs_specific_force`) |
| `sensor_kinds.jl` | the built-in kinds (owned by `Components`): `RateGyro`, `StarTracker`, `Magnetometer`, `GNSSReceiver`, `SunSensor`, `Accelerometer` |
| `mounting.jl` | `SensorMounting`: `add_sensor!`, `mounted_sensors`, `sensor_measurement`, the per-sensor random streams |

The engine samples each sensor at the run's start and every period after it,
from the accepted state (`src/orchestration/simulation/callbacks/sensor_callbacks.jl`),
which also writes the sensor columns. Sensors do not act on the dynamics;
navigation effectors read their latest measurements.
