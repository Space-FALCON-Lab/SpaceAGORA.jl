# Output Columns

Use this page to look up a column of the saved results table: what it means,
its units and its frame.

What to read next:

- [Simulation Outputs](outputs.md)
- [The Integrated State](integrated_state.md)

## CSV column reference

All per-spacecraft columns use the prefix `sc{N}_` where `N` is the
one-based position in the spacecraft list passed to `DynamicsModel`.
`SpacecraftModel.id` does not set this prefix. A single-spacecraft run
produces `sc1_*` columns even when its spacecraft has a different `id`.

### Time

| Column | Unit | Description |
|---|---|---|
| `time` | s | Elapsed simulation time from epoch |

### Position and velocity (inertial frame)

| Column | Unit | Description |
|---|---|---|
| `sc1_pos_1` | m | Inertial position X component |
| `sc1_pos_2` | m | Inertial position Y component |
| `sc1_pos_3` | m | Inertial position Z component |
| `sc1_vel_1` | m/s | Inertial velocity X component |
| `sc1_vel_2` | m/s | Inertial velocity Y component |
| `sc1_vel_3` | m/s | Inertial velocity Z component |

### Geodetic state

| Column | Unit | Description |
|---|---|---|
| `sc1_altitude` | m | Instantaneous geodetic altitude above the reference ellipsoid |
| `sc1_latitude_deg` | deg | Geodetic latitude |
| `sc1_longitude_deg` | deg | Longitude |

### Osculating orbit diagnostic

| Column | Unit | Description |
|---|---|---|
| `sc1_periapsis_altitude` | m | Osculating spherical periapsis altitude, `a * (1 - e) - planet.Rp_e` |

Here `a` and `e` are the semimajor axis and eccentricity derived from the current
inertial position and velocity. The periapsis column uses a reference sphere
with the planet's equatorial radius, `planet.Rp_e`. It describes the
instantaneous osculating Keplerian orbit; it is not a prediction of the minimum
geodetic altitude reached by the propagated trajectory. The two altitude columns
use different reference surfaces and need not agree at a periapsis event.

The aerobraking examples' separate `periapsis_events.csv` tables instead contain
geodetic estimates. Depending on the example and available solver output, these
are evaluated at event-located radial minima or selected from sampled geodetic
altitude minima.
Simulation apoapsis plots may likewise use interpolated sample crossings and an
initial sample. Those event and sampled values are distinct from
`sc1_periapsis_altitude`.

### Mass

| Column | Unit | Description |
|---|---|---|
| `sc1_mass` | kg | Total spacecraft mass (structure + remaining propellant) |

### Atmosphere and aerodynamics

| Column | Unit | Description |
|---|---|---|
| `sc1_wind_1` | m/s | Atmospheric wind vector X |
| `sc1_wind_2` | m/s | Atmospheric wind vector Y |
| `sc1_wind_3` | m/s | Atmospheric wind vector Z |
| `sc1_drag_1` | N | Aerodynamic drag force X (planet frame) |
| `sc1_drag_2` | N | Aerodynamic drag force Y |
| `sc1_drag_3` | N | Aerodynamic drag force Z |
| `sc1_lift_1` | N | Aerodynamic lift force X |
| `sc1_lift_2` | N | Aerodynamic lift force Y |
| `sc1_lift_3` | N | Aerodynamic lift force Z |
| `sc1_cross_1` | N | Aerodynamic cross force X |
| `sc1_cross_2` | N | Aerodynamic cross force Y |
| `sc1_cross_3` | N | Aerodynamic cross force Z |

### Thermal

| Column | Unit | Description |
|---|---|---|
| `sc1_heat_rate` | W/cm² | Largest of the per-link stagnation heat rates (the built-in Maxwellian model returns W/cm²) |
| `sc1_heat_load` | J/cm² | Largest of the per-link accumulated heat loads (each link integrates its own rate, without an area-unit conversion); not the sum over links |

### Attitude (orientation_sim only)

These columns are present only when `mission_configuration.orientation_sim = true`:

| Column | Unit | Description |
|---|---|---|
| `sc1_q_1` | — | Attitude quaternion component `x` (scalar-last convention `[x, y, z, w]`, inertial to body) |
| `sc1_q_2` | — | Attitude quaternion component `y` |
| `sc1_q_3` | — | Attitude quaternion component `z` |
| `sc1_q_4` | — | Attitude quaternion scalar component `w`; `[0, 0, 0, 1]` is the identity attitude |

#### Plot saved raw quaternion history

From the repository root, load the existing example helpers and explicitly
choose the spacecraft and destination:

```julia
include("examples/common.jl")
include("examples/support/aerobraking_mission_plot_utils.jl")
quaternion_history = _save_quaternion_history_plot(
    "results/my_run/simulation_results.csv";
    spacecraft_index=1,
    output_path="results/my_run/raw_quaternion.png",
)
```

`spacecraft_index` is the positive, one-based position in the spacecraft list,
not `SpacecraftModel.id`. The helper requires `time` (seconds) and all four
`scN_q_1` through `scN_q_4` columns for that spacecraft. Missing columns, empty
tables, and missing, nonnumeric or nonfinite selected values produce an error;
columns for other spacecraft are not inspected. Default quaternion save fields
are available only with `orientation_sim=true`.

The four traces retain the saved signs, component values and row order,
including non-unit and zero quaternions. Time is shown in minutes (`time / 60`);
single rows, repeated times and nonmonotone times are accepted. No normalization,
sign adjustment, interpolation or attitude conversion is applied. The component
labels are `q1`–`q4`, colored red, green, blue and orange. The result contains
`time_min`, `components`, `labels`, `colors`, `spacecraft_index`, `plot` and
`output_path`, so the prepared data can be inspected without inferring values
from the image. For an existing `DataFrame`,
`_quaternion_history_samples(df; spacecraft_index=1)` prepares the same data
without plotting or writing a file.

The call saves only the explicitly requested image and does not open a plot
window. Supply a supported filename extension (for example `.png` or `.pdf`)
and an existing destination directory; the input CSV cannot be the destination.
This opt-in helper does not run a simulation or change any existing example plots.

### Body angular-rate history (explicit save fields)

Body rates are not included by `default_save_fields`, even for orientation runs.
For an existing orientation-enabled configuration `args`, explicitly record the
body angular velocity state in radians per second with this save-field recipe:

```julia
using SpaceAGORA
using StaticArrays
SC = SpaceAGORA.SimulationModel.SimulationCallbacks
body_rate_field = SC.SaveField(
    :body_rate_rad_s,
    (u, t, integrator) -> [SVector{3, Float64}(sc.ω) for sc in u.sc];
    per_satellite=true,
    column_prefix="body_rate",
)
fields = [SC.default_save_fields(args)..., body_rate_field]
# Supply fields to your already configured run:
# run_simulation(args; save_fields=fields)
```

This records `scN_body_rate_1..3` in body-axis component order. `N` is the
spacecraft's one-based list position, not its model ID. The CYGNSS slew demo
already records the same schema for its single spacecraft. Existing custom
`:omega` save fields have a different schema and are not inferred as aliases.

To plot an existing CSV, use the example helpers explicitly:

```julia
include("examples/common.jl")
include("examples/support/aerobraking_mission_plot_utils.jl")
body_rate_history = _save_body_rate_history_plot(
    "results/my_run/simulation_results.csv";
    spacecraft_index=1,
    output_path="results/my_run/body_rate.png",
)
```

The helper requires `time` in seconds and all three selected body-rate columns
in rad/s. It plots time in minutes and angular velocity in deg/min, using the
legacy conversion `rad2deg.(rates) .* 60` and labels `ω1`, `ω2`, `ω3` in red,
green and blue. It preserves signs and row order, including repeated or
nonmonotone times; it neither derives rates from quaternions nor rotates frames.

The returned `rates_rad_s` contains copies of the saved values; `components`
contains the displayed deg/min values. The other fields are `time_min`,
`labels`, `colors`, `spacecraft_index`, `plot` and `output_path`.
`_body_rate_history_samples(df; spacecraft_index=1)` prepares the arrays without
rendering. Empty input, missing columns, invalid selected values or nonfinite
converted values produce an error. Other spacecraft's columns are not inspected.
Supply a filename extension and an existing destination directory. The input
CSV, including filesystem aliases, cannot be the destination. No plot window
or simulation is started by the plotting helper.

### Visualization fields

These columns are added automatically when `simulation_settings.save_visualization_scene = true`
and the run can supply them; see the [sidecar section](outputs.md#Visualization-scene-sidecar-(opt-in))
for the conditions:

| Column | Unit | Description |
|---|---|---|
| `sc1_link_pose_1` .. `sc1_link_pose_7n` | m, — | `[rx, ry, rz, qx, qy, qz, qw]` per non-root link, relative to the root bus, in the sidecar's link order |
| `sc1_thruster_level_1` .. `sc1_thruster_level_k` | — | Firing level of each thruster, 0 to 1, in the sidecar's thruster order |
| `sc1_arm_pose_1` .. `sc1_arm_pose_7m` | m, — | `[rx, ry, rz, qx, qy, qz, qw]` per robot-arm link: centre of mass relative to the spacecraft (inertial) and inertial attitude |
| `sc1_density` | kg/m³ | Atmospheric density at the spacecraft |
| `sun_dir_1` .. `sun_dir_3` | — | Unit vector from the planet's centre to the Sun, inertial frame, one per row (not per spacecraft) |

### Plume interaction

These columns are written by the plume-surface effector
([Plume interaction](plume_interaction.md)) for every spacecraft in the run,
when using the default save-field list, whether or not the visualization sidecar
is enabled. An explicit list must include them to write them:

| Column | Unit | Description |
|---|---|---|
| `sc1_plume_height_m` | m | Slant height of the vehicle reference point above the terrain along the engine axis (radial clearance divided by the downward exhaust cosine); infinite when the exhaust does not point toward the surface |
| `sc1_plume_pressure_pa` | Pa | Peak surface pressure under the plume |
| `sc1_plume_shear_pa` | Pa | Peak wall shear stress under the plume |
| `sc1_plume_erosion_kg_s` | kg/s | Regolith mass erosion rate |
| `sc1_plume_eroded_kg` | kg | Eroded mass integrated over accepted solver steps (trapezoidal), reset at the start of a run |
| `sc1_plume_ejecta_mps` | m/s | Bounded speed at which eroded grains leave the surface |
| `sc1_plume_ground_effect_n` | N | Modeled thrust augmentation in ground effect |

The plume quantities are diagnostics of an illustrative closure, as described on
the plume page; they are not a calibrated flight reconstruction.

### Mounted actuators

These columns are written with the default save-field list when a control
effector of the run is in `:mounted` actuation (see
[Extensibility](../extensibility.md#Add-an-actuator)), for every spacecraft
that carries mounted actuators. An explicit list must include the
`actuator_level` and `actuator_power` fields to write them:

| Column | Unit | Description |
|---|---|---|
| `sc1_actuator_level_1` .. `sc1_actuator_level_k` | the kind's | The level each mounted actuator applies, in mounting order: its command within its limits. A thruster with firing dynamics reports the average thrust its last control tick scheduled, a wheel at its momentum limit reports zero when the command would store more, a hinge or thruster gimbal reports its current angle, and a control moment gyro its gimbal rate (zero at a gimbal stop the command would pass) |
| `sc1_actuator_power` | W | Electrical power the mounted actuators draw at those levels, the sum of what each kind reports through `actuator_consumption` |

The levels come from the first control effector in `:mounted` actuation that
holds levels for the spacecraft. A row saved before that effector's first
control tick reports the levels it held then. Among the built-in kinds, laser
terminals draw `level / efficiency`; reaction wheels draw their base power and
their motor's power from the integrated wheel momentum, when they are given a
spin inertia; torque rods draw their power at full dipole times the square of
the dipole fraction. Hinges and thrusters report none.

Runs with control moment gyros also write their gimbal angles, from the
integrated state:

| Column | Unit | Description |
|---|---|---|
| `sc1_cmg_angle_1` .. `sc1_cmg_angle_n` | rad | Gimbal angle of each control moment gyro, in mounting order |

### Power budget

These columns are written with the default save-field list when a spacecraft
of the run carries power-system parts (solar arrays, a battery or constant
loads; see [Powering a Spacecraft](power_budget.md)). They are NaN for a
spacecraft without such parts, and the battery columns are NaN for one without
a battery:

| Column | Unit | Description |
|---|---|---|
| `sc1_battery_wh` | Wh | Energy stored in the battery (an integrated state) |
| `sc1_battery_soc` | — | State of charge, the stored energy over the capacity |
| `sc1_power_generated_w` | W | Power the solar arrays generate: lit cells facing the Sun, through the eclipse model |
| `sc1_power_load_w` | W | Load on the bus: constant loads plus the mounted actuators, with any shedding applied. With more than one control effector in `:mounted` actuation, the actuators each one commands, as the dynamics apply each one's |
| `sc1_power_shed_w` | W | Actuator power shed because the battery is empty |
| `sc1_power_net_w` | W | Generation minus load; negative while the battery is discharging, or while it is empty and the load exceeds the generation |

### Actuator faults

Written with the default save-field list when a spacecraft carries actuator
faults (`add_fault!`, see [Extensibility](../extensibility.md#Add-an-actuator)):

| Column | Unit | Description |
|---|---|---|
| `sc1_fault_active_1` .. `sc1_fault_active_m` | — | 1 while each fault (in the order it was added) is active, 0 otherwise |

The `actuator_level` columns show the levels as the faults leave them.

### Sensors

Written with the default save-field list when a spacecraft carries sensors
(`add_sensor!`, see [Extensibility](../extensibility.md#Add-a-sensor)):

| Column | Unit | Description |
|---|---|---|
| `sc1_sensor_1` .. `sc1_sensor_n` | the sensor's | The latest measurement of every sensor, sensors in mounting order and their components in order: a rate gyro's three rates (rad/s) and a magnetometer's three field components (T), in the sensor's frame; a star tracker's estimate of the body's attitude (a quaternion, inertial to body, scalar last); NaN before a sensor's first sample |
