# Output Columns

This page lists the columns of the saved results table with their units and
frames. [Simulation Outputs](outputs.md) describes the files that hold them,
and [The Integrated State](integrated_state.md) says which values are
integrated and which are derived.

## CSV column reference

Every per-spacecraft column starts with `sc{N}_`, where `N` is the spacecraft's
one-based position in the list passed to `DynamicsModel`. `SpacecraftModel.id`
does not set the prefix, so a single-spacecraft run produces `sc1_*` columns
whatever its spacecraft's `id`. The tables below show spacecraft 1.

With `results_layout = :long`, the spacecraft table has the same columns
without the prefix (`pos_1` for `sc{N}_pos_1`), one row per spacecraft, and
the run-wide columns (those without a prefix) are in the run table. See
[Long layout for many spacecraft](outputs.md#Long-layout-for-many-spacecraft).

### Time

| Column | Unit | Description |
|---|---|---|
| `time` | s | Elapsed simulation time from the epoch |

### Position and velocity (inertial frame)

| Column | Unit | Description |
|---|---|---|
| `sc1_pos_1` | m | Inertial position, X component |
| `sc1_pos_2` | m | Inertial position, Y component |
| `sc1_pos_3` | m | Inertial position, Z component |
| `sc1_vel_1` | m/s | Inertial velocity, X component |
| `sc1_vel_2` | m/s | Inertial velocity, Y component |
| `sc1_vel_3` | m/s | Inertial velocity, Z component |

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

Here `a` and `e` are the semimajor axis and eccentricity of the osculating
Keplerian orbit derived from the current inertial position and velocity. The
reference surface is a sphere with the planet's equatorial radius,
`planet.Rp_e`. The column describes the instantaneous osculating orbit; it does
not predict the lowest geodetic altitude the propagated trajectory reaches.
Because `sc1_altitude` and `sc1_periapsis_altitude` use different reference
surfaces, they need not agree at a periapsis.

The aerobraking examples write separate `periapsis_events.csv` tables, which
hold geodetic estimates instead. Depending on the example and on the solver
output available, these come from event-located radial minima or from the
lowest sampled geodetic altitudes. The examples' apoapsis plots may likewise
use interpolated crossings between samples and an initial sample. None of these
event or sampled values is `sc1_periapsis_altitude`.

### Mass

| Column | Unit | Description |
|---|---|---|
| `sc1_mass` | kg | Total spacecraft mass (structure plus remaining propellant) |

### Atmosphere and aerodynamics

| Column | Unit | Description |
|---|---|---|
| `sc1_wind_1` | m/s | Local wind, eastward component |
| `sc1_wind_2` | m/s | Local wind, northward component |
| `sc1_wind_3` | m/s | Local wind, upward component |
| `sc1_drag_1` | N | Aerodynamic drag force, X component (inertial frame) |
| `sc1_drag_2` | N | Aerodynamic drag force, Y component (inertial frame) |
| `sc1_drag_3` | N | Aerodynamic drag force, Z component (inertial frame) |
| `sc1_lift_1` | N | Aerodynamic lift force, X component (inertial frame) |
| `sc1_lift_2` | N | Aerodynamic lift force, Y component (inertial frame) |
| `sc1_lift_3` | N | Aerodynamic lift force, Z component (inertial frame) |
| `sc1_cross_1` | N | Aerodynamic cross force, X component (inertial frame) |
| `sc1_cross_2` | N | Aerodynamic cross force, Y component (inertial frame) |
| `sc1_cross_3` | N | Aerodynamic cross force, Z component (inertial frame) |

### Thermal

| Column | Unit | Description |
|---|---|---|
| `sc1_heat_rate` | W/cm² | Largest of the per-link stagnation heat rates (the built-in Maxwellian model returns W/cm²) |
| `sc1_heat_load` | J/cm² | Largest of the per-link accumulated heat loads (each link integrates its own rate, without an area-unit conversion); not the sum over links |

### Attitude (orientation_sim only)

These columns are written only when `mission_configuration.orientation_sim = true`:

| Column | Unit | Description |
|---|---|---|
| `sc1_q_1` | — | Attitude quaternion component `x` (scalar-last convention `[x, y, z, w]`, inertial to body) |
| `sc1_q_2` | — | Attitude quaternion component `y` |
| `sc1_q_3` | — | Attitude quaternion component `z` |
| `sc1_q_4` | — | Attitude quaternion scalar component `w`; `[0, 0, 0, 1]` is the identity attitude |

#### Plot saved raw quaternion history

The example helpers can plot these columns from a saved CSV. From the
repository root, load them and name the spacecraft and the image file
explicitly. The paths below assume an orientation run that wrote to `output/`:

```julia
include("examples/common.jl")
include("examples/support/aerobraking_mission_plot_utils.jl")
quaternion_history = _save_quaternion_history_plot(
    "output/simulation_results.csv";
    spacecraft_index=1,
    output_path="output/raw_quaternion.png",
)
```

`spacecraft_index` is the positive, one-based position in the spacecraft list,
not `SpacecraftModel.id`. The helper needs `time` (seconds) and all four
`scN_q_1` to `scN_q_4` columns for that spacecraft, and the default save fields
write those only with `orientation_sim=true`. Missing columns, an empty table,
and missing, non-numeric or non-finite values in the selected columns raise an
error; columns of other spacecraft are not read.

The plot keeps the saved signs, component values and row order, including
non-unit and zero quaternions. Time is shown in minutes (`time / 60`), and a
single row, repeated times or times out of order are accepted. The helper does
no normalization, sign adjustment, interpolation or attitude conversion. The
traces are labeled `q1` to `q4` and colored red, green, blue and orange. The
returned value holds `time_min`, `components`, `labels`, `colors`,
`spacecraft_index`, `plot` and `output_path`, so you can inspect the plotted
data without reading values off the image. For a `DataFrame` already in
memory, `_quaternion_history_samples(df; spacecraft_index=1)` prepares the same
data without plotting or writing a file.

The call writes only the requested image and opens no plot window. Give the
file name an extension Plots supports (`.png` or `.pdf`, for example) and an
existing directory; the input CSV cannot be the destination. The helper runs
no simulation and leaves the example plots unchanged.

### Body angular-rate history (explicit save fields)

`default_save_fields` does not save body rates, even for orientation runs. For
an orientation configuration `args` you have already built, this save field
records the body angular velocity state in radians per second:

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

This writes `scN_body_rate_1..3` in body-axis component order, where `N` is the
spacecraft's one-based position in the list, not its `id`. The CYGNSS slew
study in the separate SpaceAGORAMissionReconstructions repository records the
same columns for its single spacecraft. Custom `:omega` save fields written
elsewhere use a different layout, and the helpers do not treat them as the
same columns.

To plot a saved CSV, use the example helpers:

```julia
include("examples/common.jl")
include("examples/support/aerobraking_mission_plot_utils.jl")
body_rate_history = _save_body_rate_history_plot(
    "output/simulation_results.csv";
    spacecraft_index=1,
    output_path="output/body_rate.png",
)
```

The helper needs `time` in seconds and all three body-rate columns of the
selected spacecraft in rad/s. It plots time in minutes and angular velocity in
deg/min, converted as `rad2deg.(rates) .* 60`, with the labels `ω1`, `ω2` and
`ω3` in red, green and blue. It keeps the signs and row order, including
repeated times or times out of order, and it neither derives rates from
quaternions nor rotates frames.

The returned `rates_rad_s` holds copies of the saved values and `components`
the plotted deg/min values; the other fields are `time_min`, `labels`,
`colors`, `spacecraft_index`, `plot` and `output_path`.
`_body_rate_history_samples(df; spacecraft_index=1)` prepares the arrays
without plotting. An empty table, missing columns, invalid selected values, or
values that are not finite after conversion raise an error, and columns of
other spacecraft are not read. Give the file name an extension and an existing
directory. The input CSV cannot be the destination, also through a file-system
alias. The helper opens no plot window and runs no simulation.

### Visualization fields

These columns are added when `simulation_settings.save_visualization_scene = true`
and the run can supply them; the
[sidecar section](outputs.md#Visualization-scene-sidecar-(opt-in)) gives the
conditions for each:

| Column | Unit | Description |
|---|---|---|
| `sc1_link_pose_1` .. `sc1_link_pose_7n` | m, — | `[rx, ry, rz, qx, qy, qz, qw]` per non-root link, relative to the root bus, in the sidecar's link order |
| `sc1_thruster_level_1` .. `sc1_thruster_level_k` | — | Firing level of each thruster, 0 to 1, in the sidecar's thruster order |
| `sc1_thruster_direction_1` .. `sc1_thruster_direction_3k` | — | Exhaust direction of each thruster (a unit vector in its link's frame), in the sidecar's thruster order; written when a spacecraft carries a thruster gimbal, which turns it |
| `sc1_arm_pose_1` .. `sc1_arm_pose_7m` | m, — | `[rx, ry, rz, qx, qy, qz, qw]` per robot-arm link: center of mass relative to the spacecraft (inertial axes) and inertial attitude |
| `sc1_density` | kg/m³ | Atmospheric density at the spacecraft; zero outside the atmosphere |
| `sun_dir_1` .. `sun_dir_3` | — | Unit vector from the planet's center to the Sun, inertial frame, one per row (not per spacecraft) |

### Articulated spacecraft and attachments

These columns are written with the default save-field list, with or without the
visualization sidecar. Articulated spacecraft (any non-fixed `Joint`) and
spacecraft with compliant attachments get them; other spacecraft have no
columns of these fields. Quaternions are scalar-last and body-to-inertial.

| Column | Unit | Description |
|---|---|---|
| `sc1_joint_q_1` .. | rad, m | Coordinates of the non-fixed joints in `spacecraft.joints` order: an angle for a `:hinge`, a displacement for a `:slide`, four quaternion components for a `:ball` |
| `sc1_joint_qd_1` .. | rad/s, m/s | Rates of the same joints; a `:ball` joint takes three columns |
| `sc1_articulated_link_pose_1` .. `sc1_articulated_link_pose_7n` | m, — | `[rx, ry, rz, qx, qy, qz, qw]` per link in `spacecraft.links` order, from the joint kinematics: the link's center of mass in inertial meters and its attitude, including links merged through `:fixed` joints |
| `sc1_attachment_pose_1` .. `sc1_attachment_pose_7n` | m, — | `[rx, ry, rz, qx, qy, qz, qw]` per attachment body in attachment order: the body's center of mass in inertial meters (`pos + att_r`) and its attitude |
| `sc1_system_com_1` .. `sc1_system_com_3` | m | Inertial center of mass of the whole spacecraft, including attachment bodies |

`system_com` comes from the articulated fields for an articulated spacecraft
and from the attachment fields for a rigid spacecraft with attachments. The
`link_pose` visualization field still holds the configured rigid geometry, not
these poses.

### Plume interaction

The plume-surface effector ([Plume interaction](plume_interaction.md)) writes
these columns for every spacecraft in the run when the default save-field list
is used, whether or not the visualization sidecar is on. An explicit list must
include them to write them:

| Column | Unit | Description |
|---|---|---|
| `sc1_plume_height_m` | m | Slant height of the vehicle reference point above the terrain along the engine axis (radial clearance divided by the downward exhaust cosine); infinite when the exhaust does not point toward the surface |
| `sc1_plume_pressure_pa` | Pa | Peak surface pressure under the plume |
| `sc1_plume_shear_pa` | Pa | Peak wall shear stress under the plume |
| `sc1_plume_erosion_kg_s` | kg/s | Regolith mass erosion rate |
| `sc1_plume_eroded_kg` | kg | Eroded mass integrated over accepted solver steps (trapezoidal rule), reset at the start of each run; a resumed run starts a new integral at its resume time |
| `sc1_plume_ejecta_mps` | m/s | Bounded speed at which eroded grains leave the surface |
| `sc1_plume_ground_effect_n` | N | Modeled thrust augmentation in ground effect |

The plume quantities are diagnostics of an illustrative closure, as the plume
page explains; they are not a calibrated flight reconstruction.

### Mounted actuators

These columns are written with the default save-field list when a control
effector of the run is in `:mounted` actuation (see
[Extensibility](../extensibility.md#Add-an-actuator)) or an interlink drives
terminal hardware, for every spacecraft that carries mounted actuators. An
explicit list must include the `actuator_level` and `actuator_power` fields to
write them:

| Column | Unit | Description |
|---|---|---|
| `sc1_actuator_level_1` .. `sc1_actuator_level_k` | the kind's | The level each mounted actuator applies, in mounting order: its command within its limits. A thruster with firing dynamics reports the average thrust its last control tick scheduled, a wheel at its momentum limit reports zero when the command would store more, a hinge or thruster gimbal reports its current angle, a control moment gyro its gimbal rate (zero at a gimbal stop the command would pass), and a terminal an interlink drives the level its link sets |
| `sc1_actuator_power` | W | Electrical power the mounted actuators draw at those levels, the sum of what each kind reports through `actuator_consumption` |

The levels come from the first control effector in `:mounted` actuation that
holds levels for the spacecraft; without one, they are zero. A row saved
before that effector's first control tick reports the levels it held then. A
terminal that an interlink drives is the exception: whatever the effector
holds, it shows the level its active link sets, and zero while no link
drives it or while the battery sheds it. `sc1_actuator_power` likewise adds
the power those terminals draw to the effector's actuators.

Among the built-in kinds, laser and radio terminals draw `level / efficiency`.
Reaction wheels draw their base power, plus their motor's power computed from
the integrated wheel momentum when they are given a spin inertia. Torque rods
draw their power at full dipole times the square of the dipole fraction.
Hinges, gimbals, control moment gyros and thrusters draw none.

Runs with control moment gyros and `orientation_sim = true` also write the
gimbal angles from the integrated state:

| Column | Unit | Description |
|---|---|---|
| `sc1_cmg_angle_1` .. `sc1_cmg_angle_n` | rad | Gimbal angle of each control moment gyro, in mounting order, then the inner gimbal angle of each double-gimbal CMG |

### Power budget

These columns are written with the default save-field list when a spacecraft
in the run carries power-system parts (solar arrays, a battery, constant loads
or a `PowerReceiver`; see [Powering a Spacecraft](power_budget.md)). They are
NaN for a spacecraft without such parts, and the battery columns are NaN for a
spacecraft without a battery:

| Column | Unit | Description |
|---|---|---|
| `sc1_battery_wh` | Wh | Energy stored in the battery (an integrated state) |
| `sc1_battery_soc` | — | State of charge: the stored energy over the capacity |
| `sc1_power_generated_w` | W | Power the solar arrays generate (lit cells facing the Sun, through the eclipse model), plus any power received from power beams |
| `sc1_power_load_w` | W | Load on the bus, with any shedding applied: constant loads, the mounted actuators of every control effector in `:mounted` actuation (as the dynamics apply each one's actuators), and terminals an interlink drives |
| `sc1_power_shed_w` | W | Actuator power shed because the battery is empty |
| `sc1_power_net_w` | W | Generation minus load; negative while the battery is discharging, or while it is empty and the load exceeds the generation |
| `sc1_power_received_w` | W | Electrical power the spacecraft's incoming power beams deliver (part of the generation); only in a run with a `PowerBeamLink` |
| `sc1_power_beamed_w` | W | Optical power the spacecraft's outgoing power beams carry; only in a run with a `PowerBeamLink` |

### Data links

These columns are written with the default save-field list in a run with a
`DataLink` (see [Interlinks](interlinks.md#Data-links)):

| Column | Unit | Description |
|---|---|---|
| `sc1_data_rate_tx_bps` | bit/s | Data rate the spacecraft's active outgoing links carry |
| `sc1_data_rate_rx_bps` | bit/s | Data rate its active incoming links carry |
| `sc1_data_sent_bits` | bit | Bits sent since the run started (an integrated state) |
| `sc1_data_received_bits` | bit | Bits received since the run started (an integrated state) |

With data stores, the rates are those the links carry as the stores limit them.

With ground stations (see [Interlinks](interlinks.md#Ground-stations)), each
station `s` adds two columns after the spacecraft columns:

| Column | Unit | Description |
|---|---|---|
| `gs1_received_bits` | bit | Bits station 1 has received since the run started (an integrated state) |
| `gs1_sent_bits` | bit | Bits station 1 has sent since the run started (an integrated state) |

### Data storage

These columns are written with the default save-field list when a spacecraft
in the run carries a `DataStore` (see
[Interlinks](interlinks.md#Data-storage)). They are NaN for a spacecraft
without a store:

| Column | Unit | Description |
|---|---|---|
| `sc1_data_stored_bits` | bit | Data the store holds (an integrated state) |
| `sc1_data_store_fill` | — | That data over the store's capacity |
| `sc1_data_produced_bits` | bit | Data the spacecraft's sources have produced since the run started (an integrated state) |
| `sc1_data_lost_bits` | bit | Production a full store could not keep, since the run started (an integrated state) |

### Interlinks

These columns are written with the default save-field list in a run with an
interlink model (see [Interlinks](interlinks.md)), whatever its link kinds.
They stay zero for links that apply no force, such as power beams and data
links. An explicit list must include the `interlink_dv` and
`interlink_delta_sma` fields to write them:

| Column | Unit | Description |
|---|---|---|
| `sc1_interlink_dv_1` .. `sc1_interlink_dv_3` | m/s | Δv the selected links' forces have given since the run started, inertial frame (an integrated state) |
| `sc1_interlink_delta_sma` | m | The selected links' contribution to the osculating semimajor axis since the run started: the integral of `2a²(v·F)/(μm)` (an integrated state), not a comparison against an unforced trajectory |
| `sc1_laser_dv_1` .. `sc1_laser_dv_3`, `sc1_laser_delta_sma` | as above | Deprecated aliases with the same values, under the names used before 0.6.0; they will be removed in 1.0. A column cannot warn, so leave the `laser_dv` and `laser_delta_sma` fields out of an explicit list to stop writing them |

### Link schedule

A run with an interlink model writes its link schedule as a separate table,
`simulation_results_link_schedule.csv` (with `save_csv`) and
`simulation_results_link_schedule.feather`, beside the results.
`link_schedule` returns the same table (see
[The link schedule](interlinks.md#The-link-schedule)). It has one row per
window, sorted by `start_s`, then active rows before contact rows, then by
link key:

| Column | Unit | Description |
|---|---|---|
| `window` | — | `active`: a span over which the link was selected; `contact`: a span over which it was available, present only when the history records the available links |
| `link` | — | The link's number, its position in registration order |
| `kind` | — | The link kind's type name with its type parameters, for example `LaserMomentumLink` or `DataLink{ShannonBudget}` |
| `directed` | — | Whether the link is directed |
| `a_spacecraft`, `a_terminal` | — | The first endpoint of the link's key: the transmitter of a directed link, the lower `(spacecraft, terminal)` of an undirected one |
| `b_spacecraft`, `b_terminal` | — | The second endpoint; a ground station appears as minus its station number |
| `a_name`, `b_name` | — | Only in runs with stations: the station's registered name, or `spacecraft N`, for each endpoint |
| `start_s` | s | The span's first accepted step, on the clock of the results' `time` column |
| `end_s` | s | The next accepted step without the link, or the run's last one |
| `duration_s` | s | `end_s - start_s` |
| `steps` | — | The accepted steps in the span |
| `score_start`, `score_end` | the policy's | The score that selected the link at the span's first and last step; NaN for a contact row |
| `score_min`, `score_max` | the policy's | The score's extremes over the span; NaN for a contact row |
| `term1_start`, `term2_start`, ... | each term's | With `InterlinkModel(; record_terms=true)` and a `WeightedPolicy`: each term's own score at the span's first step |

The policy's unit is m/s for `SemimajorAxisRatePolicy`, 1/s for
`EccentricityRatePolicy`, W for `PowerDeliveryPolicy` and
`EnergyBalancingPolicy`, bit/s for `DataThroughputPolicy` and
`DataBacklogPolicy`, and 1 for `WeightedPolicy`. The results manifest names
the policy and its unit under `[interlinks]` (see
[Output files](outputs.md#Output-files) for its keys).

The schedule has the scheduler's resolution, the run's accepted steps, so a
window's ends are accurate to one step. It covers the span of the interlink
history, from `history_start_s` to `history_end_s` in `[interlinks]` and in the
table's metadata. A resumed run restores the usable history kept beside its
checkpoint. Without usable history, a resumed run records only the interval
after the resume, and a resume that solves nothing returns an empty schedule.

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
| `sc1_sensor_1` .. `sc1_sensor_n` | the sensor's | The latest measurement of every sensor: sensors in mounting order, each sensor's components in order |

Each sensor kind contributes these components:

| Sensor | Components | Unit | Frame |
|---|---|---|---|
| `RateGyro` | three body rates | rad/s | its link's frame |
| `Magnetometer` | three magnetic field components | T | its link's frame |
| `SunSensor` | unit vector toward the Sun | — | its link's frame |
| `Accelerometer` | three specific-force components | m/s² | its link's frame |
| `StarTracker` | attitude quaternion, inertial to body, scalar last | — | — |
| `GNSSReceiver` | position, then velocity | m, m/s | inertial, or planet-fixed with `frame=:planet_fixed` |

A sensor's columns are NaN before its first sample. A sun sensor's are also NaN
whenever the Sun is outside its field of view or eclipsed, and a GNSS
receiver's whenever it is above its `max_altitude_m`.
