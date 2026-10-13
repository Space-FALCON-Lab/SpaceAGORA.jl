# Legacy plotting saved-data contracts

The maintained raw-quaternion and body-rate helpers are explicit saved-result
consumers. They do not replace the full legacy plotting dispatcher. All 19
functions in `scripts/plotting/plot_data.jl` remain retained. A replacement for
one set of traces does not establish that an untracked legacy consumer can be
retired.

## Supported body-rate input

`_body_rate_history_samples(df; spacecraft_index)` and
`_save_body_rate_history_plot(csv_path; spacecraft_index, output_path)` consume
`time` in seconds and exactly `scN_body_rate_1`, `scN_body_rate_2`,
`scN_body_rate_3` in radians per second. The components are angular velocity
expressed in the spacecraft body frame. `N` is the one-based position in the
spacecraft list, not `SpacecraftModel.id`.

This schema has a tracked producer: `slew_save_fields` in
`scripts/dev/viewer_demos/cygnss_slew.jl` records `u.sc[1].ω` through
`SaveField(:body_rate_rad_s, ...; per_satellite=true, column_prefix="body_rate")`.
The same script reads `sc1_body_rate_1..3` and labels them `rad/s`.
`_assign_orientation_rhs!` in `src/simulation/engine/dynamics_rhs.jl` uses
`sc_view.ω` as `omega_body`; the external-propagator state contract also states
body angular velocity in rad/s. The result writer in
`src/io/outputs/io_outputs.jl` expands these saved three-vectors in component
order and names spacecraft by list position.

These are **custom, opt-in fields**. `default_save_fields` saves quaternions
when `orientation_sim=true`, but does not save body rates. The `:omega` fields
in several dynamics tests use a different explicit prefix. The plotting helper
does not silently accept `scN_omega_*`, unprefixed legacy `omega_*`, telemetry
names, or a quaternion-only table. Column names alone cannot certify an
external file's units or provenance; callers must supply the stated schema.

The maintained plot preserves the legacy angular-velocity plot's operation
order: `time / 60`, then `rad2deg.(rates) .* 60`. Its labels are `ω1`, `ω2`,
`ω3`, with red, green and blue traces, `Time [min]` and
`Angular Velocity (deg/min)` axes. It preserves row order, repeated/nonmonotone
times and rate signs. It copies raw rates into `rates_rad_s` without Float64
coercion and returns displayed arrays as `components`. No frame rotation,
smoothing, normalization, Euler conversion or quaternion differentiation occurs.

All selected columns must be present, nonempty, finite and real; Boolean values
are rejected. Unselected spacecraft columns are not validated. Singleton input
is valid. Nonfinite converted values are rejected before rendering. The facade
requires an explicit output extension, retains the established same-file guard
including hard links, symbolic links and `~` expansion, and never opens a plot
window. The raw-quaternion helper's behavior is unchanged.

## Remaining diagnostic boundaries

| Legacy diagnostic | Located producer or input | Remaining contract or capability gap |
| --- | --- | --- |
| `quaternion_plot` | Default `scN_q_1..4` for orientation runs; maintained raw-quaternion helper | Replacement traces are available; legacy retirement still requires a consumer decision. |
| `angular_velocity_plot` | Explicit `scN_body_rate_1..3` above; maintained body-rate helper | Default runs do not contain these fields. Other explicit schemas require a separate decision, not name guessing. |
| `attitude_plot` | Raw quaternion columns; legacy calls `qToEulerAngles` | No saved Euler history or agreed Euler convention/branch contract. A plot would add a mathematical transformation and needs separate approval and tests. |
| `reaction_wheel_h_plot` | Legacy nested `solution.physical_properties.rw_h`; current integrated state `h_wheels` and explicit test `:wheel_h` | Wheel identity/order and per-wheel momentum versus body-vector momentum need an output contract. The CYGNSS demo's `scN_wheel_momentum_nms_1..3` is aggregate body-frame momentum, not the legacy per-wheel history. |
| `reaction_wheel_torque_plot`, `total_reaction_wheel_torque_plot` | Legacy `rw_tau_i`, `rw_torque_ii_1..3`; CYGNSS explicit wheel speed, scalar reaction-torque magnitude and command-body torque fields | These are different quantities. Scalar magnitude or body torque cannot replace per-wheel or inertial torque components without an explicit producer/frame mapping. |
| `torque_plot`, `inertia_plot` | Legacy `tau_ii_1..3` and `J_ii_1..9`; runtime body torques/model inertia exist | No matching generic saved inertial histories located. RPO's `rpo_control_torque_body` is a held control command, not total inertial torque. Do not substitute model constants or zeros for a history. |
| `lambda_plot` | T_EDG appends matched `lambda_switch_list` and `time_switch_list` values | No tracked writes to costate histories `lamv_list`/`time_list` located and no standard saved-table fields. `constraint_tracking.jl` separately uses unsplatted vector pushes and `time_list2`, absent from the inspected runtime schema; retain that producer concern without changing it in a plotting packet. |
| `drag_passage_plot`, `performance_plots` | Modern heat maxima, orbit and mission plots cover only some panels | Per-link heat histories, limits, energy/fuel and orbit grouping are not equivalent to maxima. |
| `angle_of_attack_plot`, `wind_relative_attitude_plot` | Some mission-specific control/link fields | Generic per-link AoA/sideslip and normalized-orbit schemas remain undecided; wind-relative code also references free variable `m`. |
| `closed_form_solution_plot` | Nested guidance histories | Numerical/closed-form pairing, free `args`, units and passage segmentation need a separate guidance contract. |
| `traj_2D`, `traj_3D`, `ABM_periapsis`, `ground_track` | Maintained trajectory, geography and extrema plots cover parts | Orbital-plane transforms, passage surfaces, maneuver-magnitude axes and time-versus-orbit traces retain distinct capabilities. |
| `plots` | Legacy model/options/nested-solution dispatcher | No complete supported replacement dispatcher is established. |

The next safe work is a deliberate choice of one remaining producer/schema and
its consumers. Runtime storage alone, an optional custom field, or a successful
body-rate plot does not establish that all legacy diagnostics are supported.
