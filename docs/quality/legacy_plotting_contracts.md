# Legacy plotting saved-data contracts

`scripts/plotting/plot_data.jl` holds the legacy plotting dispatcher `plots` and
18 individual renderers; all 19 functions are still in the file. Two maintained
helpers now read current saved results directly: the raw-quaternion plot and the
body-rate plot, in `examples/support/aerobraking_result_tables.jl` and
`examples/support/aerobraking_mission_plot_utils.jl`. They do not replace the
legacy dispatcher, and covering some of its traces does not mean any legacy
renderer can be retired. This page states what each saved quantity means, which
code produces it, and what the plots may assume about it.

The tests are under `test/unit/`: `quaternion_history_plot_tests.jl`,
`body_rate_history_plot_tests.jl`, `legacy_euler_contract_tests.jl`,
`legacy_wheel_contract_tests.jl` and `legacy_guidance_diagnostic_contract_tests.jl`.

## Body-rate input

`_body_rate_history_samples(df; spacecraft_index)` and
`_save_body_rate_history_plot(csv_path; spacecraft_index, output_path)` read
`time` in seconds and exactly `scN_body_rate_1`, `scN_body_rate_2` and
`scN_body_rate_3` in rad/s: the angular velocity in the spacecraft body frame. `N`
is the one-based position in the spacecraft list, not `SpacecraftModel.id`.

These are custom save fields. `default_save_fields` saves the quaternion when
`orientation_sim=true` but saves no body rate. A run gets these columns from a
field such as
`SaveField(:body_rate_rad_s, …; per_satellite=true, column_prefix="body_rate")`
that saves `u.sc[i].ω`. The CYGNSS slew study in SpaceAGORAMissionReconstructions
uses that field; a copy of its save fields, `slew_save_fields` in
`test/fixtures/research_contracts/save_field_recipes.jl`, is what the body-rate
and wheel tests run. The engine's orientation RHS (`_assign_orientation_rhs!` in
`src/orchestration/simulation/engine/dynamics_rhs.jl`) treats `sc_view.ω` as the
body angular velocity. The result writer in
`src/interfaces/io/outputs/io_outputs.jl` expands each saved three-vector in
component order and names spacecraft by list position.

The helper does not accept `scN_omega_*`, the unprefixed legacy `omega_*`,
telemetry names or a quaternion-only table as substitutes. The `:omega` fields in
several dynamics tests use a different prefix. Column names cannot certify an
external file's units or origin; the caller must supply the stated schema.

The plot keeps the legacy angular-velocity plot's arithmetic: `time / 60`, then
`rad2deg.(rates) .* 60`. Its traces are labelled `ω1`, `ω2` and `ω3` in red,
green and blue, with axes `Time [min]` and `Angular Velocity (deg/min)`. Row
order, repeated or non-monotone times and signs are kept. The raw rates are
copied into `rates_rad_s` without conversion to Float64, and the displayed arrays
are returned as `components`. There is no frame rotation, smoothing,
normalization, Euler conversion or differentiation of quaternions.

The selected columns must be present and non-empty, and every value must be a
finite real number; Boolean values are rejected. Columns of other spacecraft are
not checked, and a single row is valid. Values that become non-finite after
conversion are rejected before plotting. `output_path` must have a file
extension and must not be the input CSV, including through hard links, symbolic
links or `~` expansion. The helper never opens a plot window. The raw-quaternion
helper, `_quaternion_history_samples` and `_save_quaternion_history_plot`, applies
the same input and output-path rules to `scN_q_1..4`.

## Euler angles

The default saved quaternion is scalar-last `[x, y, z, w]` and describes the
inertial-to-body rotation. `_save_quaternion` and the result writer keep its four
components as stored; no Euler angles are saved.

The legacy `attitude_plot` reads the unprefixed `q_1..4`, not the current
`scN_q_1..4`, and calls `qToEulerAngles` (in
`src/foundation/core/numerics/quaternion_utils.jl`) on every row. That function
normalizes each quaternion with `LinearAlgebra.normalize`, without the canonical
near-zero identity fallback or any input validation. A tiny non-zero quaternion
is normalized; zero or non-finite input can give NaN angles. It does not clamp
the pitch intermediate to its domain, so rounding near gimbal lock can raise a
`DomainError` even for finite input.

The triple it returns is a legacy conversion, not a certified roll-pitch-yaw
convention. The identity quaternion gives `[pi/2, 0, pi/2]`: both `atan(y, x)`
calls take their arguments in the opposite order to the usual roll and yaw
formulas, and the pitch expression has the opposite sign. Rescaling a quaternion
or reversing its sign leaves the result unchanged, but there is no branch
unwrapping, singularity handling or continuity guarantee. The renderer keeps
sample order and plots the values in degrees against minutes, labelled `Roll`,
`Pitch` and `Yaw` in red, green and blue. The labels do not make the angles a
conventional Euler sequence.

Keep this behavior as a tested compatibility boundary, and use the raw-quaternion
plot to inspect saved orientation. A conventional Euler plot would need a
separately named function with an explicit axis sequence, frame direction,
active or passive convention, branch interval and singularity policy. Changing
`qToEulerAngles` in place would change numbers outside plotting. No maintained
quaternion-to-Euler alternative exists yet.

## Wheel momentum and torque

### Integrated per-wheel momentum

The root reaction-wheel assembly's integrated state is `u.sc[i].h_wheels`, in
N m s (kg m²/s). Component `j` follows column `j` of that assembly's `J_rw`
matrix. It does not follow `SpacecraftModel.id`, a physical wheel identifier, or
a wheel count summed over child links. Store physical wheel names and the axis
matrix separately when they matter. This describes the root assembly only; it
does not establish support for integrating wheels on every articulated link.

A caller can save a per-wheel history with the existing `SaveField` interface,
after checking that every spacecraft in the run has orientation enabled and a
root assembly with the expected wheel count, fixed for the whole run:

```julia
wheel_h_field = SpaceAGORA.SimulationModel.SaveField(
    :wheel_h_nms,
    (u, t, integrator) -> [
        copy(u.sc[i].h_wheels)
        for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
    ];
    per_satellite=true,
    column_prefix="wheel_h_nms",
)
```

Add it to the run's save fields. The writer then produces
`scN_wheel_h_nms_1..K` at each saved `time` in seconds, where `N` is the list
position and `K` that spacecraft's root wheel count. The `copy` keeps earlier
snapshots from changing with the state. `legacy_wheel_contract_tests.jl` runs this
field through the snapshot collector and the table writer, with unequal,
non-consecutive spacecraft IDs. It is not a default field, and no plot uses it
yet. A future consumer must be given the wheel order and count explicitly; a
vector with three components does not imply three wheels.

The engine computes the body-frame aggregate momentum as `J_rw * h_wheels` and
the integrated wheel momentum rate as `J_rw_pinv * (-rw_torque_body)`, so the
wheel momentum rate has the opposite sign to the requested wheel body torque. For
a redundant assembly the pseudoinverse is the engine's allocation; for a
rank-deficient one, multiplying back by `J_rw` can project the requested vector. A
plot must not present a command as realized torque or as a statement about limits
or feasibility.

### Other wheel quantities

The test field `:wheel_h` in `test/unit/dynamics/reaction_wheel_momentum_tests.jl`
evaluates `wheel_momentum_body(model, t)`: the aggregate body-frame momentum
vector, not the integrated per-wheel state. The slew recipe's
`scN_wheel_momentum_nms_1..3` means the same. The prescribed-wheel model sums
signed speed, inertia and axis contributions in its input column order and returns
the sum in the body frame. Its reaction torque, `wheel_reaction_torque`, is
`-Ḣ_w - ω × H_w`. The magnetic momentum manager's `h_wheels` is yet another
body-frame accumulator. Similar names do not make these histories interchangeable.

| Quantity | Where it comes from | What a plot may do with it |
| --- | --- | --- |
| Integrated per-wheel momentum | `u.sc[i].h_wheels`, through the opt-in field above | Plot per wheel in root-assembly order; keep signs and sample order. |
| Prescribed aggregate momentum | `wheel_momentum_body(model, t)` | A body-frame vector; never label its components as three wheels. |
| Per-wheel torque or momentum rate | `du.sc[i].h_wheels`, computed during RHS evaluation | No save-time producer exists. Do not read the assembly's stale `h_dot_wheels` as a history. |
| Requested wheel body torque | `rw_torque_body`, passed from the control readout to the orientation RHS | May be an RPO held command (`RPOHeldActuation.rw_torque_body`). It is the request before the pseudoinverse, not realized, allocated or total spacecraft torque. |
| Prescribed reaction torque | The prescribed model's force and torque | Includes the transport term; not a wheel motor torque vector. |
| Wheel torque magnitude in the slew recipe | `scN_wheel_torque_nm`, the norm of the prescribed reaction torque | A magnitude only; it cannot supply signed components or a frame. |
| RPO held control torque | `RPOHeldActuation.torque_body` on `RPOMPCControlModel.held` (thruster torque plus requested wheel torque) | A held command, not total inertial torque. |
| Total torque, inertial inertia | Legacy `tau_ii_1..3`, `rw_torque_ii_1..3`, `J_ii_1..9` | No current producer saves them with matching times, and there is no agreed tensor flattening. Do not substitute static model inertia, body components or zeros. |

The legacy `reaction_wheel_h_plot` expects a nested
`solution.physical_properties.rw_h` with one history per wheel.
`reaction_wheel_torque_plot` expects unprefixed `rw_tau_i` columns and the legacy
model's wheel count. Both use seconds and keep trace order and signs. Both pick
colours with `i % 9` from a nine-entry palette, which gives index zero for the
ninth wheel, so nine-wheel input throws a `BoundsError`; the test records this.
`total_reaction_wheel_torque_plot`, `torque_plot` and `inertia_plot` also read
their unprefixed legacy columns directly, and their names and labels do not prove
a producer's frame or matrix order.

The opt-in per-wheel momentum field and the prescribed aggregate momentum can each
be plotted on their own terms. Per-wheel torque and inertial torque or inertia
plots first need save-time producers and explicit frame and tensor metadata. The
allocator, the dynamics and the legacy renderers are unchanged.

## Costates and switching times

The legacy `lambda_plot(name, cnf, args)` draws two trace pairs straight from the
runtime configuration: `time_switch_list` against `lambda_switch_list`, and
`time_list` against `lamv_list`. It labels them in seconds and as lambda, and does
no conversion, shape checking, sorting, normalization or splitting by passage.
Lambda units and scaling depend on the producer; the labels do not establish them.
`src/orchestration/simulation/state/runtime_types.jl` initializes all four as
empty `Float64` vectors.

No code in this repository writes or resets these arrays, so empty arrays mean
that no costates were recorded, not a zero costate. The code that used to fill
them was removed with the older E-EDG/T-EDG strategy routes. Those arrays carried
no spacecraft or passage identity, so a history saved by older code is an
unsegmented log in elapsed simulation seconds.

Typed EDG keeps its predicted lambda and alpha data local to each call, and can
save a predicted targeting switch time. The Odyssey example
`examples/AGORA_Odyssey_Energy_Depletion_Control_Test.jl` defines
`targeting_switch_save_field`, which saves `state.targeting_switch_s` as
`scN_edg_targeting_switch_s`; `time` is the saved sample time and `N` the
spacecraft list position. Targeting brackets candidate switch times between the
current elapsed simulation time and that time plus the predicted duration, with a
one-second upper margin. Separately, the typed heat-load prediction tracks use
relative times, and SpaceAGORAEDG's `_edg_low_alpha_switch_window` adds the
current simulation time to the selected track times; that heat-load window does
not produce the saved targeting field.

The field is a predicted switch time in elapsed simulation seconds. It is not
lambda, a time-to-go, an actual switch event or a costate history. The non-finite
values written on reset or when there is no prediction must stay distinguishable
from finite predictions; never replace them with zero. A diagnostic should label
those samples unavailable and treat the values as a held prediction only if that
matches how they are updated. On a pass reset, `_edg_reset_pass!` sets that
spacecraft's `targeting_switch_s` to `Inf`, `heat_load_switches_s` to
`(Inf, Inf)` and `last_switch_solve_t` to `-Inf`, and keeps the cumulative solve
counts and the last telemetry. `_edg_initialize_state!` resets the whole state,
counts and telemetry included. `legacy_guidance_diagnostic_contract_tests.jl`
exercises both through the reset functions and the saver. Combining old solver
samples with typed targeting predictions would need an explicit passage and
spacecraft mapping.

Keep the legacy lambda renderer, and treat its input as unpopulated. Use the
typed targeting-switch field only as the prediction diagnostic described above.
Rebuilding a costate history needs a supported producer, a costate definition and
scaling, an independent variable, passage identity and a reset policy; those are
scientific decisions for the EDG owners, and no plot should synthesize the
history.

## Other legacy renderers

| Legacy renderer | Producer or input found | What is still missing |
| --- | --- | --- |
| `quaternion_plot` | The default `scN_q_1..4` of orientation runs; the maintained raw-quaternion helper | The replacement traces exist; retiring the renderer still needs a decision from its users. |
| `angular_velocity_plot` | The explicit `scN_body_rate_1..3` above; the maintained body-rate helper | Default runs do not contain these fields, and the helper does not guess a schema. |
| `drag_passage_plot`, `performance_plots` | Current heat-maximum, orbit and mission plots cover some panels | Per-link heat histories, limits, energy and fuel, and orbit grouping are not the same as maxima. |
| `angle_of_attack_plot`, `wind_relative_attitude_plot` | Some mission-specific control and link fields | Generic per-link angle of attack and sideslip, and normalized-orbit schemas, are undecided; `wind_relative_attitude_plot` also uses an undefined variable `m`. |
| `closed_form_solution_plot` | Nested guidance histories | Pairing numerical with closed-form results, its use of an undefined `args`, units and passage splitting need a separate guidance contract. |
| `traj_2D`, `traj_3D`, `ABM_periapsis`, `ground_track` | Current trajectory, geography and extrema plots cover parts | Orbital-plane transforms, passage surfaces, maneuver-magnitude axes and time-versus-orbit traces have no replacement. |
| `plots` | The legacy dispatcher over model, options and nested solution | No complete supported replacement exists. |

These contracts separate what the producers do today from any approved
mathematical correction or maintained replacement. They settle what the Euler,
wheel and torque, and costate and switching data mean, without declaring the
legacy diagnostics supported and without retiring any of the 19 functions.
