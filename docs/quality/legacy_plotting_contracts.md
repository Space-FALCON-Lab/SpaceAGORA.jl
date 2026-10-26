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
the CYGNSS slew study (SpaceAGORAMissionReconstructions) records `u.sc[1].ω` through
`SaveField(:body_rate_rad_s, ...; per_satellite=true, column_prefix="body_rate")`.
The same script reads `sc1_body_rate_1..3` and labels them `rad/s`.
`_assign_orientation_rhs!` in `src/orchestration/simulation/engine/dynamics_rhs.jl` uses
`sc_view.ω` as `omega_body`; the external-propagator state contract also states
body angular velocity in rad/s. The result writer in
`src/interfaces/io/outputs/io_outputs.jl` expands these saved three-vectors in component
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

## Euler compatibility contract

The default saved quaternion is scalar-last `[x, y, z, w]`, with the
inertial-to-body rotation convention. `_save_quaternion` and the result writer
preserve its four stored components. They do not save Euler angles.

The retained `attitude_plot` consumes **unprefixed** `q_1..4`, not the current
`scN_q_1..4` schema, and calls the actual `qToEulerAngles` helper for every row.
That helper normalizes each quaternion with `LinearAlgebra.normalize`. It has
neither the canonical near-zero identity fallback nor input validation. A tiny
nonzero quaternion is normalized; zero and nonfinite inputs can yield NaN
angles. It does not clamp the pitch intermediate to its mathematical domain,
so rounding near gimbal lock can cause a `DomainError` even for finite input.

Its returned triple is a **legacy conversion**, not a certified roll/pitch/yaw
convention. In particular, the identity quaternion returns `[pi/2, 0, pi/2]`:
the two `atan(y, x)` calls exchange the arguments expected for conventional
roll and yaw, and the pitch expression has the opposite sign. Ordinary finite
rescaling and quaternion sign reversal preserve the represented result, but
there is no branch unwrapping, singularity stabilization or continuity promise.
The legacy renderer preserves sample order and plots these values in degrees
against minutes, using the labels `Roll`, `Pitch`, `Yaw` and red, green, blue.
Those labels do not establish a conventional Euler interpretation.

**Disposition:** keep this behavior as a tested compatibility boundary. Use the
maintained raw-quaternion plot when inspecting current saved orientation.
A conventional Euler plot requires a separately named, explicitly chosen axis
sequence, frame direction, active/passive convention, branch interval and
singularity policy. Correcting `qToEulerAngles` in place would change numerical
behavior outside plotting and is not part of this cleanup contract. No matching
maintained quaternion-to-Euler alternative has been established.

## Wheel momentum and torque contracts

### Integrated per-wheel momentum

The root reaction-wheel assembly's integrated state is `u.sc[i].h_wheels` in
N m s (equivalently kg m^2/s). Component `j` follows column `j` of that root
assembly's `J_rw` matrix. It does not follow `SpacecraftModel.id`, a physical
wheel identifier, or the total wheel count summed across child links. Persist
physical wheel names and the axis matrix separately when those identities
matter. This is the root assembly contract; it does not establish support for
integrating wheels on every articulated link.

A caller can opt in to a saved per-wheel history using the existing `SaveField`
interface after checking that each spacecraft in the run has orientation enabled
and a root assembly with the expected wheel count, fixed across saved rows:

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

Include this field in the run's configured save fields. The writer produces
`scN_wheel_h_nms_1..K` at the saved `time` in seconds; `N` is list position and
`K` is that spacecraft's root wheel count. The copy protects earlier snapshots
from subsequent state mutation. This opt-in schema is exercised through the
actual snapshot collector and table writer, including unequal/nonconsecutive
spacecraft IDs. It is not a default field or a new plotting facade. A future
consumer must receive the expected wheel order/count explicitly and must not
infer per-wheel semantics merely because a vector has three components.

The engine computes body-frame aggregate momentum as `J_rw * h_wheels` and
integrated wheel momentum rate as `J_rw_pinv * (-rw_torque_body)`. Thus wheel
momentum rate has the opposite reaction sign to the requested wheel body torque.
For a redundant assembly, the pseudoinverse gives the engine's allocation;
for a rank-deficient assembly, multiplication back by `J_rw` can project the
requested vector. A plot must not silently turn a command into a claim about
realized torque, limits or actuator feasibility.

### Other wheel quantities are separate

The test field named `:wheel_h` in `reaction_wheel_momentum_tests.jl` evaluates
`wheel_momentum_body(model, t)`. It is an **aggregate body-frame three-vector**,
not the integrated per-wheel state. The CYGNSS demo's
`scN_wheel_momentum_nms_1..3` has the same aggregate interpretation. The
prescribed-wheel model computes signed speed/inertia/axis contributions in its
input column order and returns their sum in the body frame; its reaction torque
includes both the negative momentum derivative and `-cross(omega_body, H)`.
The magnetic momentum manager's `h_wheels` is another body-frame three-vector
accumulator. Similar field names do not make these histories interchangeable.

| Quantity | Actual boundary | Decision for plotting |
| --- | --- | --- |
| Integrated per-wheel momentum | `u.sc[i].h_wheels`, opt-in field above | Usable producer with explicit root wheel ordering; preserve signs and saved sample order. |
| Prescribed aggregate momentum | `wheel_momentum_body(model, t)` | Body-frame vector; never label its three components as three wheels. |
| Per-wheel torque or momentum rate | Engine computes `du.sc[i].h_wheels` during RHS evaluation | No corresponding generic save-time producer is established. Do not read stale assembly `h_dot_wheels` as an output history. |
| Requested reaction-wheel body torque | `rw_torque_body` passed from control readout to the orientation RHS | May itself be an RPO held command; it is the pre-pseudoinverse request, not proof of realized/allocated torque or total spacecraft torque. |
| Prescribed reaction torque | Prescribed model's actual force/torque boundary | Includes transport term; not interchangeable with a wheel motor torque vector. |
| Demo reaction-torque magnitude | Scalar norm in CYGNSS demo | Magnitude only; cannot supply signed components or a frame transformation. |
| RPO held control torque | `rpo_control_torque_body` | Held control command, not total inertial torque. |
| Total torque / inertial inertia | Legacy `tau_ii_1..3`, `rw_torque_ii_1..3`, `J_ii_1..9` | No matching generic time-aligned saved producer or tensor flattening contract is established. Do not substitute static model inertia, body components or zeros. |

The retained wheel momentum renderer expects nested
`solution.physical_properties.rw_h` with one history per wheel. The retained
per-wheel torque renderer expects unprefixed `rw_tau_i` columns and the legacy
model's wheel count. Both use seconds and preserve trace order and signs. Their
nine-entry palette indexes zero on the ninth wheel, so nine-wheel input is a
known compatibility failure, not supported by the retained renderer. The total
torque and inertia renderers also take their unprefixed legacy columns directly;
their names and labels alone do not prove a producer's frame or matrix ordering.

**Disposition:** the opt-in integrated momentum schema and existing prescribed
aggregate schema are separately usable. A maintained wheel plot can follow this
contract in a later bounded packet. Per-wheel torque and inertial torque/inertia
plots need save-time producers and explicit frame/tensor metadata first. This
packet preserves the allocator, dynamics and retained renderers.

## Costate and switching contracts

The retained `lambda_plot` draws two independent trace pairs directly from
runtime configuration: `time_switch_list` / `lambda_switch_list` and
`time_list` / `lamv_list`. Its seconds and lambda labels perform no conversion,
shape checking, sorting, normalization or passage segmentation. Lambda units
and scaling are producer-specific; the labels alone do not establish them. The
removed legacy T-EDG formula kept the sign of `CD_slope`, while the typed
threshold uses its absolute value, so histories from the two cannot be silently
aliased. The runtime constructor initializes all four as empty `Float64` vectors.

On the inspected main revision, no tracked code writes or resets any of the four
arrays. Empty arrays are absence of recorded costates, not a zero costate
solution. The two legacy producers, the T-EDG plotting producer and the
constraint-tracking producer, were removed with the older E-EDG/T-EDG strategy
routes. Neither had a tracked caller, and the second failed on its own pushes. Their arrays carried no spacecraft or
passage identity, so a history recorded by older code is an unsegmented append
log in elapsed simulation seconds.

The typed EDG strategy has separate, call-local predicted lambda/alpha data and
an optional saved targeting switch time. The Odyssey example's
`targeting_switch_save_field` records `state.targeting_switch_s` in
`scN_edg_targeting_switch_s`, with `time` giving the saved simulation sample time and `N` the spacecraft list
position. Targeting brackets candidate switch times between the current elapsed
simulation time and that time plus the predicted duration (with its existing
one-second upper margin). Separately, typed heat-load prediction tracks use
relative times, and SpaceAGORAEDG's `_edg_low_alpha_switch_window` adds the current simulation
time to its selected track times. That heat-load window is not the saved
targeting field's producer.
The field is a **predicted switch timestamp in elapsed simulation seconds**, not lambda, time-to-go, an actual
switch event or a costate history. Nonfinite reset/no-prediction sentinels must
remain distinguishable from finite predictions; never replace them with zero.
A downstream diagnostic should label those samples unavailable and use a
step/held-prediction interpretation only when that matches its actual update
lifecycle. On pass reset, `_edg_reset_pass!` resets the affected spacecraft
`targeting_switch_s` to `Inf`, `heat_load_switches_s` to `(Inf, Inf)` and
`last_switch_solve_t` to `-Inf`; it retains cumulative solve counts and last
telemetry. `_edg_initialize_state!` resets the whole state, including those
counts and telemetry. These lifetimes are exercised through the actual reset
functions and saver. An explicit passage and spacecraft mapping is still needed
before combining legacy solver samples with typed targeting predictions.

**Disposition:** retain the legacy lambda renderer and explicitly classify its
costate input as unpopulated on main. Use the typed targeting-switch saved field
only for its documented prediction diagnostic. Reconstructing a costate history
requires a supported producer, costate definition/scaling, independent variable,
passage identity and reset policy. Existing active EDG owner work remains the
place for those scientific decisions; this cleanup does not synthesize histories
or change targeting behavior.

## Other retained capabilities

| Legacy diagnostic | Located producer or input | Remaining capability boundary |
| --- | --- | --- |
| `quaternion_plot` | Default `scN_q_1..4` for orientation runs; maintained raw-quaternion helper | Replacement traces are available; legacy retirement still requires a consumer decision. |
| `angular_velocity_plot` | Explicit `scN_body_rate_1..3` above; maintained body-rate helper | Default runs do not contain these fields; no schema guessing. |
| `drag_passage_plot`, `performance_plots` | Modern heat maxima, orbit and mission plots cover only some panels | Per-link heat histories, limits, energy/fuel and orbit grouping are not equivalent to maxima. |
| `angle_of_attack_plot`, `wind_relative_attitude_plot` | Some mission-specific control/link fields | Generic per-link AoA/sideslip and normalized-orbit schemas remain undecided; wind-relative code also references free variable `m`. |
| `closed_form_solution_plot` | Nested guidance histories | Numerical/closed-form pairing, free `args`, units and passage segmentation need a separate guidance contract. |
| `traj_2D`, `traj_3D`, `ABM_periapsis`, `ground_track` | Maintained trajectory, geography and extrema plots cover parts | Orbital-plane transforms, passage surfaces, maneuver-magnitude axes and time-versus-orbit traces retain distinct capabilities. |
| `plots` | Legacy model/options/nested-solution dispatcher | No complete supported replacement dispatcher is established. |

These contracts distinguish current producer behavior from an approved
mathematical correction or a maintained plotting replacement. They resolve the
Euler, wheel/torque and costate/switching ownership questions without declaring
all legacy diagnostics supported or retiring any of the 19 functions.
