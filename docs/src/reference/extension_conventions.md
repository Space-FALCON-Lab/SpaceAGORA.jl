# Extension Conventions

Every extension hook follows the same rules: one argument order, one name for
each kind of state, stated frames, SI units and positive consumption. This page
states each rule once. The [Extension Hook Reference](../generated/extension_hook_reference.md)
lists every hook with its shape, signature and fallback, and
[Migrating Extension Hooks](migrating_extension_hooks.md) maps the old names
to their replacements.

Use these rules when you write a method of a hook, and when you review a new
public name.

## Hook shapes

Every hook has one of eight shapes. The shape fixes where the hook runs and the
order of its arguments.

| Shape | Signature | Where it runs | Examples |
|---|---|---|---|
| T, tick | `update_<stage>!(model, u, p, t, spacecraft_index) -> nothing` | periodic and event callbacks, between integration steps | `update_control!`, `update_guidance!`, `update_navigation!` |
| R, right-hand side | `<quantity>(model, x, context, t) -> value` | every right-hand-side evaluation, one spacecraft | `wrench`, `control_wrench`, `control_propellant_rate`, `control_wheel_torque`, `gravity_backbone_acceleration_ii` |
| Q, query | `<quantity>(model, query...[, t][, flags...][, p])` | environment queries, no spacecraft | `getDensity`, `getDensityBatch!`, `getHeatRate`, `terrain_height`, `ephemerides_frame_lpi` |
| D, declaration | `<property>(model) -> value` | before the run, at setup | `environment_requirements`, `solver_partition`, `actuation_mode`, `actuator_limits`, `tick_schedule`, `planner_capabilities` |
| A, accessor | `<quantity>(model, spacecraft_index)` | at save time and in engine bookkeeping | `actuator_levels`, `control_thruster_levels`, `touchdown_spec`, `bind_spacecraft` |
| H, hardware | `actuator_<quantity>(actuator, level, ctx)`; `sensor_measure(sensor, truth, ctx, rng)` | inside the mounted wrench; at sensor samples | `actuator_wrench`, `actuator_consumption`, `actuator_saturate`, `sensor_measure` |
| P, planner or algorithm | `verb!(state, algorithm, request, rng)` | planning ticks | `plan_rpo!`, `retime_rpo!`, `propose!` |
| L, configuration | `verb(model, args; keywords...)` | before the run | `with_density_model_epoch`, `initialize_planner`, `gnc_callbacks` |

## Argument order

```julia
name(mutated, model, state, context, t, spacecraft_index, flags...; keywords...)
```

1. **The mutated argument first:** output buffers (`getDensityBatch!`) or run
   state that is not the model (a planner's `state`). When the model itself
   changes, as in a tick, the model is first.
2. **The model:** the object whose type your method dispatches on.
3. **The state:** `u` (all spacecraft) or `x` (one spacecraft). See below.
4. **The context:** `p` for GNC hooks, `env` for physical effectors, `ctx` for
   hardware. A query carries its query inputs here and keeps an optional `p`
   last.
5. **The time `t`.**
6. **The spacecraft index `spacecraft_index`,** only when the state is `u`.
7. **Flags** such as `wind::Bool`, then keywords.

**Time before index** is the part to remember. A hook whose order changes always
gets a new name, so an old method never receives its arguments in a new order.
Inside a versioned facade, a major contract version stands in for the new name.

## The state a hook receives

- **`u`** is the whole integrated state. `u.sc[spacecraft_index]` selects one
  spacecraft, with `pos` (m), `vel` (m/s) and `mass` (kg), and with attitude
  propagation `q` and `ω`. Among the model hooks, only the GNC ticks receive
  `u`. A GNC algorithm binding also receives it, in
  `build_algorithm_request(binding, u, p, t, reason, request_id)`: the binding's
  tick passes it on, so that the request is built from the state at that tick.
- **`x::StateSample`** is one spacecraft's state, for every right-hand-side
  hook:

  | Field | Meaning |
  |---|---|
  | `pos_ii`, `vel_ii` | position (m) and velocity (m/s) relative to the planet's centre, inertial frame |
  | `mass_kg` | the spacecraft's mass |
  | `q_ib` | the attitude quaternion, or `nothing` without attitude propagation |
  | `ω_body` | the body angular velocity in rad/s, or `nothing` without attitude propagation |
  | `spacecraft` | the spacecraft model, for its geometry and inertia |
  | `spacecraft_index` | the spacecraft's position in the run; 0 for a sample built outside a run |

- The raw per-spacecraft view `u.sc[i]` is not passed to any new hook.

### The supported sample argument

Type the sample argument `x::StateSample`, or leave it untyped; the same holds
for `env::EnvironmentSample`. The engine looks your method up with these
abstract types. A method typed on a concrete `StateSample{S}`, for example
`x::StateSample{SpacecraftModel}`, is not routed to. The pre-run check then
refuses the run: it names the hook, and says that the model has methods of that
function, but none accepts the arguments the run passes. The templates and the
walkthroughs use the supported form.

Ticks are looked up with the run's own state and parameter types. Type `u` as
`u`, `u::AbstractVector` or `u::ComponentVector`; a method typed
`u::Vector{Float64}` misses the state the engine passes and is refused before
the run in the same way.

## The spacecraft index

- `spacecraft_index` is the spacecraft's position in the run, from 1, in the
  order of `dynamics_model.spacecraft`. It is an argument of the ticks and the
  accessors, and the field `x.spacecraft_index` in the right-hand-side hooks.
- `*_id` names the user's stable identifier, such as `SpacecraftModel(id=...)`,
  or `chaser_id` and `target_id` on an RPO request. An index and an identifier
  are never interchangeable.
- Models that select their vehicle by a field use the same words. The built-in
  models' fields `sat_idx`, `spacecraft_idx`, `chaser_idx` and `target_idx`
  became `spacecraft_index`, `chaser_index` and `target_index` in 0.3.0, and the
  `spacecraft_idx` of `RobotArmReactionEffector` in 0.3.1; the old names were
  removed in 0.5.0, and the engine reads only the new ones on a model of your
  own. `StateAnchor`'s `sat_idx` and the `sat_idx` keyword of
  `get_periapsis_pulse_callback` became `spacecraft_index` in 0.5.0, with a
  warning for the old name until 1.0. `SchedulingPolicyModel`'s `target_idx` and
  `target_idxs` became `target_index` and `target_indices` in 0.5.0, and all four
  were removed in 0.6.0 with the scheduling policies chosen by a symbol: a policy
  type holds its own targets (see the [migration page](migrating_extension_hooks.md)).
- In the external-propagation facade the owner-local position of a spacecraft is
  called `slot`, because it counts within the owner's list, not the run.

## Context

- **GNC hooks** (ticks and control outputs) take `p`, the run's parameters;
  `p.args` is the `SimulationConfiguration`. There is no separate GNC context
  object.
- **Physical effectors** take `env::EnvironmentSample`: the planet, the
  planet-fixed frame, the atmosphere, the Sun and the third bodies, each filled
  when the model's `environment_requirements` asks for it.
- **Actuators** take `ctx::ActuatorContext`: the mounting link's pose in the body
  frame, the vehicle's mass, attitude and rate, the magnetic field and the
  sampled environment the kind requests.
- **Sensors** take `truth::SensorTruth` (inertial position and velocity,
  attitude, rate, mass, time and the field) and `ctx`.
- **Planners and algorithms** take a `request` and an `rng`.

## Time

`t` is the elapsed time since the run's start epoch, in seconds (`Float64`).
`getDensity`'s `el_time` is the same elapsed time. Only the ephemerides hooks take
`et`, ephemeris seconds.

## Frames

Names carry their frame as a suffix: `_ii` inertial, `_body` body (the bus
axes), `_pp` planet-fixed, `_rtn` the target's radial, along-track and
cross-track frame, `_site` a landing site's frame.

| Level | Hooks | Force | Torque |
|---|---|---|---|
| Vehicle | `wrench`, `control_wrench` | inertial frame, N | body frame, about the centre of mass, N·m |
| Hardware | `actuator_wrench` | body frame, N; the engine rotates it to the inertial frame | body frame, about the body reference point, N·m |
| Wheels | `control_wheel_torque` | none | body frame, N·m |
| Interlinks | `interlink_wrench` | inertial frame, N | body frame, about the body reference point, N·m |

The attitude quaternion is `q = (x, y, z, w)`, scalar last, with
`v_body = rot(q) * v_inertial`.

## Units

- SI everywhere: metres, seconds, kilograms, newtons, radians, watts. A rate is
  per second.
- A unit that is not obvious goes in the name: `alt_m`, `thrust_n`, `lat_deg`.
- Degrees appear only where a run is configured, such as `InitialCondition` or
  `ApolloDescentConfig(site_lat_deg=...)`.

Three documented exceptions stay until 1.0 or later:

- the terrain hooks take planetocentric latitude and longitude in degrees
  (`terrain_height(model, lat_deg, lon_deg)`);
- `getHeatRate` returns W/cm²;
- `getDensity` takes geodetic altitude in metres and geodetic latitude in
  radians.

## Signs

- **Consumption is positive.** `control_propellant_rate` returns the propellant
  used in kg/s, and `actuator_consumption` the mass rate and power used. The
  engine owns the sign of the mass derivative.
- A thruster's direction is the force it puts on the vehicle.
- Wind is the air's velocity relative to the rotating planet, as east, north and
  up components.
- **The sign trap.** `calcControlMassFlowRate`, removed in 0.5.0, returned the
  mass derivative, negative while propellant is used. A method renamed to
  `control_propellant_rate` without flipping its sign makes the vehicle gain
  mass.

## Naming

- Names are snake_case, with `!` exactly when the hook mutates an argument or
  the model's own state.
- New names have no `calc` or `get` prefix: the name says what is returned
  (`control_wrench`) or what happens (`update_control!`).
- A family prefix groups related hooks: `control_`, `actuator_`, `sensor_`,
  `gravity_backbone_`, `terrain_`, `planner_`.
- A hook that returns one vector in a fixed frame ends with the frame suffix
  (`gravity_backbone_acceleration_ii`).
- No export takes a name that packages often loaded beside SpaceAGORA export: a
  plain `density` would collide with `Plots.density`. A test checks the hook
  names against Plots, OrdinaryDiffEq, ComponentArrays, StaticArrays and
  DataFrames.
- `getDensity`, `getDensityBatch!` and `getHeatRate` keep their names until 1.0.

## Fallbacks and method typing

- A hook's default method types nothing, or only the model, with an abstract type
  you subtype. A method you write for your own model type is therefore more
  specific and never ambiguous with the default.
- A method you write never needs a type annotation only to avoid an ambiguity.
- The ticks and the right-hand-side hooks have no catch-all method, except
  `control_propellant_rate` and `control_wheel_torque`, whose untyped defaults
  return `0.0` and `nothing`: an effector that uses no propellant, or drives no
  wheel, needs no method of them. The pre-run check refuses a model without a
  method of a hook the run calls.

## Deprecations

When a hook is renamed or retired, its old name keeps working for a deprecation
window:

- **Before the run.** The pre-run check warns once per model type and old name,
  before integration, never during the solve, and names the replacement.

- **Direct calls.** Calling an old name on a model that implements the new hook
  runs the new hook and warns once per type.
- **As an error.** Under `julia --depwarn=error` these warnings are errors: the
  run stops before integration, and a direct call throws. CI runs the examples
  and the probe scripts this way.
- **Removal.** A deprecated name is removed in a later release, and not before a
  tagged release has carried the warnings and the companion packages have moved
  to the new names. SpaceAGORA 0.5.0 removed the hook names 0.3.0 deprecated.
- **Version numbers** follow `version` in `Project.toml`. The older `v1.0.0` and
  `v2.0.0` git tags predate it.

The other deprecations of the 0.3.0 cycle, the old index field names,
`RPOGuidanceModel` and selecting a comparison planner by name, were removed in
0.5.0. [Migrating Extension Hooks](migrating_extension_hooks.md) has the full
table and the schedule.

## Versioned facades

`EDGServices`, `ExternalPropagationServices` and `HYPRServices` are named,
versioned boundaries for optional packages. A companion checks the facade's
`CONTRACT_VERSION` exactly, and a major version marks a change of any
signature. `EDGServices` and `HYPRServices` are at contract 3.0.0 and
`ExternalPropagationServices` at 2.0.0; all three follow these conventions.
Names in a facade's inventory keep their exact spelling, underscores included.
See the [Optional Adapter Standard](adapter_contract.md) and each facade's page.

## Run state and randomness

- Keep run state in the model, or in the state your algorithm declares, never in
  globals.
- Every model must survive `deepcopy`: `run_simulation` runs a copy of the
  configuration unless `isolate_state=false`.
- A model declared on one spacecraft implements
  `bind_spacecraft(model, spacecraft_index)` to return a copy for that
  spacecraft, with run state of its own.
- Draw random numbers only from the `rng` you are given. Planning functions take
  an explicit stream; there is no default one.
