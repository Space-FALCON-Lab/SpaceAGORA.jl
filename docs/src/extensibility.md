# Extensibility

Use this page when you are adding new models or control hooks to SpaceAGORA.

This page is for extension authors who want to stay on the stable root package
surface instead of binding new work to package implementation details.

Shortest successful starting point:

```text
templates/force_torque_model_template.jl
```

What to read next:

- [Adding a Force or Torque of Your Own](user/custom_effector.md) (a complete, tested walkthrough of the force/torque hook)
- [Adding an Actuator of Your Own](user/custom_actuator.md) (a complete, tested walkthrough of the actuator contract)
- [Public API](generated/public_api.md)
- [Concepts](user/concepts.md)

SpaceAGORA now exposes a minimal stable extension contract from the root package.

The rule is simple:

1. subtype a stable abstract interface exported by `SpaceAGORA`
2. implement the matching root hook function
3. wire the model into the typed runtime configuration

Do not build new extensions by importing internal modules directly unless you
are intentionally working on SpaceAGORA internals.

Optional packages with a broader integration use a named, versioned facade.
See the [Optional Adapter Standard](reference/adapter_contract.md) for shared
ownership, compatibility, installation and conformance requirements, and
[External Propagation Services](reference/external_propagation_services.md) for
external spacecraft integrators. HYPR uses its existing
[versioned planning services](reference/hypr_services.md).

## Add a force or torque model

Stable interface:

- `SpaceAGORA.AbstractForceTorqueModel`
- `SpaceAGORA.wrench`
- `SpaceAGORA.environment_requirements`
- `SpaceAGORA.solver_partition`
- `SpaceAGORA.gravity_backbone_structure`
- `SpaceAGORA.gravity_backbone_acceleration_ii`
- `SpaceAGORA.gravity_backbone_kick_structure`
- `SpaceAGORA.gravity_backbone_kick_acceleration_ii`
- `SpaceAGORA.calcForceTorque` (deprecated, removed in 0.5.0)

Preferred additive methods:

```julia
SpaceAGORA.environment_requirements(model) -> SpaceAGORA.EffectorEnvironmentRequirements
SpaceAGORA.solver_partition(model) -> :explicit | :implicit
SpaceAGORA.gravity_backbone_structure(model) -> :unsupported | :position_only_static_gravity
SpaceAGORA.gravity_backbone_acceleration_ii(model, x::SpaceAGORA.StateSample, env::SpaceAGORA.EnvironmentSample, t::Float64) -> accel_ii
SpaceAGORA.gravity_backbone_kick_structure(model) -> :unsupported | :velocity_kick_explicit
SpaceAGORA.gravity_backbone_kick_acceleration_ii(model, x::SpaceAGORA.StateSample, env::SpaceAGORA.EnvironmentSample, t::Float64) -> accel_ii
SpaceAGORA.wrench(model, x::SpaceAGORA.StateSample, env::SpaceAGORA.EnvironmentSample, t::Float64) -> (force_ii, torque_body)
```

Deprecated method, removed in SpaceAGORA 0.5.0:

```julia
SpaceAGORA.calcForceTorque(model, x, p, i) -> (force_n, torque_n_m)
```

Implement `wrench`. The engine owns stage-consistent sampling and passes typed
state/environment data into the effector; `x.spacecraft_index` is the
spacecraft's position in the run. A model with only `calcForceTorque` still runs
through it, and the run warns before integration (an error under
`julia --depwarn=error`). Calling `calcForceTorque` on a model that implements
`wrench`, as every built-in model does, returns that model's `wrench` for the
state `x` it is given.
Override `solver_partition` only if the effector should move to the implicit
side of `split_imex`; the default is `:explicit`.
Override the gravity-backbone hooks only if the effector should participate in
`gravity_backbone_split`, either as a strict position-only static-gravity core
term or as an explicit translational velocity kick.

Solver partition contract:

- `split_imex` is the atmosphere-implicit IMEX path
- `multirate` remains the control-focused split path
- `gravity_backbone_split` is a fixed-step gravity-backbone split mode with a symplectic gravity core plus explicit SRP / N-body velocity kicks
- `gravity_backbone_split` remains translational-only; it is not a fully symplectic whole-system solve
- force returned by `wrench` is inertial-frame and torque is body-frame
- the current heat-load state is an accumulated heat-rate integral and stays explicit

Registration:

- add the model instance to `dynamic_effectors` in the simulation configuration

Template:

- `templates/force_torque_model_template.jl`

## Add a density model

Stable interface:

- `SpaceAGORA.AbstractDensityModel`
- `SpaceAGORA.getDensity`
- `SpaceAGORA.getDensityBatch!`

Required scalar method:

```julia
SpaceAGORA.getDensity(model, h, lat, lon, el_time, wind[, p]) -> (rho, temperature, wind_vec)
```

`h` is height in metres, `lat` and `lon` are in radians, and `el_time` is elapsed
seconds from the scenario epoch. `wind` is a `Bool` saying whether the caller
needs the wind vector: the force and heating paths pass `true`, and density-only
callers such as the visualization pass `false`. A model may return a zero
`wind_vec` when `wind` is `false`. `wind_vec` holds the east, north and up wind
components in m/s.

Optional batch method:

```julia
SpaceAGORA.getDensityBatch!(rhos, Ts, winds, model, hs, lats, lons, el_time, wind, p)
```

Registration:

- set the model as `environment_model.density_model`

Template:

- `templates/density_model_template.jl`

## Add a control effector or hook

Stable interface:

- `SpaceAGORA.AbstractControlEffectorModel`
- `SpaceAGORA.calcControlEffect!`
- `SpaceAGORA.control_wrench`
- `SpaceAGORA.control_propellant_rate`
- `SpaceAGORA.control_wheel_torque`
- `SpaceAGORA.calcControlForceTorque`, `SpaceAGORA.calcControlMassFlowRate` and
  `SpaceAGORA.calcReactionWheelTorque` (deprecated, removed in 0.5.0)

Typical methods:

```julia
SpaceAGORA.calcControlEffect!(model, u, p, t, spacecraft_index)
SpaceAGORA.control_wrench(model, x::SpaceAGORA.StateSample, p, t::Float64) -> (force_ii, torque_body)
SpaceAGORA.control_propellant_rate(model, x::SpaceAGORA.StateSample, p, t::Float64) -> kg/s, positive
SpaceAGORA.control_wheel_torque(model, x::SpaceAGORA.StateSample, p, t::Float64) -> nothing or torque_body
```

`calcControlEffect!` is the control tick: it computes the command and holds it.
The three output hooks return the held command at every right-hand-side
evaluation, and `x.spacecraft_index` is the spacecraft's position in the run.
An effector that uses no propellant, or drives no wheel, leaves out
`control_propellant_rate` or `control_wheel_torque`.

The deprecated output hooks, removed in SpaceAGORA 0.5.0:

| Deprecated | Replacement | What changes |
|---|---|---|
| `calcControlForceTorque(model, u, p, i, t)` | `control_wrench(model, x, p, t)` | a state sample, the time before the index, the index in `x.spacecraft_index` |
| `calcControlMassFlowRate(model, u, p, i, t)` | `control_propellant_rate(model, x, p, t)` | the same, and **the sign**: the new hook returns the propellant used, positive, and the old one the mass derivative, negative |
| `calcReactionWheelTorque(model, u, p, i, t)` | `control_wheel_torque(model, x, p, t)` | the same as for the wrench |

An effector with only the old hooks still runs through them, and the run warns
before integration (an error under `julia --depwarn=error`). Calling an old
name on an effector that implements its replacement, as every built-in
effector does, returns the replacement's result, the mass rate negated, and
warns once per effector type.

Registration:

- add the model instance to `ControlModel(control_effectors=(...))`

Template:

- `templates/control_hook_template.jl`

### Per-spacecraft declaration

A `SpacecraftModel` can carry its own GNC (`SpacecraftModel(...; control=ControlModel(...))`,
likewise `guidance` and `navigation`). At the start of `run_simulation` each such effector is
passed through `SpaceAGORA.bind_spacecraft(effector, sat_idx)` with the position of its
spacecraft and appended to the configuration-level tuples. An effector that selects its
vehicle by an index field implements it by returning a copy with that field set:

```julia
struct MyPush <: SpaceAGORA.AbstractControlEffectorModel
    sat_idx::Int
    force_n::SVector{3, Float64}
end

SpaceAGORA.bind_spacecraft(m::MyPush, sat_idx::Int) = MyPush(sat_idx, m.force_n)
```

The fallback throws an `ArgumentError`: an effector without an index acts on every
spacecraft, and an effector that couples several vehicles (the RPO chaser and target
types) cannot belong to one, so both stay at configuration level. Built in: the
momentum manager, the robot-arm control effector and the Apollo descent guidance and control
models bind.

## Add an actuator

An actuator is hardware mounted on a spacecraft link that acts on the vehicle
when commanded. You describe the hardware once; the engine derives the force,
the torque, the propellant use and the result columns (each actuator's applied
level, the actuators' power draw and the thrusters' firing levels) from the
actuators mounted on the vehicle.

Stable interface:

- `SpaceAGORA.AbstractActuator`
- `SpaceAGORA.actuator_wrench` and `SpaceAGORA.actuator_limits` (required)
- `SpaceAGORA.actuator_consumption`, `SpaceAGORA.actuator_saturate`,
  `SpaceAGORA.actuator_fraction`, `SpaceAGORA.actuator_unit_wrench`,
  `SpaceAGORA.actuator_requirements` (optional, with defaults)
- `SpaceAGORA.add_actuator!`, `SpaceAGORA.mounted_actuators`,
  `SpaceAGORA.allocation_matrix`, `SpaceAGORA.allocate_actuators`

Typical methods:

```julia
SpaceAGORA.actuator_limits(a::MyActuator) = (lower, upper)          # level bounds, the kind's unit
SpaceAGORA.actuator_wrench(a::MyActuator, level, ctx) = (force_body, torque_body)
SpaceAGORA.actuator_consumption(a::MyActuator, level, ctx) = (mass_rate_kg_s = ..., power_w = ...)
```

The wrench is in the spacecraft body frame, the torque about the body
reference point; `ctx` (an `ActuatorContext`) gives the mounting link's pose
in the body frame and the vehicle's mass, attitude and rate.

A kind that needs more of its surroundings declares it with
`SpaceAGORA.actuator_requirements(a) = ActuatorRequirements(...)`: the
planet-fixed frame, the atmosphere, the Sun's position or the local magnetic
field. The engine samples what the mounted actuators ask for once per
evaluation and passes it in the context: the field as
`ctx.magnetic_field_body`, the rest as `ctx.environment`, an
`EnvironmentSample` (for example `ctx.environment.atmosphere.rho_kg_m3`).
Without a request, `ctx.environment` is `nothing`.

The built-in kinds:

- `Thruster`: its level is a thrust in newtons and its `direction` is the
  thrust direction. With `firing_dynamics=true` the command goes through the
  thruster's minimum on-time, Schmitt trigger and thrust ramp at each control
  tick, and the vehicle feels the resulting average thrust until the next one.
- `ReactionWheel`: its level is the torque on the body about its spin axis; it
  mounts on the root link and its momentum is integrated, with torque and
  momentum limits. Given a spin inertia, a base power and a motor efficiency,
  it reports the power it draws.
- `Magnetorquer`: its level is the commanded dipole along its axis; it acts in
  the run's magnetic field. Given its power at full dipole, it reports a power
  that grows with the square of the dipole.
- `LaserTerminal`: its level is the transmitted optical power; the vehicle
  recoils with power over the speed of light, and it draws power over its
  efficiency.
- `PanelHinge`: its level is the angle of the non-root link it turns, with
  angle and rate limits; the engine applies it at the commanding effector's
  control ticks. The energy-depletion controller also moves its panels through
  a hinge mounted on them.
- `ThrusterGimbal`: one axis of a gimbal on a thruster mounted on the same
  link; its level is the gimbal angle, with angle and rate limits, applied at
  control ticks like a hinge's. The thrust direction turns about the gimbal
  axis (the nozzle stays where it is), so the thruster's force, torque,
  propellant use and allocation column follow it; two gimbals on one
  thruster make a two-axis gimbal, turning in mounting order.
  `thrust_direction(thruster)` gives the turned direction.
- `ControlMomentGyro`: a single-gimbal CMG on the root link, a rotor of
  constant momentum whose spin axis turns about the gimbal axis; its level is
  the gimbal rate, within its rate limit, and its gimbal stops at
  `min_gimbal_angle_rad`/`max_gimbal_angle_rad` (none by default). The gimbal angle is integrated with
  the vehicle's state (`cmg_angles`, only when a CMG is mounted); the vehicle
  feels `-δ̇ (g × h)` and the rotor's momentum in its gyroscopic term.
  `cmg_steering(spacecraft, angles, torque)` turns a desired torque into
  gimbal rates with a singularity-robust inverse. The gravity-backbone mode
  and articulated spacecraft refuse CMGs.

The actuators' power feeds the spacecraft's power budget when it carries one:
solar arrays, a battery and constant loads mounted with `add_power!`. See
[Powering a Spacecraft](user/power_budget.md).

A mounted actuator can fail: `add_fault!(spacecraft, ActuatorFault(actuator;
kind, start_s, stop_s, value))` makes it apply level zero (`:failed`; a failed
hinge or gimbal holds its angle), hold a level (`:stuck`, by default the one it was
commanded to when the fault began), act at a fraction of its level and limits
(`:degraded`) or with an offset (`:biased`) between `start_s` and `stop_s`. The
fault switches exactly at those times, and acts on what the vehicle feels, the
propellant, the power and the `actuator_level` columns; the controller keeps
seeing its own commands. On a thruster with firing dynamics the fault acts on
the thrust it delivers, so it takes effect at its start, between control
ticks. Runs with faults write `sc{i}_fault_active_{j}`. A
checkpointed run keeps the faults' state from one segment to the next, but a
`:stuck` fault without a `value` refuses `resume_from_checkpoint`: the
checkpoint file does not hold the level it captured at onset.

Registration: `add_actuator!(spacecraft, link, actuator)`.

A control effector drives the mounted actuators by declaring
`SpaceAGORA.actuation_mode(c) = :mounted` and returning one level per mounted
actuator, in mounting order, from `SpaceAGORA.actuator_levels(c, i)`; it
computes the levels in `calcControlEffect!` and writes no force or mass-flow
hook. Effectors in the default `:ideal` mode keep applying their own hooks.
The built-in RPO and Apollo descent controllers and the magnetic momentum
manager accept `actuation=:mounted`. So does the orbit-maneuver burn scheduler
`BaseThrusterModel`: it fires a mounted engine during the burn window and
points it along the burn with the mounted reaction wheels; the engine switches
exactly at the window's ends, between integration steps. A run with a mounted
effector writes the `sc{i}_actuator_level_{k}` and `sc{i}_actuator_power` columns (see
[Output Columns](user/output_columns.md#Mounted-actuators)); `applied_actuator_levels` and
`mounted_power` compute the same values from a spacecraft and its levels.

Template:

- `templates/actuator_template.jl`

Walkthrough: [Adding an Actuator of Your Own](user/custom_actuator.md).
Example: `examples/Earth_CubeSat_Mounted_Actuators.jl`, a CubeSat whose
mounted wheels detumble it while its torque rods unload them.

## Add a sensor

A sensor measures the state of the vehicle for navigation; it does not act on
the dynamics. A sensor kind is a mutable struct subtyping `AbstractSensor`
with its parameters and five fields the engine fills or reads: `period_s`
(the sampling period, finite and positive), `seed` (selects its random
stream), `rng`, `latest` (the measurement's components as a vector, or a
number for a one-component sensor) and `latest_time_s`. It implements:

- `sensor_measure(sensor, truth::SensorTruth, ctx::ActuatorContext, rng)`:
  one measurement of `truth` (the inertial position and velocity, the
  attitude, the body rate, the mass, the time and, when the sensor asks for
  it, the magnetic field in the body frame), with any noise drawn from `rng`.
  `ctx.link_to_body` maps the sensor's frame, its link's, to the body frame.
- `sensor_dimension(sensor)`: the number of components.
- optionally `sensor_requirements(sensor)`: `ActuatorRequirements(magnetic_field=true)`
  to receive the field (from the run's `MagneticTorqueRodModel`, or the
  tilted dipole); no other environment is available to sensors.

Mount it with `add_sensor!(spacecraft, link, sensor)`; a sensor on a hinged
panel turns with it. The built-in kinds are `RateGyro` (body rate with white
noise, a bias and a bias random walk), `StarTracker` (attitude with cross- and
about-boresight noise, boresight along its link's +z) and `Magnetometer`
(field with white noise and a bias).

The engine samples every sensor at the run's start and every period after
it, from the accepted state; the samples are solver stops, and a checkpointed
run samples at the same times as a run in one piece. Each sensor's random
stream comes from its `seed`, the spacecraft's ID and its mounting position,
so a run repeats with the same seeds. Seeds are not varied for you: the
built-in kinds default to `seed=1`, so in a Monte Carlo campaign set each
sensor's `seed` from the sample's seed where the sample's spacecraft is
built, or the samples share one noise stream. A navigation effector reads the latest measurement in
`calcNavigationEffect!` with `sensor_measurement(p.args.dynamics_model.spacecraft[i], k)`
(and its time with `sensor_measurement_time`); sensors sample before
navigation effectors at the same time. Runs with sensors write the
`sc{i}_sensor_{k}` columns (see [Output Columns](user/output_columns.md#Sensors))
and cannot resume from a checkpoint file; the gravity-backbone and multirate
solver modes refuse mounted sensors.

## Planet, ephemerides, thermal, thruster, and guidance extensions

Stable interfaces:

- `SpaceAGORA.AbstractPlanet`
- `SpaceAGORA.AbstractEphemeridesModel`
- `SpaceAGORA.AbstractThermalModel`
- `SpaceAGORA.AbstractThrusterModel`
- `SpaceAGORA.AbstractGuidanceModel`

These interfaces are public because no-GRAM onboarding and future user-defined
body/back-end, thermal, thruster, and guidance work need a stable contract,
but the practical extension recipes for these types are still thinner than
the force/density/control paths above: there is no single dedicated root
hook function to implement (thermal models are consumed inside the
force/torque and heat-load machinery, thruster models inside control
effectors, and guidance models inside the GNC control hooks). Treat them as
advanced extension points — subtype the abstract type and follow an existing
implementation in `src/models/vehicle/thermal/`, `src/models/vehicle/actuators/`,
or `src/models/gnc/guidance/` as a model rather than a documented recipe.

## Design guidance

When adding a new model:

1. keep state in the model object, not in global variables
2. make units explicit in field names or docstrings
3. keep expensive caches out of hot dispatch paths unless they are typed
4. add a small typed constructor and one smoke test
5. document the model through the root `SpaceAGORA` hook surface, not by telling
   users to import internal modules

## Templates and runnable examples

The four files under `templates/` define different extension contracts: a force
or torque, an atmosphere density model, a control hook and an actuator. They define model
types and methods; they do not construct or run a mission. Copy the relevant
file into your model's owning package and add tests for the behavior you implement.
The force/torque template implements `wrench`;
[Adding a Force or Torque of Your Own](user/custom_effector.md) walks through
one in detail.

Use the [Examples Catalog](user/examples_catalog.md) when you want a complete
scenario to run or adapt. Helpers under `examples/support/` serve those scenarios
and are not extension contracts.
