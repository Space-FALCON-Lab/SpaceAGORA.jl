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
- [Writing a GNC Algorithm](user/gnc_algorithm.md) (controllers, estimators and planners, with worked examples run in CI)
- [Building Your Own Research Repository](user/research_repository.md) (starters, conformance checks and CI for work kept in a repository of your own)
- [Public API](generated/public_api.md)
- [Concepts](user/concepts.md)

SpaceAGORA exposes a stable extension contract from the root package. The rule
is simple:

1. subtype a stable abstract interface exported by `SpaceAGORA`
2. implement the matching root hook function
3. wire the model into the typed runtime configuration

Every hook follows the same [Extension Conventions](reference/extension_conventions.md):
the argument order (time before the spacecraft index), which state a hook
receives, frames, units, signs and names. The
[Extension Hook Reference](generated/extension_hook_reference.md) lists each
hook with its shape, signature and fallback, and
[Migrating Extension Hooks](reference/migrating_extension_hooks.md) maps the
old names to their replacements. Each section below names its hooks'
shape and links its rows in the reference.

Do not build new extensions by importing internal modules directly unless you
are intentionally working on SpaceAGORA internals.

Optional packages with a broader integration use a named, versioned facade.
See the [Optional Adapter Standard](reference/adapter_contract.md) for shared
ownership, compatibility, installation and conformance requirements, and
[External Propagation Services](reference/external_propagation_services.md) for
external spacecraft integrators. HYPR uses its existing
[versioned planning services](reference/hypr_services.md).

## Add a force or torque model

Shape: right-hand side (R), plus declarations (D). Rows:
[Physical effectors](generated/extension_hook_reference.md#Physical-effectors).

Stable interface:

- `SpaceAGORA.AbstractForceTorqueModel`
- `SpaceAGORA.wrench`
- `SpaceAGORA.environment_requirements`
- `SpaceAGORA.solver_partition`
- `SpaceAGORA.gravity_backbone_structure`
- `SpaceAGORA.gravity_backbone_acceleration_ii`
- `SpaceAGORA.gravity_backbone_kick_structure`
- `SpaceAGORA.gravity_backbone_kick_acceleration_ii`

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

Implement `wrench`. The engine owns stage-consistent sampling and passes typed
state/environment data into the effector; `x.spacecraft_index` is the
spacecraft's position in the run. To have an effector act on some spacecraft
only, give its entry as `scoped(effector, scope)`
([Effectors on part of a constellation](user/constellations.md#Effectors-on-part-of-a-constellation)). Type the sample `x::SpaceAGORA.StateSample`,
or leave it untyped
([the supported sample argument](reference/extension_conventions.md#The-supported-sample-argument)).
Override `solver_partition` only if the effector should move to the implicit
side of `split_imex`; the default is `:explicit`.
Override the gravity-backbone hooks only if the effector should participate in
`gravity_backbone_split`, either as a strict position-only static-gravity core
term or as an explicit translational velocity kick.

Solver partition contract:

- `split_imex` is the atmosphere-implicit IMEX path
- `multirate` remains the control-focused split path
- `multirate_mri` uses a fixed coupled partition: custom dynamic effectors stay slow, and `solver_partition` does not change their MRI placement
- `gravity_backbone_split` is a fixed-step gravity-backbone split mode with a symplectic gravity core plus explicit SRP / N-body velocity kicks
- `gravity_backbone_split` remains translational-only; it is not a fully symplectic whole-system solve
- force returned by `wrench` is inertial-frame and torque is body-frame
- the current heat-load state is an accumulated heat-rate integral and stays explicit

Registration:

- add the model instance to `dynamic_effectors` in the simulation configuration

Template:

- `templates/force_torque_model_template.jl`

## Add a density model

Shape: query (Q). Rows:
[Environment queries](generated/extension_hook_reference.md#Environment-queries).

Stable interface:

- `SpaceAGORA.AbstractDensityModel`
- `SpaceAGORA.getDensity`
- `SpaceAGORA.getDensityBatch!`

Required scalar method:

```julia
SpaceAGORA.getDensity(model, h, lat, lon, el_time, wind[, p]) -> (rho, temperature, wind_vec)
```

`h` is geodetic altitude in metres, `lat` (geodetic) and `lon` are in radians,
and `el_time` is elapsed seconds from the scenario epoch. `wind` is a `Bool` saying whether the caller
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

Shape: tick (T) for `update_control!`, right-hand side (R) for the outputs.
Rows: [GNC ticks](generated/extension_hook_reference.md#GNC-ticks) and
[Control outputs](generated/extension_hook_reference.md#Control-outputs).

Stable interface:

- `SpaceAGORA.AbstractControlEffectorModel`
- `SpaceAGORA.update_control!`
- `SpaceAGORA.control_wrench`
- `SpaceAGORA.control_propellant_rate`
- `SpaceAGORA.control_wheel_torque`

Typical methods:

```julia
SpaceAGORA.update_control!(model, u, p, t, spacecraft_index)
SpaceAGORA.control_wrench(model, x::SpaceAGORA.StateSample, p, t::Float64) -> (force_ii, torque_body)
SpaceAGORA.control_propellant_rate(model, x::SpaceAGORA.StateSample, p, t::Float64) -> kg/s, positive
SpaceAGORA.control_wheel_torque(model, x::SpaceAGORA.StateSample, p, t::Float64) -> nothing or torque_body
```

`update_control!` is the control tick: it computes the command and holds it.
The three output hooks return the held command at every right-hand-side
evaluation, and `x.spacecraft_index` is the spacecraft's position in the run.
`control_propellant_rate` returns the propellant used, **positive**. An effector
that uses no propellant, or drives no wheel, leaves out
`control_propellant_rate` or `control_wheel_torque`: their defaults return
`0.0` and `nothing`.

A controller also declares what it acts on and what it assumes
(`spacecraft_scope`, `gnc_period_s`, `skip_inactive`, `held_command`), can
share an object with an estimator or a planner through `gnc_links`, and resets
its run state at each run's start in `initialize_control!`.
[Writing a GNC Algorithm](user/gnc_algorithm.md#The-feedback-profile) explains
these declarations with a worked controller run in CI.

Registration:

- add the model instance to `ControlModel(control_effectors=(...))`

Template:

- `templates/control_hook_template.jl`

### Per-spacecraft declaration

A `SpacecraftModel` can carry its own GNC (`SpacecraftModel(...; control=ControlModel(...))`,
likewise `guidance` and `navigation`). At the start of `run_simulation` each such effector is
passed through `SpaceAGORA.bind_spacecraft(effector, spacecraft_index)` with the position
of its spacecraft and appended to the configuration-level tuples. An effector that selects
its vehicle by an index field implements it by returning a copy with that field set, and
declares the spacecraft it acts on with `spacecraft_scope`, so the engine ticks it and
calls its hooks for that spacecraft only:

```julia
struct MyPush <: SpaceAGORA.AbstractControlEffectorModel
    spacecraft_index::Int
    force_n::SVector{3, Float64}
end

SpaceAGORA.bind_spacecraft(m::MyPush, spacecraft_index::Int) = MyPush(spacecraft_index, m.force_n)
SpaceAGORA.spacecraft_scope(m::MyPush) = BoundSpacecraft((m.spacecraft_index,))
```

A bound copy whose `spacecraft_scope` is still `AllSpacecraft()` acts on every spacecraft,
as at configuration level; each run warns once for each such type, naming the run position
and id of the spacecraft. The fallback throws an
`ArgumentError`: an effector without an index acts on every spacecraft, and an effector
that couples several vehicles (the RPO chaser and target types) cannot belong to one, so
both stay at configuration level. So does the burn scheduler, `BaseThrusterModel`, one
model with one entry per spacecraft. Built in: the momentum manager, the robot-arm
control effector, the Apollo descent guidance and control models, the algorithm guidance
binding (`AlgorithmGuidance`) and, from SpaceAGORA 0.7.0, the two maneuver guidance models
(`ApoapsisTargetPeriapsisRaiseGuidanceModel` and
`AerobrakingCampaignPropulsiveManeuverGuidanceModel`) bind.

## Add a guidance or navigation effector

Shape: tick (T). Rows:
[GNC ticks](generated/extension_hook_reference.md#GNC-ticks).

Stable interface:

- `SpaceAGORA.AbstractGuidanceModel` and `SpaceAGORA.update_guidance!`
- `SpaceAGORA.AbstractNavigationModel` and `SpaceAGORA.update_navigation!`

Typical methods:

```julia
SpaceAGORA.update_guidance!(model, u, p, t, spacecraft_index)
SpaceAGORA.update_navigation!(model, u, p, t, spacecraft_index)
```

The ticks take the arguments of `update_control!`. A navigation effector takes
its measurements or updates its estimate; a guidance effector updates a plan,
a target or a mode; each keeps the result in its own state for the control
effectors to read. When ticks coincide, navigation runs first, then control,
then guidance; energy depletion guidance with targeting runs before control
instead.

Registration:

- add guidance effectors to `GuidanceModel(guidance_effectors=(...), guidance_rates=[...])`
- add navigation effectors to `NavigationModel(navigation_effectors=(...), navigation_rates=[...])`

[Writing a GNC Algorithm](user/gnc_algorithm.md#Worked-example-B:-an-estimator-and-a-link)
has a worked estimator that a controller reads through a declared link.

A planner, whose output is used later or can be checked before use, follows the
planning profile instead: it proposes, and its family validates and installs.
See [Writing a GNC Algorithm](user/gnc_algorithm.md#The-planning-profile), the
`templates/gnc_planner_template.jl` template, and for rendezvous and proximity
operations the [RPO planner contract](reference/rpo_planner_contract.md).

## Add an actuator

Shape: hardware (H), plus declarations (D) and accessors (A). Rows:
[Actuators](generated/extension_hook_reference.md#Actuators).

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
- `CMGInnerGimbal`: the inner gimbal of a double-gimbal CMG. The CMG's gimbal
  is the outer one, fixed in the body, and the inner one turns the spin axis
  about an axis the outer gimbal carries round; its level is the inner gimbal
  rate, with its own rate limit and stops. Mount it on the root link after its
  CMG. The inner angles follow every outer angle in `cmg_angles`, the inner
  gimbals follow the CMGs in the mounting order, `cmg_torque_axes` gives a
  double-gimbal CMG's two torque axes and `cmg_steering` solves for every
  gimbal.

The hinges' and gimbals' angles and the firing state of thrusters with firing
dynamics are run state the control ticks change. Like a battery's flags and
the actuators' faults, they start each run as configured, so two runs of one
configuration with `isolate_state=false` start equal. A value changed between
two runs, after the last control tick, is kept, and later runs start from it.
A run resumed from a checkpoint continues the checkpoint's values.

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
checkpointed run keeps the faults' state from one segment to the next, and a
run resumed from its checkpoint file takes it from the file, a stuck
actuator's level at onset included.

Registration: `add_actuator!(spacecraft, link, actuator)`.

A control effector drives the mounted actuators by declaring
`SpaceAGORA.actuation_mode(c) = :mounted` and returning one level per mounted
actuator, in mounting order, from `SpaceAGORA.actuator_levels(c, spacecraft_index)`; it
computes the levels in `update_control!` and writes no force or mass-flow
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

Shape: hardware (H), plus declarations (D). Rows:
[Sensors](generated/extension_hook_reference.md#Sensors).

A sensor measures the state of the vehicle for navigation; it does not act on
the dynamics. A sensor kind is a mutable struct subtyping `AbstractSensor`
with its parameters and five fields the engine fills or reads: `period_s`
(the sampling period, finite and positive), `seed` (selects its random
stream), `rng`, `latest` (the measurement's components as a vector, or a
number for a one-component sensor) and `latest_time_s`. It implements:

- `sensor_measure(sensor, truth::SensorTruth, ctx::ActuatorContext, rng)`:
  one measurement of `truth` (the inertial position and velocity, the
  attitude, the body rate, the mass, the time and, when the sensor asks for
  them, the magnetic field, the Sun's direction and lit fraction with the
  planet's radius, the specific force with the body's angular acceleration,
  all in the body frame, and the planet-fixed position, velocity and
  altitude), with any noise drawn from `rng`.
  `ctx.link_to_body` maps the sensor's frame, its link's, to the body frame.
- `sensor_dimension(sensor)`: the number of components.
- optionally `sensor_requirements(sensor)`: `ActuatorRequirements(magnetic_field=true)`
  to receive the field (from the run's `MagneticTorqueRodModel`, or the
  tilted dipole) and `ActuatorRequirements(solar=true)` to receive the Sun's
  direction and the lit fraction of its disk (the run's ephemerides and the
  eclipse model), and `ActuatorRequirements(planet_frame=true)` to receive
  the planet-fixed position, velocity and altitude; the atmosphere is not
  available to sensors.
- optionally `sensor_needs_specific_force(sensor) = true` to receive the
  specific force: the vehicle's acceleration at the sample, from the solver,
  less the acceleration of the run's gravitational effectors (those for which
  `is_gravitational` is true; a user-defined gravity model should add a
  method returning true), and the body's angular acceleration, from the
  solver as well (an accelerometer away from the centre of mass adds the
  rotation's terms).
- optionally `sensor_period_s(sensor)`, the sampling period (by default its
  `period_s` field), and `initialize_sensor!(sensor, truth, ctx, rng)`, which
  resets at the start of each run any state the sensor keeps beyond the
  engine's fields, such as a filter. It runs once per run, just before the
  first measurement, with that measurement's `truth` and `ctx` and the fresh
  `rng`.

Mount it with `add_sensor!(spacecraft, link, sensor)`; a sensor on a hinged
panel turns with it. The built-in kinds are `RateGyro` (body rate with white
noise, a bias and a bias random walk), `StarTracker` (attitude with cross- and
about-boresight noise, boresight along its link's +z), `Magnetometer`
(field with white noise and a bias), `GNSSReceiver` (inertial position and
velocity with white noise), `SunSensor` (the Sun's direction with angular
noise, boresight along its link's +z, no reading outside its field of view or
in eclipse) and `Accelerometer` (specific force with white noise and a bias).

The engine samples every sensor at the run's start and every period after
it, from the accepted state; the samples are solver stops, and a checkpointed
run samples at the same times as a run in one piece. Each sensor's random
stream comes from its `seed`, the spacecraft's ID and its mounting position,
so a run repeats with the same seeds. Seeds are not varied for you: the
built-in kinds default to `seed=1`, so in a Monte Carlo campaign set each
sensor's `seed` from the sample's seed where the sample's spacecraft is
built, or the samples share one noise stream. A navigation effector reads the latest measurement in
`update_navigation!` with `sensor_measurement(p.args.dynamics_model.spacecraft[spacecraft_index], k)`
(and its time with `sensor_measurement_time`); sensors sample before
navigation effectors at the same time. Runs with sensors write the
`sc{i}_sensor_{k}` columns (see [Output Columns](user/output_columns.md#Sensors)).
A run resumed from a checkpoint file takes the sensors' streams, latest
measurements and gyro biases from the file and samples on the same grid; a
kind's own state beyond those fields is in the file when the kind implements
`checkpoint_state` (see [Run state in checkpoint files](#Run-state-in-checkpoint-files)). The
gravity-backbone and multirate solver modes refuse mounted sensors.

## Run state in checkpoint files

Shape: configuration (L). Rows:
[GNC lifecycle](generated/extension_hook_reference.md#GNC-lifecycle).

A checkpoint file holds the state vector and the run state the engine keeps
(see [Checkpointing](user/outputs.md#Checkpointing)). A guidance, navigation or
control model, or a sensor kind, that keeps run state of its own (a held
command, a plan, a filter's estimate, a counter) adds it with two methods:

```julia
SpaceAGORA.checkpoint_state(model) -> data      # plain data, or nothing
SpaceAGORA.restore_checkpoint_state!(model, data)
```

The engine calls `checkpoint_state` at every checkpoint and stores a copy of
the data. A run resumed from the file first calls the model's run-start
initializer (`initialize_navigation!`, `initialize_guidance!` or
`initialize_control!`), then `restore_checkpoint_state!` with the saved data,
which replaces the reset; a sensor's state is restored before its sampling
resumes, and a resumed run does not call `initialize_sensor!`. Return plain
data (numbers, arrays, named tuples, random streams) rather than references to
other objects of the run, and return `nothing` for the whole run or data for
the whole run.

```julia
mutable struct TickCounter <: SpaceAGORA.AbstractNavigationModel
    ticks::Int
end
SpaceAGORA.update_navigation!(m::TickCounter, u, p, t, spacecraft_index) = (m.ticks += 1; nothing)
SpaceAGORA.initialize_navigation!(m::TickCounter, u, p, t) = (m.ticks = 0; nothing)
SpaceAGORA.checkpoint_state(m::TickCounter) = (ticks = m.ticks,)
SpaceAGORA.restore_checkpoint_state!(m::TickCounter, data) = (m.ticks = data.ticks; nothing)
```

A checkpointed run refuses at setup a model whose `checkpoint_state` returns
data and which has no `restore_checkpoint_state!` method for it. A resumed run
refuses a checkpoint that holds no run state for a model that keeps some. When
the checkpoint holds model state, it also refuses one whose guidance, navigation
or control models or sensors differ from its own in type name or position;
without model state, those models' types are not compared. Every checkpoint also
records the spacecraft each guidance, navigation and control model acts on, after
the run-start initializers, and a resumed run refuses models acting on other
spacecraft before any `restore_checkpoint_state!`, and a checkpoint written by
SpaceAGORA 0.5.1 or earlier, which holds no such record (see
[Checkpointing](user/outputs.md#Checkpointing)); a model that acts on some
spacecraft only should say so with `spacecraft_scope` for that check to cover it.
Mounted sensors are checked
separately as well, with the rest of the run's callback state: each sensor's type
name, mounting link, period and seed must match, with or without model state.
A type name is the short name, without the module or type parameters, so a
same-named type from another module, or another parametrization of the same
type, is not told apart.
Without the two methods,
a resumed run starts a guidance, navigation or control model's own state from
its reset, and a sensor's from what the sensor object holds. The conformance
check `checkpoint_round_trip` (added in 0.4.0) compares a run resumed from a
checkpoint with an uninterrupted one, so it finds run state the two methods
leave out; it is required from 0.5.0. The control and actuator templates
(`templates/control_hook_template.jl`, `templates/actuator_template.jl`)
implement both methods for the command they hold.

The built-in models that keep run state implement them: the momentum manager,
the Apollo descent guidance and control, the periapsis-raise guidance, the burn
scheduler (`BaseThrusterModel`), the RPO planner binding with the plan buffer it
writes, the RPO MPC controller (its held command, warm start and solver iterate)
and the robot-arm controller. The
energy-depletion guidance and the GNC algorithm adapter (`AlgorithmGuidance`)
refuse checkpointed runs.

## Planet, ephemerides, thermal and thruster extensions

Stable interfaces:

- `SpaceAGORA.AbstractPlanet`
- `SpaceAGORA.AbstractEphemeridesModel`
- `SpaceAGORA.AbstractThermalModel` (its hook is `getHeatRate`, in W/cm²)
- `SpaceAGORA.AbstractThrusterModel`

These interfaces are public because no-GRAM onboarding and future user-defined
body/back-end, thermal and thruster work need a stable contract, but their
extension recipes are still thinner than the families above. An ephemerides
model implements the model-first hooks listed under
[Ephemerides models](user/simulation_configuration.md#Ephemerides-models)
(`ephemerides_frame_lpi(model, planet, et)` and the rest); thermal models are
consumed inside the heat-load machinery, and thruster models inside control
effectors. Treat them as advanced extension points: subtype the abstract type
and follow an existing implementation in `src/models/environment/`,
`src/models/vehicle/thermal/` or `src/models/vehicle/actuators/` rather than a
documented recipe.

## Design guidance

When adding a new model:

1. keep state in the model object, not in global variables
2. make units explicit in field names or docstrings
3. keep expensive caches out of hot dispatch paths unless they are typed
4. add a small typed constructor and one smoke test
5. document the model through the root `SpaceAGORA` hook surface, not by telling
   users to import internal modules

## Templates and runnable examples

The five files under `templates/` define different extension contracts: a force
or torque, an atmosphere density model, a controller in the feedback profile
(`templates/control_hook_template.jl`), an actuator, and a planner in the planning
profile (`templates/gnc_planner_template.jl`). They define model types and
methods; they do not construct or run a mission. Copy the relevant file into
your model's owning package and add tests for the behavior you implement, with
the [conformance checks](user/gnc_algorithm.md#Test-your-algorithm). The
force/torque template implements `wrench`;
[Adding a Force or Torque of Your Own](user/custom_effector.md) walks through
one in detail, and [Writing a GNC Algorithm](user/gnc_algorithm.md) through a
controller, an estimator and a planner. For a repository of your own, the
starters under `templates/starters/` add tests and CI
([Building Your Own Research Repository](user/research_repository.md)).

Use the [Examples Catalog](user/examples_catalog.md) when you want a complete
scenario to run or adapt. Helpers under `examples/support/` serve those scenarios
and are not extension contracts.
