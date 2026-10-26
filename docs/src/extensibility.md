# Extensibility

New models and controllers plug into SpaceAGORA through hook functions of the
root package. Every kind of extension follows the same pattern:

1. subtype an abstract type that `SpaceAGORA` exports;
2. add a method of the matching root hook for your type;
3. put an instance of it in the run's typed configuration.

Build on these names rather than on SpaceAGORA's internal modules, unless you are
changing SpaceAGORA itself. The quickest start is a copy of
`templates/force_torque_model_template.jl`. The walkthroughs go further:

- [Adding a Force or Torque of Your Own](user/custom_effector.md), a tested walkthrough of the force/torque hook;
- [Adding an Actuator of Your Own](user/custom_actuator.md), the same for the actuator contract;
- [Writing a GNC Algorithm](user/gnc_algorithm.md), on controllers, estimators and planners, with worked examples that CI runs;
- [Building Your Own Research Repository](user/research_repository.md), on starters, conformance checks and CI for work kept in a repository of your own;
- the [Public API](generated/public_api.md) and [Concepts](user/concepts.md) pages.

Every hook follows the [Extension Conventions](reference/extension_conventions.md):
one argument order (time before the spacecraft index), one name for each kind of
state, and stated frames, units, signs and names. The
[Extension Hook Reference](generated/extension_hook_reference.md) lists each hook
with its shape, signature and fallback, and
[Migrating Extension Hooks](reference/migrating_extension_hooks.md) maps the old
names to their replacements. Each section below gives the shape of its hooks and
links their rows in the reference.

An optional package with a wider integration reaches SpaceAGORA through a named,
versioned facade instead. The [Optional Adapter Standard](reference/adapter_contract.md)
sets the rules on ownership, compatibility, installation and conformance;
[External Propagation Services](reference/external_propagation_services.md) is
the facade for external spacecraft integrators, and HYPR uses its
[versioned planning services](reference/hypr_services.md).

## Add a force or torque model

Shape: right-hand side (R), plus declarations (D). Rows:
[Physical effectors](generated/extension_hook_reference.md#Physical-effectors).

Subtype `SpaceAGORA.AbstractForceTorqueModel` and implement `wrench`. The other
hooks are optional and have defaults. Their signatures:

```julia
SpaceAGORA.wrench(model, x::SpaceAGORA.StateSample, env::SpaceAGORA.EnvironmentSample, t::Float64) -> (force_ii, torque_body)
SpaceAGORA.environment_requirements(model) -> SpaceAGORA.EffectorEnvironmentRequirements
SpaceAGORA.solver_partition(model) -> :explicit | :implicit
SpaceAGORA.gravity_backbone_structure(model) -> :unsupported | :position_only_static_gravity
SpaceAGORA.gravity_backbone_acceleration_ii(model, x::SpaceAGORA.StateSample, env::SpaceAGORA.EnvironmentSample, t::Float64) -> accel_ii
SpaceAGORA.gravity_backbone_kick_structure(model) -> :unsupported | :velocity_kick_explicit
SpaceAGORA.gravity_backbone_kick_acceleration_ii(model, x::SpaceAGORA.StateSample, env::SpaceAGORA.EnvironmentSample, t::Float64) -> accel_ii
```

The engine samples the state and the environment consistently for each solver
stage and passes them in as `x` and `env`; `x.spacecraft_index` is the
spacecraft's position in the run. `wrench` returns the force in the inertial
frame and the torque in the body frame. Type the sample `x::SpaceAGORA.StateSample`,
or leave it untyped
([the supported sample argument](reference/extension_conventions.md#The-supported-sample-argument)).
`environment_requirements` declares what the engine samples into `env`; by
default it samples nothing. To have an effector act on some spacecraft only,
give its entry as `scoped(effector, scope)`
([Effectors on part of a constellation](user/constellations.md#Effectors-on-part-of-a-constellation)).

The other declarations matter only in particular solver modes:

- `solver_partition` places the effector under `split_imex`, the
  atmosphere-implicit IMEX mode: `:implicit`, or `:explicit` (the default). The
  other split modes do not read it. `multirate` is the control-focused split,
  and `multirate_mri` uses a fixed partition in which custom dynamic effectors
  are always on the slow side.
- The gravity-backbone hooks let the effector take part in
  `gravity_backbone_split`, a fixed-step mode with a symplectic gravity core and
  explicit velocity kicks for terms such as SRP and N-body gravity. An effector
  joins as a position-only static-gravity core term or as an explicit
  translational kick; one that declares neither is refused in this mode. The
  mode is translational only, not a symplectic solve of the whole system.

The heat load in the state is the integral of the heat rate and always stays on
the explicit side.

Registration: add the model to the `dynamic_effectors` of the `DynamicsModel`.

Template: `templates/force_torque_model_template.jl`.

## Add a density model

Shape: query (Q). Rows:
[Environment queries](generated/extension_hook_reference.md#Environment-queries).

Subtype `SpaceAGORA.AbstractDensityModel` and implement `getDensity`:

```julia
SpaceAGORA.getDensity(model, h, lat, lon, el_time, wind[, p]) -> (rho, temperature, wind_vec)
```

`h` is geodetic altitude in metres, `lat` geodetic latitude and `lon` east
longitude in radians, and `el_time` the elapsed seconds since the scenario epoch.
`wind` is a `Bool` saying whether the caller needs the wind vector: the force and
heating paths pass the run's wind setting (`EnvironmentModel(wind=...)`), and
density-only callers such as the visualization pass `false`. A model may return a
zero `wind_vec` when `wind` is `false`. The results are density in kg/m³,
temperature in K and `wind_vec`, the air's velocity relative to the rotating
planet as east, north and up components in m/s.

A run calls the seven-argument form, with the run's parameters `p` last, and
refuses before integration a density model that has no method of it. A model that
needs the calendar date reads it from `p.args.initial_time`. The template defines
the six-argument form and a seven-argument form that forwards to it.

The batch form is optional. Its default fills each output from the scalar form,
query by query, and `el_time` is one elapsed time for every query or one per
query:

```julia
SpaceAGORA.getDensityBatch!(rhos, Ts, winds, model, hs, lats, lons, el_time, wind, p)
```

Registration: set the model as `environment_model.density_model`.

Template: `templates/density_model_template.jl`.

## Add a control effector or hook

Shape: tick (T) for `update_control!`, right-hand side (R) for the outputs.
Rows: [GNC ticks](generated/extension_hook_reference.md#GNC-ticks) and
[Control outputs](generated/extension_hook_reference.md#Control-outputs).

Subtype `SpaceAGORA.AbstractControlEffectorModel`. The tick and the output hooks
have these signatures:

```julia
SpaceAGORA.update_control!(model, u, p, t, spacecraft_index)
SpaceAGORA.control_wrench(model, x::SpaceAGORA.StateSample, p, t::Float64) -> (force_ii, torque_body)
SpaceAGORA.control_propellant_rate(model, x::SpaceAGORA.StateSample, p, t::Float64) -> kg/s, positive
SpaceAGORA.control_wheel_torque(model, x::SpaceAGORA.StateSample, p, t::Float64) -> nothing or torque_body
```

`update_control!` is the control tick: it computes the command and holds it. The
three output hooks return the held command at every right-hand-side evaluation,
for each spacecraft in the model's scope, and `x.spacecraft_index` is the
spacecraft's position in the run. `control_propellant_rate` returns the
propellant used, as a positive number. An effector that uses no propellant, or
drives no wheel, needs no method of `control_propellant_rate` or
`control_wheel_torque`: their defaults return `0.0` and `nothing`.

A controller also declares what it acts on and how it runs (`spacecraft_scope`,
`gnc_period_s`, `skip_inactive`, `held_command`), can share an object with an
estimator or a planner through `gnc_links`, and resets its run state at the start
of each run in `initialize_control!`.
[Writing a GNC Algorithm](user/gnc_algorithm.md#The-feedback-profile) explains
these declarations with a worked controller that CI runs.

Registration: add the model to `ControlModel(control_effectors=(...), control_rates=[...])`,
with one control period in seconds for each effector.

Template: `templates/control_hook_template.jl`.

### Per-spacecraft declaration

A `SpacecraftModel` can carry its own GNC: `SpacecraftModel(...; control=ControlModel(...))`,
and likewise `guidance` and `navigation`. When `run_simulation` starts, it passes
each such effector through `SpaceAGORA.bind_spacecraft(effector, spacecraft_index)`
with the position of its spacecraft, and appends the result to the
configuration-level tuples. An effector that selects its vehicle by an index field
implements `bind_spacecraft` by returning a copy with that field set. It also
declares the spacecraft it acts on with `spacecraft_scope`, so that the engine
ticks it and calls its hooks for that spacecraft only. The fragment below shows
both methods; it assumes `using SpaceAGORA, StaticArrays`, and the effector's
hooks are left out:

```julia
struct MyPush <: SpaceAGORA.AbstractControlEffectorModel
    spacecraft_index::Int
    force_n::SVector{3, Float64}
end

SpaceAGORA.bind_spacecraft(m::MyPush, spacecraft_index::Int) = MyPush(spacecraft_index, m.force_n)
SpaceAGORA.spacecraft_scope(m::MyPush) = BoundSpacecraft((m.spacecraft_index,))
```

A bound copy whose `spacecraft_scope` is still `AllSpacecraft()` acts on every
spacecraft, as it would at configuration level, and each run warns once for each
such type, naming the run position and id of the spacecraft. The default
`bind_spacecraft` throws an `ArgumentError`: an effector without an index acts on
every spacecraft, and one that couples several vehicles (the RPO chaser and
target types) cannot belong to one, so both stay at configuration level. So does
the burn scheduler, `BaseThrusterModel`, which is one model with one entry per
spacecraft. The built-in models that bind are the momentum manager, the
robot-arm control effector, the Apollo descent guidance and control models, the
algorithm guidance binding (`AlgorithmGuidance`) and, from SpaceAGORA 0.7.0, the
two maneuver guidance models, `ApoapsisTargetPeriapsisRaiseGuidanceModel` and
`AerobrakingCampaignPropulsiveManeuverGuidanceModel`.

## Add a guidance or navigation effector

Shape: tick (T). Rows:
[GNC ticks](generated/extension_hook_reference.md#GNC-ticks).

Subtype `SpaceAGORA.AbstractGuidanceModel` and implement `SpaceAGORA.update_guidance!`,
or subtype `SpaceAGORA.AbstractNavigationModel` and implement
`SpaceAGORA.update_navigation!`. Both ticks take the arguments of
`update_control!`:

```julia
SpaceAGORA.update_guidance!(model, u, p, t, spacecraft_index)
SpaceAGORA.update_navigation!(model, u, p, t, spacecraft_index)
```

A navigation effector takes its measurements or updates its estimate; a guidance
effector updates a plan, a target or a mode. Each keeps the result in its own
state for the control effectors to read. When ticks coincide, navigation runs
first, then control, then guidance. Two kinds of guidance run before control
instead: energy-depletion guidance with targeting, and the algorithm binding
`AlgorithmGuidance`, so that the controller sees a product at the tick it is
installed.

Registration:

- add guidance effectors to `GuidanceModel(guidance_effectors=(...), guidance_rates=[...])`;
- add navigation effectors to `NavigationModel(navigation_effectors=(...), navigation_rates=[...])`.

[Writing a GNC Algorithm](user/gnc_algorithm.md#Worked-example-B:-an-estimator-and-a-link)
has a worked estimator that a controller reads through a declared link.

A planner, whose output is used later or can be checked before use, follows the
planning profile instead: it proposes, and its family validates and installs the
proposal. See [Writing a GNC Algorithm](user/gnc_algorithm.md#The-planning-profile),
the template `templates/gnc_planner_template.jl` and, for rendezvous and
proximity operations, the [RPO planner contract](reference/rpo_planner_contract.md).

## Add an actuator

Shape: hardware (H), plus declarations (D) and accessors (A). Rows:
[Actuators](generated/extension_hook_reference.md#Actuators).

An actuator is hardware mounted on a spacecraft link that acts on the vehicle
when commanded. You describe the hardware once. The engine derives the force, the
torque, the propellant use and the result columns (each actuator's applied level,
the actuators' power draw and the thrusters' firing levels) from the actuators
mounted on the vehicle.

The hooks and helpers:

- `SpaceAGORA.AbstractActuator`, the type to subtype;
- `SpaceAGORA.actuator_wrench` and `SpaceAGORA.actuator_limits`, required;
- `SpaceAGORA.actuator_consumption`, `SpaceAGORA.actuator_saturate`,
  `SpaceAGORA.actuator_fraction`, `SpaceAGORA.actuator_unit_wrench` and
  `SpaceAGORA.actuator_requirements`, optional, with defaults;
- `SpaceAGORA.add_actuator!`, `SpaceAGORA.mounted_actuators`,
  `SpaceAGORA.allocation_matrix` and `SpaceAGORA.allocate_actuators`, for mounting
  and allocation.

Typical methods, for an actuator type `MyActuator`:

```julia
SpaceAGORA.actuator_limits(a::MyActuator) = (lower, upper)          # level bounds, the kind's unit
SpaceAGORA.actuator_wrench(a::MyActuator, level, ctx) = (force_body, torque_body)
SpaceAGORA.actuator_consumption(a::MyActuator, level, ctx) = (mass_rate_kg_s = ..., power_w = ...)
```

The wrench is in the spacecraft body frame, with the torque about the body
reference point. `ctx`, an `ActuatorContext`, gives the mounting link's pose in
the body frame and the vehicle's mass, attitude and rate.

A kind that needs more of its surroundings declares it with
`SpaceAGORA.actuator_requirements(a) = ActuatorRequirements(...)`: the
planet-fixed frame, the atmosphere, the Sun's position or the local magnetic
field. The engine samples what the mounted actuators ask for once per evaluation
and passes it in the context: the field as `ctx.magnetic_field_body`, the rest
as `ctx.environment`, an `EnvironmentSample` (for example
`ctx.environment.atmosphere.rho_kg_m3`). Without a request, `ctx.environment` is
`nothing`.

The built-in kinds:

- `Thruster`: its level is a thrust in newtons, and its `direction` is the
  direction of the thrust on the vehicle. With `firing_dynamics=true` the command
  goes through the thruster's minimum on-time, Schmitt trigger and thrust ramp at
  each control tick, and the vehicle feels the resulting average thrust until the
  next one.
- `ReactionWheel`: its level is the torque on the body about its spin axis. It
  mounts on the root link, and its momentum is integrated, with torque and
  momentum limits. Given a spin inertia, a base power and a motor efficiency, it
  reports the power it draws.
- `Magnetorquer`: its level is the commanded dipole along its axis, and it acts in
  the run's magnetic field. Given its power at full dipole, it reports a power
  that grows with the square of the dipole.
- `LaserTerminal`: its level is the transmitted optical power. The vehicle recoils
  with the power over the speed of light, and the terminal draws the power over
  its efficiency.
- `PanelHinge`: its level is the angle of the non-root link it turns, with angle
  and rate limits. The engine applies it at the commanding effector's control
  ticks. The energy-depletion controller also moves its panels through a hinge
  mounted on them.
- `ThrusterGimbal`: one axis of a gimbal on a thruster mounted on the same link.
  Its level is the gimbal angle, with angle and rate limits, applied at control
  ticks as a hinge's is. The thrust direction turns about the gimbal axis (the
  nozzle stays where it is), so the thruster's force, torque, propellant use and
  allocation column follow it. Two gimbals on one thruster make a two-axis
  gimbal, turning in mounting order. `thrust_direction(thruster)` gives the turned
  direction.
- `ControlMomentGyro`: a single-gimbal CMG on the root link, a rotor of constant
  momentum whose spin axis turns about the gimbal axis. Its level is the gimbal
  rate, within its rate limit, and its gimbal stops at
  `min_gimbal_angle_rad`/`max_gimbal_angle_rad` (none by default). The gimbal
  angle is integrated with the vehicle's state (`cmg_angles`, present only when a
  CMG is mounted); the vehicle feels `-δ̇ (g × h)` and the rotor's momentum in its
  gyroscopic term. `cmg_steering(spacecraft, angles, torque)` turns a desired
  torque into gimbal rates with a singularity-robust inverse. The gravity-backbone
  mode and articulated spacecraft refuse CMGs.
- `CMGInnerGimbal`: the inner gimbal of a double-gimbal CMG. The CMG's own gimbal
  is the outer one, fixed in the body; the inner one turns the spin axis about an
  axis that the outer gimbal carries round. Its level is the inner gimbal rate,
  with its own rate limit and stops. Mount it on the root link after its CMG. In
  `cmg_angles` the inner angles follow all the outer angles, in mounting order.
  `cmg_torque_axes` gives a double-gimbal CMG's two torque axes, and
  `cmg_steering` solves for every gimbal.

The angles of hinges and gimbals, and the firing state of thrusters with firing
dynamics, are run state that the control ticks change. Like a battery's flags and
the actuators' faults, they start each run as configured, so two runs of one
configuration with `isolate_state=false` start equal. A value changed between two
runs, after the last control tick, is kept, and later runs start from it. A run
resumed from a checkpoint continues the checkpoint's values.

The actuators' power feeds the spacecraft's power budget when it carries one:
solar arrays, a battery and constant loads mounted with `add_power!`. See
[Powering a Spacecraft](user/power_budget.md).

A mounted actuator can fail. `add_fault!(spacecraft, ActuatorFault(actuator;
kind, start_s, stop_s, value))` makes it, between `start_s` and `stop_s`, apply
level zero (`:failed`; a failed hinge or gimbal holds its angle), hold a level
(`:stuck`, by default the level commanded when the fault began), act at a
fraction of its level and limits (`:degraded`) or act with an offset (`:biased`).
The fault switches exactly at those times. It acts on what the vehicle feels, on
the propellant, the power and the `actuator_level` columns; the controller keeps
seeing its own commands. On a thruster with firing dynamics the fault acts on the
thrust it delivers, so it takes effect at its start, between control ticks. Runs
with faults write `sc{i}_fault_active_{j}`. A checkpointed run keeps the faults'
state from one segment to the next, and a run resumed from its checkpoint file
takes it from the file, including a stuck actuator's level at onset.

Registration: `add_actuator!(spacecraft, link, actuator)`.

A control effector drives the mounted actuators by declaring
`SpaceAGORA.actuation_mode(c) = :mounted` and returning one level per mounted
actuator, in mounting order, from `SpaceAGORA.actuator_levels(c, spacecraft_index)`.
It computes the levels in `update_control!` and needs no `control_wrench` or
propellant-rate method. Effectors in the default `:ideal` mode keep applying their
own hooks. The built-in RPO and Apollo descent controllers and the magnetic
momentum manager accept `actuation=:mounted`. So does the orbit-maneuver burn
scheduler `BaseThrusterModel`: it fires a mounted engine during the burn window
and points it along the burn with the mounted reaction wheels, and the engine
switches exactly at the window's ends, between integration steps. A run with a
mounted effector writes the `sc{i}_actuator_level_{k}` and `sc{i}_actuator_power`
columns (see [Output Columns](user/output_columns.md#Mounted-actuators));
`applied_actuator_levels` and `mounted_power` compute the same values from a
spacecraft and its levels.

Template: `templates/actuator_template.jl`.

Walkthrough: [Adding an Actuator of Your Own](user/custom_actuator.md).
Example: `examples/Earth_CubeSat_Mounted_Actuators.jl`, a CubeSat whose mounted
wheels detumble it while its torque rods unload them.

## Add a sensor

Shape: hardware (H), plus declarations (D). Rows:
[Sensors](generated/extension_hook_reference.md#Sensors).

A sensor measures the state of the vehicle for navigation; it does not act on the
dynamics. A sensor kind is a mutable struct subtyping `AbstractSensor`. Besides
its parameters, it has five fields that the engine fills or reads: `period_s`
(the sampling period, finite and positive), `seed` (which selects its random
stream), `rng`, `latest` (the measurement's components as a vector, or a number
for a one-component sensor) and `latest_time_s`. It implements:

- `sensor_measure(sensor, truth::SensorTruth, ctx::ActuatorContext, rng)`: one
  measurement of `truth`, with any noise drawn from `rng`. `truth` holds the
  inertial position and velocity, the attitude, the body rate, the mass and the
  time and, when the sensor asks for them, the magnetic field, the Sun's
  direction and lit fraction with the planet's radius, the specific force with
  the body's angular acceleration (all in the body frame), and the planet-fixed
  position, velocity and altitude. `ctx.link_to_body` maps the sensor's frame,
  which is its link's, to the body frame.
- `sensor_dimension(sensor)`: the number of components.
- optionally `sensor_requirements(sensor)`, to receive more of the truth:
  `ActuatorRequirements(magnetic_field=true)` for the field (from the run's
  `MagneticTorqueRodModel`, or the tilted dipole), `ActuatorRequirements(solar=true)`
  for the Sun's direction and the lit fraction of its disk (from the run's
  ephemerides and the eclipse model), and `ActuatorRequirements(planet_frame=true)`
  for the planet-fixed position, velocity and altitude. The atmosphere is not
  available to sensors.
- optionally `sensor_needs_specific_force(sensor) = true`, to receive the
  specific force and the body's angular acceleration. The specific force is the
  vehicle's acceleration at the sample, from the solver, less the acceleration of
  the run's gravitational effectors, those for which `is_gravitational` is true
  (a user-defined gravity model should add a method that returns true). The
  angular acceleration also comes from the solver; an accelerometer away from the
  centre of mass adds the rotation's terms itself.
- optionally `sensor_period_s(sensor)`, the sampling period (by default the
  `period_s` field), and `initialize_sensor!(sensor, truth, ctx, rng)`, which
  resets any state the sensor keeps beyond the engine's fields, such as a filter,
  at the start of each run. It runs once per run, immediately before the first
  measurement, with that measurement's `truth` and `ctx` and the fresh `rng`.

Mount a sensor with `add_sensor!(spacecraft, link, sensor)`; a sensor on a hinged
panel turns with it. The built-in kinds are `RateGyro` (body rate with white
noise, a bias and a bias random walk), `StarTracker` (attitude with cross- and
about-boresight noise, boresight along its link's +z), `Magnetometer` (field with
white noise and a bias), `GNSSReceiver` (inertial or planet-fixed position and
velocity with white noise), `SunSensor` (the Sun's direction with angular noise,
boresight along its link's +z, no reading outside its field of view or in
eclipse) and `Accelerometer` (specific force with white noise and a bias).

The engine samples every sensor at the run's start and every period after it,
from the accepted state. The sample times are solver stops, and a checkpointed run
samples at the same times as a run in one piece. Each sensor's random stream comes
from its `seed`, the spacecraft's ID and its mounting position, so a run repeats
with the same seeds. Seeds are not varied for you, and the built-in kinds default
to `seed=1`. In a Monte Carlo campaign, set each sensor's `seed` from the sample's
seed where the sample's spacecraft is built; otherwise the samples share one noise
stream.

A navigation effector reads the latest measurement in `update_navigation!` with
`sensor_measurement(p.args.dynamics_model.spacecraft[spacecraft_index], k)`, and
its time with `sensor_measurement_time`. Sensors sample before the navigation
effectors due at the same time. Runs with sensors write the `sc{i}_sensor_{k}`
columns (see [Output Columns](user/output_columns.md#Sensors)). A run resumed from
a checkpoint file takes the sensors' streams, latest measurements and gyro biases
from the file and samples on the same grid. A kind's own state beyond those
fields is in the file when the kind implements `checkpoint_state` (see
[Run state in checkpoint files](#Run-state-in-checkpoint-files)). The
`:gravity_backbone_split` and `:multirate` solver modes refuse mounted sensors.

## Run state in checkpoint files

Shape: configuration (L). Rows:
[GNC lifecycle](generated/extension_hook_reference.md#GNC-lifecycle).

A checkpoint file holds the state vector and the run state the engine keeps (see
[Checkpointing](user/outputs.md#Checkpointing)). A guidance, navigation or
control model, or a sensor kind, that keeps run state of its own (a held command,
a plan, a filter's estimate, a counter) adds it with two methods:

```julia
SpaceAGORA.checkpoint_state(model) -> data      # plain data, or nothing
SpaceAGORA.restore_checkpoint_state!(model, data)
```

The engine calls `checkpoint_state` at every checkpoint and stores a copy of the
data. A run resumed from the file first calls the model's run-start initializer
(`initialize_navigation!`, `initialize_guidance!` or `initialize_control!`), then
`restore_checkpoint_state!` with the saved data, which replaces the reset. A
sensor's state is restored before its sampling resumes, and a resumed run does not
call `initialize_sensor!`. Return plain data (numbers, arrays, named tuples,
random streams) rather than references to other objects of the run. A model
returns `nothing` for the whole run, or data for the whole run, not one at some
checkpoints and the other at others.

For example, a navigation model that counts its ticks:

```julia
mutable struct TickCounter <: SpaceAGORA.AbstractNavigationModel
    ticks::Int
end
SpaceAGORA.update_navigation!(m::TickCounter, u, p, t, spacecraft_index) = (m.ticks += 1; nothing)
SpaceAGORA.initialize_navigation!(m::TickCounter, u, p, t) = (m.ticks = 0; nothing)
SpaceAGORA.checkpoint_state(m::TickCounter) = (ticks = m.ticks,)
SpaceAGORA.restore_checkpoint_state!(m::TickCounter, data) = (m.ticks = data.ticks; nothing)
```

The engine checks the models against the file:

- A checkpointed run refuses at setup a model whose `checkpoint_state` returns
  data but which has no `restore_checkpoint_state!` method for it.
- A resumed run refuses a checkpoint that holds no run state for a model that
  keeps some.
- When the checkpoint holds model state, a resumed run refuses it if the run's
  guidance, navigation or control models or sensors differ from the checkpoint's
  in type name or position. Without model state, those models' types are not
  compared.
- Every checkpoint records the spacecraft that each guidance, navigation and
  control model acts on, taken after the run-start initializers. Before any
  `restore_checkpoint_state!`, a resumed run refuses models that act on other
  spacecraft than the record says, and a run with such models refuses a
  checkpoint written by SpaceAGORA 0.5.1 or earlier, which has no record (see
  [Checkpointing](user/outputs.md#Checkpointing)). A model that acts on some
  spacecraft only should say so with `spacecraft_scope`, so that this check
  covers it.
- Mounted sensors are checked separately, with the rest of the run's callback
  state, with or without model state: each sensor's type name, mounting link,
  period and seed must match.

A type name here is the short name, without the module or type parameters, so a
same-named type from another module, or another parametrization of the same
type, is not told apart.

Without the two methods, a resumed run starts a guidance, navigation or control
model's own state from its reset, and a sensor's from what the sensor object
holds. The conformance check `checkpoint_round_trip` compares a run resumed from a
checkpoint with an uninterrupted one, so it finds run state that the two methods
leave out. It was added in 0.4.0 and has been required since 0.5.0. The control
and actuator templates (`templates/control_hook_template.jl`,
`templates/actuator_template.jl`) implement both methods for the command they
hold.

The built-in models that keep run state implement them: the momentum manager,
the Apollo descent guidance and control, the periapsis-raise guidance, the burn
scheduler (`BaseThrusterModel`), the RPO planner binding with the plan buffer it
writes, the RPO MPC controller (its held command, warm start and solver iterate)
and the robot-arm controller. The energy-depletion guidance and the GNC algorithm
binding (`AlgorithmGuidance`) refuse checkpointed runs.

## Planet, ephemerides, thermal and thruster extensions

`SpaceAGORA.AbstractPlanet`, `SpaceAGORA.AbstractEphemeridesModel`,
`SpaceAGORA.AbstractThermalModel` (whose hook, `getHeatRate`, returns W/cm²) and
`SpaceAGORA.AbstractThrusterModel` are public, so that runs without GRAM and new
bodies, ephemeris back-ends, thermal and thruster models have a stable type to
build on. They have no recipe of the kind above yet. An ephemerides model
implements the model-first hooks listed under
[Ephemerides models](user/simulation_configuration.md#Ephemerides-models)
(`ephemerides_frame_lpi(model, planet, et)` and the rest). Thermal models are
called inside the heat-load computation, and thruster models inside control
effectors. To write one, subtype the abstract type and follow an existing
implementation in `src/models/environment/`, `src/models/vehicle/thermal/` or
`src/models/vehicle/actuators/`.

## Design guidance

When you add a model:

1. keep its state in the model object, not in global variables;
2. make units explicit in field names or docstrings;
3. keep expensive caches out of hot dispatch paths unless they are typed;
4. add a small typed constructor and one smoke test;
5. document the model through the root `SpaceAGORA` hooks, without telling users
   to import internal modules.

## Templates and runnable examples

The five files under `templates/` cover different extension contracts: a force or
torque, an atmosphere density model, a controller in the feedback profile
(`templates/control_hook_template.jl`), an actuator, and a planner in the planning
profile (`templates/gnc_planner_template.jl`). They define model types and
methods; they do not build or run a mission. Copy the one you need into the
package that will own your model, and add tests for the behaviour you implement,
with the [conformance checks](user/gnc_algorithm.md#Test-your-algorithm). The
force/torque template implements `wrench`;
[Adding a Force or Torque of Your Own](user/custom_effector.md) walks through one
in detail, and [Writing a GNC Algorithm](user/gnc_algorithm.md) through a
controller, an estimator and a planner. For a repository of your own, the starters
under `templates/starters/` add tests and CI
([Building Your Own Research Repository](user/research_repository.md)).

For a complete scenario to run or adapt, see the
[Examples Catalog](user/examples_catalog.md). The helpers under
`examples/support/` serve those scenarios and are not extension contracts.
