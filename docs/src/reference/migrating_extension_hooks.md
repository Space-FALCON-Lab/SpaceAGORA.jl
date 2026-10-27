# Migrating Extension Hooks

SpaceAGORA's extension hooks follow one set of [conventions](extension_conventions.md):
time before the spacecraft index, a state sample for right-hand-side hooks, and
positive consumption. Ten hook names, and one form of the plume model's `wrench`,
did not. SpaceAGORA 0.3.0 deprecated them and 0.5.0 removed them. In 0.3.0 and
0.4.0 the pre-run check warned once per process about each model type's method of
an old name, a direct call of an old name warned too, and the five-argument plume
form warned only when called directly. This page lists each name with its
replacement, and the later changes to how the hooks are called.

A method of a removed name no longer runs. One that names the package, such as
`SpaceAGORA.calcForceTorque(...) = ...`, fails when its file loads, because the
name is gone. The pre-run check refuses a model left without a method of a
required replacement (a tick, `control_wrench` unless the effector's actuation is
`:mounted`, or `wrench`) and names the hook it needs. `control_propellant_rate`
and `control_wheel_torque` are optional: without a method of its own, an effector
gets `0.0` and `nothing`, and the check does not refuse it. An ephemerides model
without the model-first hooks stops with a `MethodError` naming the hook when the
run first asks it for its epoch or frame. To find the old names in your code:

```bash
grep -rnE "calcForceTorque|calcControl(ForceTorque|MassFlowRate|Effect!)|calcReactionWheelTorque|calc(Guidance|Navigation)Effect!|planet_frame_lpi|ephemerides_(sun_direction_ii|time_seconds)" src test
```

## Removed hook names

| Removed in 0.5.0 | Replacement | What changes |
|---|---|---|
| `calcForceTorque(model, x, p, i)` | `wrench(model, x::StateSample, env::EnvironmentSample, t::Float64)` | a state sample instead of the raw view; the environment the engine sampled instead of `p`; the time as an argument and the index in `x.spacecraft_index` |
| `wrench(model::PlumeSurfaceInteractionModel, x, env, t, i)` | `wrench(model, x, env, t)` | the index in `x.spacecraft_index` instead of an argument |
| `calcControlForceTorque(model, u, p, i, t)` | `control_wrench(model, x::StateSample, p, t::Float64)` | a state sample instead of the raw view, the time before the index, the index in `x.spacecraft_index` |
| `calcControlMassFlowRate(model, u, p, i, t)` | `control_propellant_rate(model, x::StateSample, p, t::Float64)` | the same, and the sign (see below) |
| `calcReactionWheelTorque(model, u, p, i, t)` | `control_wheel_torque(model, x::StateSample, p, t::Float64)` | the same as for the wrench |
| `calcControlEffect!(model, u, p, t, i)` | `update_control!(model, u, p, t, spacecraft_index)` | the name only |
| `calcGuidanceEffect!(model, u, p, t, i)` | `update_guidance!(model, u, p, t, spacecraft_index)` | the name only |
| `calcNavigationEffect!(model, u, p, t, i)` | `update_navigation!(model, u, p, t, spacecraft_index)` | the name only |
| `planet_frame_lpi(planet, et, model)` | `ephemerides_frame_lpi(model, planet, et)` | the model first |
| `ephemerides_sun_direction_ii(planet, et, model)` | `ephemerides_sun_ii(model, planet, et)` | the model first, and the name |
| `ephemerides_time_seconds(initial_time, model)` | `ephemerides_epoch_seconds(model, initial_time)` | the model first, and the name |

!!! warning "The sign changes"
    `calcControlMassFlowRate` returned the spacecraft's mass derivative, which is
    negative while propellant is used. `control_propellant_rate` returns the
    propellant used, which is positive. Flip the sign when you rename the method:
    a method renamed without it makes the vehicle gain mass.

The frames, units and values are otherwise unchanged. An effector that uses no
propellant, or drives no wheel, needs no `control_propellant_rate` or
`control_wheel_torque`: their defaults return `0.0` and `nothing`. A force or
torque model declares what it reads from the environment with
`environment_requirements(model)`. [Adding a Force or Torque of Your Own](../user/custom_effector.md)
and [Extending SpaceAGORA](../extensibility.md) show complete methods of the
replacements.

## Other deprecations in the 0.3.0 cycle

| Item | Replacement | Status | Removed in |
|---|---|---|---|
| `RPOGuidanceModel` | `make_rpo_configuration` with a planner, for example `planner=HYPRRPOPlanner(RPOPSOConfig(...))`, which configures a `PlannerGuidance` | warned before the run from 0.3.0; 0.4.0 is the last release that has it | 0.5.0 |
| the index fields and keywords `MagneticMomentumManagerModel.sat_idx`, `RobotArmControlEffector.spacecraft_idx`, `RobotArmReactionEffector.spacecraft_idx`, and `chaser_idx` and `target_idx` of `RPOMPCControlModel` and `PlannerGuidance` | `spacecraft_index`, `chaser_index`, `target_index` | warned from 0.3.0, `RobotArmReactionEffector`'s from 0.3.1; on these five types an old keyword is now a `MethodError` and an old property a `FieldError` (for the older names that still work on other types, see [Deprecated in 0.5.0](#Deprecated-in-0.5.0)) | 0.5.0 |
| selecting an RPO planner by a Symbol, such as `:rrt_star` | a planner object | warned from 0.3.0; a name is now a `MethodError` | 0.5.0 |
| a default random stream in the GNC planning functions | an explicit `rng` keyword | required now, without a window | |

A model of your own that names its spacecraft by a field must use the new names
too: from 0.5.0 the engine reads only `chaser_index`, `target_index` and
`spacecraft_index`. A model whose only index fields have the older names loses
these behaviours:

- with `chaser_idx` and `target_idx`, it is no longer seen as coupling two
  spacecraft, so the stiff solvers' Jacobian sparsity pattern leaves out the
  entries between them;
- with `chaser_idx`, `spacecraft_idx` or `sat_idx`, a control effector is no
  longer bound to its spacecraft by the external-propagation checks: it counts
  as acting on every spacecraft, and a run with an external propagator refuses
  it;
- with `spacecraft_idx`, a control effector is no longer bound to its spacecraft
  by the robot-arm and articulated-panel checks: the run does not couple an arm
  to the effector's plan, the panel-link check applies the effector to every
  spacecraft, and the viewer shows no arm plan for it.

Rename the fields; nothing else changes.

## Deprecated in 0.5.0

Older index names that were never deprecated follow the same convention now.
Each keeps working as a keyword or property until SpaceAGORA 1.0, with one
warning per type and name. An old keyword warns where it is passed, which is
before the run; an old property warns where it is read, which can be inside a
callback during the run. Under `julia --depwarn=error` either is an error at that
point.

| Item | Replacement | Status | Removed in |
|---|---|---|---|
| the `StateAnchor` property `sat_idx` | `spacecraft_index` | warns from 0.5.0 | 1.0 |
| the `sat_idx` keyword of `get_periapsis_pulse_callback` | `spacecraft_index` | warns from 0.5.0 | 1.0 |

An interlink checkpoint keeps its stored keys, so a checkpoint written before
0.5.0 resumes unchanged. The interlink names that 0.5.0 deprecates, from when
every link was a laser link, are listed on the
[Interlinks page](../user/interlinks.md#Deprecated-names).

## Changed in 0.6.0

One hook keeps its name and signature and is called for fewer spacecraft.

| Hook | 0.5.1 and earlier | 0.6.0 | What to do |
|---|---|---|---|
| `register_tstops!(integrator, model, spacecraft_index)` | called after each periodic control tick for every spacecraft index of the run | called for each spacecraft the model's `spacecraft_scope` names, in ascending order, as the tick is | a model that adds solver stops for a spacecraft outside its declared scope declares a scope that names it, or `AllSpacecraft()` |

A model without a `spacecraft_scope` method of its own has the default
`AllSpacecraft()`, which names every spacecraft, so its calls do not change. The
burn scheduler, the only built-in model with a `register_tstops!` method, is
event-driven and registers its stops from its own callbacks.

The scope itself is read once per periodic tick instead of once per spacecraft,
and twice for a control model that consumes a link: once before its producers
prepare it, and once after. So a model whose tick changes its own scope sees the
change from its next tick. A model that, ticking for spacecraft 1, retargets
itself to spacecraft 2 used to be ticked for spacecraft 2 in the same tick, and
is now ticked for it from the next one. A model that sets its spacecraft before
the run, in its constructor, in `bind_spacecraft` or in an initializer, is
unaffected. (0.7.0 goes further and stops a run whose scope changes; see
[Changed in 0.7.0](#Changed-in-0.7.0).)

## Removed in 0.6.0

The interlink scheduling policies chosen by a symbol were removed in 0.6.0,
without a deprecation period, together with `SchedulingPolicyModel`'s
`target_index` and `target_indices` and the names `target_idx` and
`target_idxs`, which had warned since 0.5.0. Scheduling policies are types now
(see [Interlinks](../user/interlinks.md#Policies)). A call with a symbol throws an
`ArgumentError` that names its typed replacement, whether or not Julia runs with
`--depwarn=error`. Interlink checkpoints written with the symbol policies still
resume, under the typed policy that schedules the same.

| Item | Replacement | Status | Removed in |
|---|---|---|---|
| `SchedulingPolicyModel(:gve_sma)`, with `target_index=k`, `target_index=nothing` or `target_indices=v` | `SemimajorAxisRatePolicy()`, `targets=k`, `targets=:both` or `targets=v` | an `ArgumentError` naming the replacement | 0.6.0 |
| `SchedulingPolicyModel(:gve_sma_first_index)` | `SemimajorAxisRatePolicy(targets=:lower_index)` | an `ArgumentError` naming the replacement | 0.6.0 |
| `SchedulingPolicyModel(:gve_eccentricity)`, with `target_index=k` or `target_index=nothing` | `EccentricityRatePolicy()`, `targets=k` or `targets=:both` | an `ArgumentError` naming the replacement | 0.6.0 |
| `SchedulingPolicyModel(:power_delivery)`, `SchedulingPolicyModel(:data_throughput)` | `PowerDeliveryPolicy()`, `DataThroughputPolicy()` | an `ArgumentError` naming the replacement | 0.6.0 |
| `SchedulingPolicyModel`'s keywords `target_idx` and `target_idxs` | as `target_index` and `target_indices` above | warned from 0.5.0; now an `ArgumentError` naming the replacement | 0.6.0 |
| `SchedulingPolicyModel`'s property `policy` | `objective`, the policy value | an `ArgumentError` naming the replacement | 0.6.0 |
| `SchedulingPolicyModel`'s properties `target_index`, `target_indices`, `target_idx` and `target_idxs` | `objective.targets` | an `ArgumentError` naming the replacement (the last two warned from 0.5.0) | 0.6.0 |
| the interlink states `laser_dv` and `laser_delta_sma` of each spacecraft's state | `interlink_dv` and `interlink_delta_sma` | an error from the state's own property lookup: the name of a state component cannot warn | 0.6.0 |

The columns of those two states keep their old names as aliases until 1.0 (see
[Names that cannot warn](#Names-that-cannot-warn)), and checkpoints keep the old
names inside, so checkpoints written before 0.6.0 resume.

## Changed in 0.7.0

A guidance, navigation or control model's per-spacecraft hooks follow its
`spacecraft_scope`, as its ticks already did. The hooks keep their names and
signatures and are called for fewer spacecraft.

| Hook | 0.6 and earlier | 0.7.0 | What to do |
|---|---|---|---|
| `control_wrench`, `control_propellant_rate`, `control_wheel_torque` | called at every evaluation of the equations of motion for every spacecraft | called for each spacecraft the model's scope names; for a `CoupledSpacecraft`, the primary only, since its partners are read, not acted on | nothing, if the hook returns zeros or `nothing` for other spacecraft, as the documentation asked; a model that acts on a spacecraft outside its declared scope, such as a coupled model that pushes its partners, declares a scope that names it, or `AllSpacecraft()`, and is ticked for it too |
| `actuator_levels`, `control_thruster_levels` | read for every spacecraft by the dynamics, the power and fault readers and the result columns | read during the run for each spacecraft the model's scope names; `control_thruster_levels` is still asked for every spacecraft before the run, when the result columns are sized | as above |
| `atmosphere_transition!` | called for every guidance and control model when a spacecraft crosses | called for the models whose scope names that spacecraft | as above |

The other changes of 0.7.0 to scopes:

- The scope is fixed for the run. The engine takes each model's scope when the
  run starts, once the run-start initializers have run, and stops the run when a
  later read differs: at each tick, and once per accepted step for an
  event-driven model. A model that changed its scope during a run, which 0.6
  ticked for its new spacecraft from its next tick (see
  [Changed in 0.6.0](#Changed-in-0.6.0)), sets its spacecraft before the run
  instead: in the model, in `bind_spacecraft` or in an initializer.
- A `SpacecraftSelection` scopes GNC models too. 0.6 refused it; it is now
  resolved when the run starts, with its `where` called once per spacecraft. Its
  `where` function is compared by identity, so keep the selection in the model
  and return that same object: a selection rebuilt at each call, with a new
  closure, stops the run.
- Direct callers are unchanged. The conformance kit and a model's own tests may
  still call the hooks for any spacecraft, so keep the zero and `nothing` returns
  outside the scope.
- The kit demo's `PDPointingController` and `AttitudeFilter` declare a narrower
  scope and return zero output outside it, so their runs do not change.
- A bound copy that acts on every spacecraft warns. A model declared on a
  spacecraft whose `bind_spacecraft` copy still declares `AllSpacecraft()` acts
  on every spacecraft, as at configuration level, and each run warns once for
  each such type (in an ensemble, each member's run). Return
  `BoundSpacecraft((k,))` from the copy's `spacecraft_scope`, for the copy bound
  to run position `k`.
- Ensembles refuse configuration-level GNC bound to run positions.
  `run_constellation_ensemble(args; allow_gnc_effectors=true)` refuses a
  configuration-level model whose scope is a `BoundSpacecraft`, a
  `CoupledSpacecraft` or a `SpacecraftSelection`. 0.6 passed such a model to every
  member, where it acted on every member's spacecraft or failed. Declare it on
  its spacecraft instead.
- The maneuver guidance models bind. `ApoapsisTargetPeriapsisRaiseGuidanceModel`
  and `AerobrakingCampaignPropulsiveManeuverGuidanceModel` gain
  `spacecraft_index` (0, meaning every spacecraft, by default) and may be
  declared on a spacecraft; their positional constructors without it still work.
  `BaseThrusterModel` stays one model for every spacecraft, and its refusal on a
  spacecraft says so.

## Schedule

- 0.3.0: the new names, with warnings for the old ones.
  `ExternalPropagationServices`, `HYPRServices` and `EDGServices` are at contract
  2.0.0, with SpaceAGORAMuJoCo 0.2.0, HYPR 0.2.0 and SpaceAGORAEDG 0.2.0.
- 0.3.1: `RobotArmReactionEffector.spacecraft_idx` becomes `spacecraft_index`,
  with a warning for the old name.
- 0.4.0: the companion packages use the new names. `skip_inactive` defaults to
  `true`: a guidance, navigation or control model stops ticking for a spacecraft
  once it is inactive, unless the model returns `false`. The pre-run
  configuration checks that warned in 0.3 refuse the run with an `ArgumentError`,
  whether or not Julia runs with `--depwarn=error`: a control period
  (`gnc_period_s`) that differs from its control rate or from its controller's
  discretization step, a spacecraft scope that names a position outside the run,
  and Apollo descent guidance and control that do not share their state and
  spacecraft.
- 0.5.0: the deprecated names are removed, and `control_propellant_rate` and
  `control_wheel_torque` gain defaults, `0.0` and `nothing`. The tagged releases
  0.3.0 to 0.4.0 carried the warnings, and the companion packages moved to the
  new names in 0.4.0.
- 0.6.0: the interlink scheduling policies chosen by a symbol, and
  `SchedulingPolicyModel`'s target names, are removed without a deprecation
  period (see [Removed in 0.6.0](#Removed-in-0.6.0)). A symbol is refused with
  its typed replacement, and old interlink checkpoints resume. The interlink
  states `laser_dv` and `laser_delta_sma` become `interlink_dv` and
  `interlink_delta_sma`; their `laser_*` columns stay as aliases until 1.0.
  `HYPRServices` and `EDGServices` move to contract 3.0.0, with HYPR 0.5.0 and
  SpaceAGORAEDG 0.5.0.
- 0.7.0: a guidance, navigation or control model's per-spacecraft hooks follow
  its `spacecraft_scope`, and a scope that changes during a run stops it (see
  [Changed in 0.7.0](#Changed-in-0.7.0)).

Release numbers follow `version` in `Project.toml`. The older `v1.0.0` and
`v2.0.0` git tags predate it.

## Checklist

1. Find each method of an old name with the command above. On 0.4.0, a test run
   with `julia --depwarn=error` also stops at each one. Search for
   `RPOGuidanceModel` too: its warning on 0.3 and 0.4 was an ordinary warning,
   which `--depwarn=error` does not turn into an error.
2. Rename each tick.
3. Move each output hook to its state-sample form, and flip the sign of the
   propellant rate.
4. Move each `calcForceTorque` method to `wrench`, with an
   `environment_requirements` method if it reads the atmosphere, the Sun or the
   planet-fixed frame.
5. Type the sample argument `x::StateSample`, or leave it untyped (see the
   [conventions](extension_conventions.md#The-supported-sample-argument)).
6. Run the tests again, and check that the pre-run check accepts every model.

## Names that cannot warn

The interlink names deprecated in 0.5.0 are listed under
[Interlinks](../user/interlinks.md#Deprecated-names). Two of them,
`ConstellationInterLinkModel` (now `InterlinkModel`) and `InterLinkParameters`
(now `LaserMomentumLink`), are aliases of the new types and cannot warn: code
that uses them runs silently, also under `julia --depwarn=error`. Search your
code for them before SpaceAGORA 1.0 removes them:

```bash
grep -rnE "ConstellationInterLinkModel|InterLinkParameters" src test
```

The columns `sc1_laser_dv_1` .. `sc1_laser_dv_3` and `sc1_laser_delta_sma`,
deprecated in 0.6.0, cannot warn either: until 1.0 a run writes them beside
`sc1_interlink_dv_1` .. `sc1_interlink_dv_3` and `sc1_interlink_delta_sma`, with
the same values. Search your analysis code for them:

```bash
grep -rnE "laser_dv|laser_delta_sma" .
```
