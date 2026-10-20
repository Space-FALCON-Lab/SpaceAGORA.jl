# Migrating Extension Hooks

SpaceAGORA's extension hooks follow one set of [conventions](extension_conventions.md):
time before the spacecraft index, a state sample for right-hand-side hooks, and
positive consumption. Ten hook names, and one form of the plume model's
`wrench`, did not. SpaceAGORA 0.3.0 deprecated them and 0.5.0 removed them.
In 0.3.0 and 0.4.0, the pre-run check warned about a model's method of an old
name once per process for each model type and name, a direct call of an old
name warned too, and the five-argument plume form warned only when called
directly. This page lists each name with its replacement.

A method of a removed name no longer runs. One that names the package, such as
`SpaceAGORA.calcForceTorque(...) = ...`, fails when its file loads, because the
name is gone. The pre-run check refuses a model left without a method of a
required replacement, a tick, `control_wrench` (unless the effector's actuation
is `:mounted`) or `wrench`, and names the hook it needs. `control_propellant_rate`
and `control_wheel_torque` are optional: without a method of its own, an
effector gets `0.0` and `nothing`, and the check does not refuse it. An
ephemerides model without the model-first hooks stops with a `MethodError`
naming the hook when the run first asks it for its epoch or frame. To find the
old names in your code:

```bash
grep -rnE "calcForceTorque|calcControl(ForceTorque|MassFlowRate|Effect!)|calcReactionWheelTorque|calc(Guidance|Navigation)Effect!|planet_frame_lpi|ephemerides_(sun_direction_ii|time_seconds)" src test
```

## Removed hook names

| Removed in 0.5.0 | Replacement | What changes |
|---|---|---|
| `calcForceTorque(model, x, p, i)` | `wrench(model, x::StateSample, env::EnvironmentSample, t::Float64)` | a state sample instead of the raw view; the environment the engine sampled instead of `p`; the time as an argument and the index in `x.spacecraft_index` |
| `wrench(model::PlumeSurfaceInteractionModel, x, env, t, i)` | `wrench(model, x, env, t)` | the index in `x.spacecraft_index` instead of an argument |
| `calcControlForceTorque(model, u, p, i, t)` | `control_wrench(model, x::StateSample, p, t::Float64)` | a state sample instead of the raw view, the time before the index, the index in `x.spacecraft_index` |
| `calcControlMassFlowRate(model, u, p, i, t)` | `control_propellant_rate(model, x::StateSample, p, t::Float64)` | the same, and **the sign** |
| `calcReactionWheelTorque(model, u, p, i, t)` | `control_wheel_torque(model, x::StateSample, p, t::Float64)` | the same as for the wrench |
| `calcControlEffect!(model, u, p, t, i)` | `update_control!(model, u, p, t, spacecraft_index)` | the name only |
| `calcGuidanceEffect!(model, u, p, t, i)` | `update_guidance!(model, u, p, t, spacecraft_index)` | the name only |
| `calcNavigationEffect!(model, u, p, t, i)` | `update_navigation!(model, u, p, t, spacecraft_index)` | the name only |
| `planet_frame_lpi(planet, et, model)` | `ephemerides_frame_lpi(model, planet, et)` | the model first |
| `ephemerides_sun_direction_ii(planet, et, model)` | `ephemerides_sun_ii(model, planet, et)` | the model first, and the name |
| `ephemerides_time_seconds(initial_time, model)` | `ephemerides_epoch_seconds(model, initial_time)` | the model first, and the name |

!!! warning "The sign changes"
    `calcControlMassFlowRate` returned the spacecraft's mass derivative,
    **negative** while propellant is used. `control_propellant_rate` returns the
    propellant used, **positive**. Flip the sign when you rename the method: a
    method renamed without it makes the vehicle gain mass.

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
`spacecraft_index`. A model whose only index fields are older names no longer
gets these behaviours:

- with `chaser_idx` and `target_idx`, it is no longer seen as coupling two
  spacecraft, so the stiff solvers' Jacobian sparsity pattern leaves out the
  entries between them;
- with `chaser_idx`, `spacecraft_idx` or `sat_idx`, a control effector is no
  longer bound to its spacecraft by the external-propagation checks: it counts as
  acting on every spacecraft, and a run with an external propagator refuses it;
- with `spacecraft_idx`, a control effector is no longer bound to its spacecraft
  by the robot-arm and articulated-panel checks: the run does not couple an arm to
  the effector's plan, the panel-link check applies the effector to every
  spacecraft, and the viewer shows no arm plan for it.

Rename the fields; nothing else changes.

## Deprecated in 0.5.0

Three older index names that were never deprecated follow the same convention.
Each keeps working as a keyword or property until SpaceAGORA 1.0, with one warning
per type and name. An old keyword warns where it is passed, which is before the
run; an old property warns where it is read, which can be inside a callback
during the run. Under `julia --depwarn=error` either is an error at that point.

| Item | Replacement | Status | Removed in |
|---|---|---|---|
| `SchedulingPolicyModel`'s keywords and properties `target_idx` and `target_idxs` | `target_index` and `target_indices` | warn from 0.5.0 | 1.0 |
| the `StateAnchor` property `sat_idx` | `spacecraft_index` | warns from 0.5.0 | 1.0 |
| the `sat_idx` keyword of `get_periapsis_pulse_callback` | `spacecraft_index` | warns from 0.5.0 | 1.0 |

An interlink checkpoint keeps its stored keys, so a checkpoint written before
0.5.0 resumes unchanged. The interlink names that 0.5.0 deprecates, from when every
link was a laser link, are listed on the
[Interlinks page](../user/interlinks.md#Deprecated-names).

## Schedule

- **0.3.0.** The new names, with warnings for the old ones.
  `ExternalPropagationServices`, `HYPRServices` and `EDGServices` are at
  contract 2.0.0, with SpaceAGORAMuJoCo 0.2.0, HYPR 0.2.0 and SpaceAGORAEDG 0.2.0.
- **0.3.1.** `RobotArmReactionEffector.spacecraft_idx` becomes `spacecraft_index`,
  with a warning for the old name.
- **0.4.0.** The companion packages use the new names. `skip_inactive` defaults
  to `true`: a guidance, navigation or control model stops ticking for a
  spacecraft once it is inactive, unless it returns `false`. The pre-run
  configuration checks that warned in 0.3 refuse the run with an `ArgumentError`,
  whether or not Julia runs with `--depwarn=error`: a control period
  (`gnc_period_s`) that differs from its control rate or from its controller's
  discretization step, a spacecraft scope that names a position outside the run,
  and Apollo descent guidance and control that do not share their state and
  spacecraft.
- **0.5.0.** The deprecated names are removed, and `control_propellant_rate` and
  `control_wheel_torque` gain defaults, `0.0` and `nothing`. The tagged releases
  0.3.0 to 0.4.0 carried the warnings, and the companion packages moved to the
  new names in 0.4.0.

Release numbers follow `version` in `Project.toml`. The older `v1.0.0` and
`v2.0.0` git tags predate it.

## Checklist

1. Find each method of an old name with the command above. On 0.4.0, a test
   run with `julia --depwarn=error` also stops at each one. Search for
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
6. Run the tests again: the pre-run check accepts every model.

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
