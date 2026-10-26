# Run an RPO case with interchangeable planners

!!! note "Optional package"
    The direct baseline planner runs with SpaceAGORA alone. The HYPR planner
    requires the optional `HYPR` package; see
    [Installing optional HYPR](#Installing-optional-HYPR).

The opt-in RPO pilot connects a planner to the existing simulation and LQ-MPC
controller. It supports truth observations, a passive target on a circular Earth
orbit, and a static station point cloud aligned with the target's rotating RTN
frame. The RPO research examples run on this path; for the guidance model they
used before, see [Moving from RPOGuidanceModel](#Moving-from-RPOGuidanceModel).

```julia
using SpaceAGORA
args = make_rpo_configuration(planner=DirectRPOPlanner(), seed=741)
solution = run_simulation(args; return_solution=true)
report = rpo_run_report(solution)
report.planners[1].records
report.commands[1].log
```

Load `SpaceAGORA` and `HYPR` in your project, then use `HYPRRPOPlanner(RPOPSOConfig(...))` as the `planner` to choose HYPR. The separate
`examples/rpo_planner_env` scripts prepare an environment outside the source
checkout and provide exact, deterministic settings for both choices through this interface. Its README documents setup and offline
execution. The default two-second case is a partial maneuver; the planners need
not return the same trajectory. It is not a completed rendezvous benchmark.

## The bounded scenario

The target is on a circular equatorial Earth orbit at 420 km. The default chaser
starts at `(3,0,0)` m and aims for `(5,0,0)` m in target RTN, with zero initial RTN
relative velocity. The run uses inverse-square gravity, no atmosphere and simple
ephemerides. The synthetic station has one point at the origin, a 0.25 m keepout
radius and a 0.15 m largest-chaser-half-extent allowance. Its geometry is a
point-sphere approximation, not watertight CAD geometry.

Reference limits are 0.2 m clearance, 0.5 m/s speed and `0.5*0.05/5.2` m/s²
acceleration. The 5 kg chaser carries 0.2 kg propellant and six 0.05 N thrusters
with 60 s specific impulse. Control runs every 0.1 s with a 12-step preview,
`Q=diag(20,20,20,2,2,2)`, `R=0.1I` and `Qf=10Q`. The default explicit Tsit5 solver
uses orbit/quaternion absolute and relative tolerances of 1e-8 and a 0.05 s cap.
Only this default `:tsit5` route is validated for the planner lifecycle. Other
solver modes may be accepted by the configuration type, but are outside this
pilot; some modes omit guidance and control callbacks. These are pilot settings,
not a certification of tracking or actuator feasibility.

Both adapters retain their reviewed 1% planning reserve and independent reference
validation. The direct baseline uses a cubic rest-to-rest law. Neither adapter
matches an arbitrary measured initial velocity as a boundary condition; the pilot
starts at zero relative velocity. Replanning can therefore introduce a reference
velocity discontinuity. General tracking acceptance remains separate.

## Initialization, replanning and expiry

The engine copies the configuration and prepares the initial plan from the actual
initial state before propagation. The planner receives a separate copied request;
the lifecycle validates against its retained snapshot. Validation checks identities,
frame, endpoints, uniform timing, finite values, configured limits and sampled
clearance. Only then is a complete plan installed for the controller.

Navigation, control and guidance retain their existing order. At a coincident
update the controller uses the preceding reference; guidance can install the next
one afterward. Ordinary and rejected dynamics evaluations only read held control
commands and do not invoke planning. A scheduled request is delivered at the first
accepted guidance tick at or after its requested time:

```julia
events = (RPOPlanningEvent(0.5; reason=:forced, token=1),
          RPOPlanningEvent(1.0; reason=:replan, token=2))
args = make_rpo_configuration(planner=DirectRPOPlanner(), planning_events=events)
```

A token names one event. Duplicate delivery is ignored; different tokens at the
same tick are distinct requests. Optional positive `replan_interval_s` or
`tracking_error_limit_m` enable background requests. These policies operate on the
configured geometry; a geometry provider supplies geometry that changes during the
run (see [Changing geometry, holds and retiming](@ref)). Periodic and scheduled requests use the same configured
absolute time tolerance at tick boundaries. There is no guidance tick at the
mission end. Events without an eligible guidance tick before the end remain
undelivered; inspect the recorded requests rather than assuming every configured
event ran.

By default, every failed initial, forced, background or retiming request stops this
route with `RPOPlanningError`, which carries copied diagnostics, spacecraft/request
IDs, and any original exception and backtrace. There is no implicit retain-old-plan
or emergency-controller fallback; `on_rejection=:install_fallback` opts into an
explicit hold instead. An explicit `:retime` event requires the planner's retiming
capability; `HYPRRPOPlanner` advertises it and `DirectRPOPlanner` does not.

Immediately before control reads a reference, the lifecycle checks its lifetime
through the entire preview. Repeating the last reference column cannot bypass
expiry. `plan_validity_s` is a request lifetime, not the duration of the maneuver.
A short plan may repeat its final sample while still valid. A long plan must fit
inside its lifetime. Reference validity does not guarantee a stationary vehicle.

The constructor accepts a lifetime as short as `12*control_dt_s`, but the first
control update is at `control_dt_s` and needs at least `13*control_dt_s` of initial
validity. For an interval aligned with guidance ticks, periodic replanning needs
at least `replan_interval_s + 12*control_dt_s`, because control reads the old plan
before the replacement. Other intervals must cover the next eligible guidance
tick plus the preview. Without replanning, validity must reach the last control
update plus the preview. These bounds describe consumption timing; a planner may
require a longer lifetime to fit its maneuver.

## Changing geometry, holds and retiming

Three opt-in keywords of `make_rpo_configuration` extend the loop above. Each
request still carries one owned geometry snapshot, and every installed reference
still passes the same validation.

`geometry_provider` supplies the geometry when it changes during the run. The
built-in `RPOReplanningSpheres` adds spheres that appear and disappear to the
station geometry:

```julia
spheres = RPOReplanningSpheres([(center_rtn=(4.0, 0.0, 0.0), radius_m=0.25,
                                 appear_time_s=0.2, disappear_time_s=0.8)])
args = make_rpo_configuration(planner=HYPRRPOPlanner(config), geometry_provider=spheres)
```

Each geometry has a revision string; equal revisions mean equal geometry. When the
revision changes, the remaining path of the installed reference is checked against
the new geometry with the same clearance query and sampling as validation. If it
still clears, the reference is kept; otherwise a `:replan` request follows. To
write a provider, subtype `AbstractRPOGeometryProvider` and implement
`rpo_geometry_revision(provider, revision, time_s)` and
`rpo_geometry_snapshot(provider, geometry, time_s)`; both receive the configured
station geometry or its revision.

`on_rejection=:install_fallback` replaces a rejected request with a hold: a
reference at the chaser's position with zero velocity, which the controller brakes
onto. A hold is not validated and does not reach the goal. It is retried once with
a `:replan` request when the chaser settles on it (`hold_position_tolerance_m`,
`hold_speed_tolerance_mps`), renewed by a request before it expires, and replanned
when the geometry revision changes. Interval and tracking-error replans are
suspended while it holds, so a planner that keeps refusing is asked again only by
these triggers or a planning event. The default, `on_rejection=:fail`, stops the
run as before.

`retime_tracking_error_m` requests `:retime` when the tracking error exceeds it and
stays within `tracking_error_limit_m`. It needs a planner that declares retiming.
`HYPRRPOPlanner` retimes the installed reference's validated path from the point
nearest the chaser, with the same limits as a new plan; a retime cannot reach a
changed goal, so it is then refused. Validation can also refuse a retimed
reference, for example where HYPR's tangential profile turns faster than the
common acceleration limit; with `on_rejection=:install_fallback` the chaser then
holds.

The run records mark each of these requests, and each kept geometry revision, with
`event=:trigger`, and each hold with `event=:fallback_installed`.

`examples/Earth_RPO_Planner_Disturbances.jl` combines them with
`DirectRPOPlanner`, in the examples environment. A sphere appears on the remaining path,
so the replan is refused and the chaser holds; once it settles, the planner is
asked once more and refuses again. When the sphere disappears, the new revision
replaces the hold with a reference to the goal. A displacement later raises the
tracking error above `tracking_error_limit_m`, and the same guidance tick replans
from the displaced state. The script prints this timeline from `rpo_run_report`.

## Ownership, reproducibility and inspection

Each controlled spacecraft owns fresh planner state, an independent controller
workspace, a reference slot and a MersenneTwister stream. Stream seed words are the
eight little-endian UInt32 words of SHA-256 over the ASCII string
`SpaceAGORA/RPO/v1/<seed>/<chaser_id>/<target_id>`. They do not depend on spacecraft
order, callbacks, threads or Julia's randomized `hash`. Reusing a configuration
starts fresh state. Deterministic parity is scoped to the recorded Julia and
package versions; wall-clock stopping and cross-version bit identity are excluded.

`chasers` accepts named tuples containing `id`, `start_rtn_m` and `goal_rtn_m`.
They share the passive target specified by `target_id` but own distinct runtime
state. The geometry checks the station only; inter-chaser collisions are not
modeled by this interface.

`rpo_run_report(solution)` returns copies of requests, validation, installed
references, planner settings, seed mapping, held-command samples and saved states.
Reading and changing this report cannot change the run. The returned solution
itself still exposes the usual advanced solver internals; do not mutate them while
running. Reference acceptance is distinct from physical tracking acceptance.

`isolate_state=false` is refused before propagation or output creation. A
checkpoint holds the planner binding's run state (its random stream, installed
request and reference, records and trigger state), the plan buffer and the
controller's held command, warm start and solver iterate, so a run resumed from
it continues as the uninterrupted run would. A planner that keeps state of its
own (`initialize_planner` returns something) runs with checkpoints only if it
declares `restart` in `planner_capabilities`; the built-in planners keep none.
A core-only installation can run the direct planner without HYPR; see the optional
installation instructions below.

### Checkpoints

The table lists what a checkpoint covers in any run, the pilot's included.
[Checkpointing](outputs.md#Checkpointing) describes the settings and the run
state a checkpoint file holds, and
[Run state in checkpoint files](../extensibility.md#Run-state-in-checkpoint-files)
how a model adds its own.

| | What |
|---|---|
| Supported: a resumed run continues as the uninterrupted run would | the integrated state; which spacecraft are active; the run state the engine keeps on the spacecraft (batteries, actuator faults, sensors, hinges and gimbals, and the ramp factor and held thrust of mounted thrusters with firing dynamics); the times of the periodic guidance, navigation and control ticks: a tick and a checkpoint boundary that differ only by rounding count as one time (a 0.1 s tick period with 0.3 s checkpoints, for example), a checkpoint records the ticks due at its time that were not yet taken, also when the run stopped early there, and the resumed run takes them first, once, and a checkpoint written by SpaceAGORA 0.4.0 just before a tick resumes with that tick taken at the resumed run's start; the joint state of articulated and compliant vehicles; interlinks, with their topology, link kinds and parameters; guidance, navigation and control models that save their run state with `checkpoint_state`: the RPO planner binding with its plan buffer, the RPO MPC controller, the burn scheduler, the momentum manager, the Apollo descent guidance and control, the periapsis-raise guidance, the robot-arm controller, and models of your own that implement the two methods |
| Refused before the run, with an `ArgumentError` | energy-depletion guidance and control; `AlgorithmGuidance`; an RPO planner that keeps state of its own and does not declare `restart`; `run_mission`; externally propagated spacecraft; resuming with a GRAM density perturbation mode on; a model whose `checkpoint_state` returns data and that has no `restore_checkpoint_state!` method; resuming from a checkpoint whose recorded configuration differs from the run's: the type name and position of each sensor and, when a model or sensor saves state with `checkpoint_state`, of each guidance, navigation and control model (the short type name, without its module or type parameters: a same-named type from another module, or another parametrization of the same type, is not told apart), the batteries, faults, hinges and gimbals, the mounting place of each thruster with firing dynamics (its position in the mounting order, its link's position and its position among that link's actuators) and its `max_thrust`, `cutoff_frequency`, `min_firing_time`, `level_on` and `level_off`, and the interlinks (a change to a field the checkpoint does not record is not detected); resuming from a checkpoint whose guidance, navigation or control models act on other spacecraft (each model's `spacecraft_scope` after the run-start initializers), or are another number in a group, whether or not they save state; resuming a run with guidance, navigation or control models from a checkpoint written by SpaceAGORA 0.5.1 or earlier, which holds no such record; resuming a run with firing-dynamics thrusters from a checkpoint without their firing state, such as one written by SpaceAGORA 0.4 or earlier |
| Outside the guarantee: not refused, but a resumed run may differ | the orbit counter, which the checkpoint does not hold, so `MissionOrbits` and guidance that counts orbits count again from the resumed run's start; a spacecraft already inactive when the run resumes, for which a model other than the burn scheduler starts from its restored command: with `skip_inactive` `true` and a `release_command!` that clears the command, a nonzero restored command can stay visible until the model's first tick releases it; the default `release_command!` does nothing, so it does not clear it, and a model with `skip_inactive` `false` keeps ticking as for an active spacecraft (the trajectory does not change, since an inactive spacecraft is not propagated); callbacks a model returns from `gnc_callbacks`, beyond the state the model saves (see [Writing a GNC Algorithm](gnc_algorithm.md#Declarations)); the mutable state of actuator kinds of your own, which have no checkpoint hook; the periodic ticks of the `:multirate` solver mode, which repeat because its sub-solves revisit intervals, with or without checkpoints; solver modes other than `:tsit5`, `:dp8` and `:gravity_backbone_split`, the three the resume tests run with (the `:gravity_backbone_split` tests resume translational runs, since that mode integrates no attitude and runs no guidance, navigation or control models) |

## Add your own planner or force

Subtype `AbstractRPOPlanner`. Implement `planner_capabilities`, optionally
`initialize_planner`, and `plan_rpo!(state, planner, request, rng)`. Return an
`RPOPlanningResult` containing an `RPOReference` for candidate status. Reference
times are uniform and relative to `origin_time_s`; vectors are SI target RTN,
including the rotating-frame transport term in velocities. The
[RPO planner contract](../reference/rpo_planner_contract.md) lists the validation
rules. Add `retime_rpo!` only when you also declare retiming.

An RPO planner is also an algorithm of SpaceAGORA's general GNC algorithm
contract, the first family of it.
[The general algorithm contract](../reference/rpo_planner_contract.md#The-general-algorithm-contract)
maps the RPO names to the general ones; the RPO names remain the supported way
to write an RPO planner.

A custom planner must support `deepcopy`, including when copied per spacecraft,
per run and into failure reports. `initialize_planner` must return fresh state
for each call. Use only the supplied `rng` for randomness; keep mutable runtime
state in that returned state, with no module-level mutable state or background
tasks. The lifecycle isolates cooperating implementations; it is not a sandbox
for arbitrary planner code.

The separate-project smoke defines a tiny test planner using only these public
methods, including an explicit failed result. Its force example uses the existing
`AbstractForceTorqueModel` and `wrench` interface
([Adding a Force or Torque](custom_effector.md)), passed through
`extra_effectors=(model,)`. No internal package module or source include is needed.

## Installing optional HYPR

`using SpaceAGORA` loads the direct planner, simulation lifecycle, shared geometry,
sampling, metrics, RRT, retiming and quintic robot-arm planner. It preserves HYPR
configuration and result types so saved configurations and qualified names retain
their defining modules. It does not load the HYPR optimizer or configured execution.

From a reviewed SpaceAGORA Git clone, prepare a separate environment:

```sh
julia --startup-file=no scripts/optional_packages.jl setup hypr /path/to/my-rpo-project
```

The helper checks out `packages/robotics/HYPR.jl` at the parent repository's
recorded revision, verifies its identity and clean state, and delegates to the
companion's `scripts/setup_env.jl` with that path and this core checkout. The pin agrees with
`packages/HYPRSource.toml`. It selects HYPR 0.5.0 alongside
the supported SpaceAGORA source. HYPR remains optional for the simulator core;
there is no separate compatibility loader to install.

The older command remains as a compatibility shortcut to the same installer:

```sh
julia --startup-file=no scripts/setup_hypr.jl /path/to/my-rpo-project
```

Without an override, this shortcut checks out the same pinned submodule. For
coordinated local development, `SPACEAGORA_HYPR_PATH` explicitly substitutes
another local HYPR checkout; that override is not published-installation evidence.
Failed resolution is an error, not a fallback to another revision. HYPR, EDG and
MuJoCo use the same shared installation policy: the destination must be outside
both source checkouts, and existing unmanaged directories are refused.

To preserve an older HYPR environment, migrate its project, manifest and local
preferences into a fresh destination with
`scripts/setup_hypr.jl --migrate-from /path/to/old /path/to/new`. This leaves the
old directory unchanged, selects the current core and HYPR together, and removes
the retired `SpaceAGORAHYPR` loader. Migration refuses other packages that still
require that loader; update those packages first. See [Optional package commands](installation_environment.md#Optional-package-commands)
for source checks and Git-clone requirements.

Run with `--project=/path/to/my-rpo-project` and explicitly load the companion.

```julia
using SpaceAGORA, HYPR
@assert hypr_available()
planner = HYPRRPOPlanner(RPOPSOConfig())
```

When the planning configuration carries an acceleration limit, also set
`retime_accel_limit_enable=true` in `RPOPSOConfig`; see
[Built-in planners](../reference/rpo_planner_contract.md#Built-in-planners).

Either package may be loaded first. Without the companion, HYPR selection raises
`HYPRUnavailableError` with installation instructions. Public RPO assembly and
simulation preflight perform this check before creating run outputs.
Robot-arm `planner=:hypr` also requires explicit loading; `:cloth_quintic` does not.

Comparison planners that request HYPR-configured objectives, refinement, sampling
or retiming need the companion even when their search algorithm is RRT, CHOMP or
STOMP. The seven comparison modes are also planner objects, such as
`RRTStarPlanner(...)`, which can be the `planner` of `make_rpo_configuration`; see
[Comparison planners](../reference/rpo_planner_contract.md#Comparison-planners).
Selecting a comparison mode by name was removed in SpaceAGORA 0.5.0. Shared RRT with explicit policies and shared scalar APIs remain independent.
The baseline LQ-MPC controller still needs OSQP; optional HYPR does not remove that
controller dependency. Package separation changes installation, not the validated
physical limits, planning reserve or numerical policies.

The `examples/rpo_planner_env/setup.jl` command above installs both components
explicitly for the two-planner demonstration in an external managed environment.
`examples/rpo_planner_env/run.jl` selects that environment for the RPO and compliant
robot-arm HYPR examples and ISS viewer demo. They load the companion before their
common helper continues in the examples environment,
which `julia --startup-file=no examples/setup.jl` prepares. See
[Examples Catalog](examples_catalog.md) for the exact commands. Baseline examples
use `--project=examples` and do not preload HYPR. Do not add the companion as a required
dependency of SpaceAGORA itself.

## Moving from RPOGuidanceModel

`RPOGuidanceModel`, the HYPR guidance model of the older RPO examples, was
removed in SpaceAGORA 0.5.0; 0.4.0 is the last release that has it. Build RPO
runs with `make_rpo_configuration` and a planner. For HYPR planning, pass
`planner=HYPRRPOPlanner(RPOPSOConfig(...))`; the configuration then uses a
`PlannerGuidance`, which validates each reference before installing it.

The RPO examples run on this path: the main example plans inside the run. The
reconstruction example, now in the SpaceAGORARPOStudies project,
flies each comparison planner's plan with `RetainedPathRPOPlanner`, which times a
given path for the run's request so that the run can validate it. The replanning
study of the reconstructed disturbance families, including a goal change and the
model's brake-then-replan policy, ran on the model and was removed with it; that
project's reconstruction guide (SpaceAGORARPOStudies) says how to reproduce it
with SpaceAGORA 0.4.0.
`examples/Earth_RPO_Planner_Disturbances.jl` shows disturbances on the planner
path.

The planner path covers the model's replanning features:

| `RPOGuidanceModel` | `make_rpo_configuration` |
| --- | --- |
| `replanning_config` spheres | `geometry_provider=RPOReplanningSpheres(spheres)` |
| safe hold when replanning fails | `on_rejection=:install_fallback`, `hold_position_tolerance_m`, `hold_speed_tolerance_mps` |
| retiming on tracking error | `retime_tracking_error_m`, or a `:retime` planning event |
| `safe_distance_m` | the clearance in `constraints` |
| `force_replan` | `RPOPlanningEvent(time_s; reason=:forced)` |

The two differ in policy. The model brakes to a hold before every replan and
retimes when clearance falls below a margin; the planner path replans at once and
holds only when a request is rejected. A run on the planner path also draws from a
different random stream, and its validation can refuse a reference the model would
install, so its trajectories are not identical to the model's.
