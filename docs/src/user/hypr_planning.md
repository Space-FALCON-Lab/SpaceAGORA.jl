# Motion planning with HYPR

!!! note "Optional package"
    Requires the optional `HYPR` package; see
    [Install and try the planner](#Install-and-try-the-planner).

HYPR plans a route through a cluttered environment and turns it into a
reference that a spacecraft controller can follow. Its rendezvous and
proximity operations (RPO) implementation in SpaceAGORA is meant for motion
near a target spacecraft or station. The planner chooses the route and its
timing; the simulation then shows how well the controlled vehicle follows it
and how much propellant the maneuver uses. These are separate questions. A low
planning cost says something about a candidate route, but it does not measure
the fuel consumed while tracking it.

Jakob Gorisek-Gazze and Giusy Falcone describe the method in *HyPR: Safety-
and Fuel-Aware Motion Planning and Replanning for Autonomous SmallSat
Servicing*, accepted to the IEEE International Conference on Space Robotics
(iSpaRo), 2026. The [Space-FALCON Lab publication record](https://www.spacefalconlab.com/publications/)
gives the citation. This page keeps that study apart from the current software
interfaces and from the small examples used to check an installation.

## From a route to a controlled maneuver

The translational method searches for a feasible route with RRT-Connect and
uses it as the initial guess for a particle swarm whose particles are candidate
Bézier paths. The swarm improves the geometric path against the configured
obstacle and fuel-related objective. In the software the RRT-Connect warm start
is optional: it runs when `RPOPSOConfig` sets `rrt_warmstart_enable=true` (the
default is `false`), and `HYPRRPOPlanner(config; rrt_on_replan=true)` also uses
it for `:replan` requests. In `:manuscript` mode the warm start's detour and
iteration count set the search effort and coefficient schedules. The returned
path is the best candidate evaluated, not the swarm's final particle positions.

Retiming assigns speeds and times to the chosen geometry. Clearance, reaction
time and acceleration limits bound how fast the vehicle should traverse each
part of the route. These acceleration-based speed caps are separate from the
optional forward and backward passes that bound tangential acceleration, which
`retime_accel_limit_enable=true` turns on; they are off by default.
SpaceAGORA's LQ-MPC controller then tracks the reference with its own
prediction horizon, weights and actuator model, so report the planner's fuel
proxy and the tracked maneuver's propellant use separately. Changing the
controller can change the propellant use even when the geometric path stays
the same.

Planning coordinates are in the target's rotating radial, transverse and
normal (RTN) frame, with distances in meters and time in seconds. Relative
velocity in that frame includes the rotation transport term. The station
geometry and the chaser's extent together determine the clearance, so the
distance from the chaser's center to a station sample is not, by itself, the
available collision margin. Keep the geometry, sampling resolution and safety
distance with the experiment settings when you compare results.

## Install and try the planner

The search implementation lives in the separate
[HYPR.jl repository](https://github.com/Space-FALCON-Lab/HYPR.jl). SpaceAGORA
provides the simulation, the shared geometry and the versioned integration
services. Load both `SpaceAGORA` and `HYPR`; loading SpaceAGORA alone does not
activate HYPR.

From a SpaceAGORA Git clone, prepare a separate environment and run the public
planner example:

```sh
julia --startup-file=no scripts/optional_packages.jl setup hypr /path/to/hypr-env
JULIA_LOAD_PATH='@:@stdlib' julia --startup-file=no --project=/path/to/hypr-env \
  examples/rpo_planner_env/smoke.jl
```

The `JULIA_LOAD_PATH` setting lets the example load the standard `Test`
library without picking up packages from your default environment. The setup
helper checks out the submodule at `packages/robotics/HYPR.jl` at the revision
recorded in `packages/HYPRSource.toml` and verifies it. HYPR sits in the
robotics folder only as a matter of source organization; it plans spacecraft
RPO as well as robot-arm motion. The
[installation guide](installation_environment.md#Optional-package-commands)
describes the source checks, and the
[RPO guide](rpo_planner_pilot.md#Installing-optional-HYPR) also covers the
`scripts/setup_hypr.jl` compatibility shortcut to the same companion installer. Record the HYPR pin and the SpaceAGORA revision
with each run.

The short example runs both the direct planner and HYPR through the public
interface and checks that their references reach the controller and produce
finite, bounded outputs. It covers only the first two seconds of a maneuver.
Its small swarm, simple geometry and short controller preview make it an
installation check, not the paper's experiment.

For a new application, load `SpaceAGORA` and `HYPR`, build an `HYPRRPOPlanner`
from an explicit `RPOPSOConfig`, and pass it as `planner` to
`make_rpo_configuration`. Set `hypr_mode=:manuscript` in the configuration when
you investigate the manuscript formulation; the default is `:legacy`. The
[RPO planner guide](rpo_planner_pilot.md) explains the complete call sequence,
reference validation, seed mapping and the information returned by
`rpo_run_report`.

## Relating a run to the paper

The accepted iSpaRo paper describes a 500-case static comparison and, for the
dynamic study, 100 baseline paths under eight additional scenarios, or 800
dynamic cases per planner. Its static comparison uses MersenneTwister seeds
`740 + i`, a 0.5 m safety distance, a 0.1 s controller update and a 60-step
MPC horizon, with controller weights `Q = diag(20,20,20,2,2,2)`, `R = 0.1I`
and `Qf = 20Q`. Competing planners get a per-case time budget based on HYPR's
measured planning time, so a different processor or thread count can change
both the reported runtime and how much optimization a competing planner
performs.

The repository's static and dynamic research entry points are
`examples/Earth_RPO_CubeSat_MPC_PlannerComparison.jl` and
`examples/Earth_RPO_CubeSat_MPC_Replanning.jl`. Their defaults are useful
configurations, but they are not a transcription of the paper's Table I. The
comparison defaults use a different MPC horizon and weights, keep the
`:legacy` planner mode and leave acceleration-limited retiming off. The batch
driver also derives each planner's seed from a master random stream instead of
passing `740 + i` to every planner. Keep the actual endpoints, sampled
geometry, random states and every effective setting; a seed number alone does
not identify an experiment.

The tracking model differs as well. The comparison helper propagates a
discrete linear Hill-Clohessy-Wiltshire model and accounts for the continuous
commanded impulse, whereas the paper describes nonlinear two-body propagation
with a 20% pulse-width-modulation activation threshold. Matching the weights
of Table I does not remove that difference. A run of this helper is useful for
checking planning and tracking, but on its own it cannot establish agreement
with the paper's fuel statistics.

The paper computes static fuel statistics over successful runs. Its matched
comparison uses the same successful cases for each competing planner and for
HYPR; comparing unrelated sets of successful cases answers a different
question. The dynamic results use a different baseline, the HYPR static-map
fuel consumption for that experiment, and neither quantity can stand in for
the other. Keep unsuccessful cases in the success-rate denominator and keep
their failure reasons.

The manuscript-mode tests, in SpaceAGORA's test suite (for example
`test/unit/gnc/rpo_hypr_manuscript_tests.jl`) and in HYPR's, check selected
equations, seeded search behavior and integration contracts. They do not
establish agreement with the paper's aggregate results. A reproduction also needs the paper's experiment inputs,
scenario definitions and complete outcome records. Until that comparison
exists, cite the paper's numbers as results reported in the manuscript, not as
measurements of the current release.

## Replanning and the supported interface

The paper and the research replanning example distinguish three responses to
a meaningful change in clearance, tracking error or goal: keep the reference,
retime it, or search for a new route. The public planner interface
(`make_rpo_configuration`) uses a static scene and explicit planning requests
by default, and it records why each request was accepted or rejected. Opt-in
keywords add geometry that changes during the run (for example
`RPOReplanningSpheres`, spheres that appear, move and disappear), background
replanning on an interval or on tracking error, retiming, and a hold when a
request is rejected; the [RPO planner guide](rpo_planner_pilot.md#Changing-geometry,-holds-and-retiming)
describes them. None of this reproduces the paper's complete dynamic
experiment.

The research replanning example decides with `rpo_replanning_decision`, using
the settings of an `RPOReplanningConfig`; both live in
`SpaceAGORA.SimulationModel.GuidanceHooks`. A requested goal
(`desired_goal_rtn`) more than `goal_change_tolerance_m` from the reference's
end, a remaining path whose clearance falls below `safe_distance_m`, or a
tracking error above `tracking_error_replan_m` leads to a replan. Clearance
below `retime_clearance_m` or a tracking error above `tracking_error_retime_m`
leads to a retime. The tracking error is the distance to the future part of the
reference (line 7 of the paper's Algorithm 1). The remaining path is resampled
segment by segment, so a coarse sampling interval cannot drop a corner and
create a shortcut. An expired reference, including a multipoint one, is not
retimed again merely because its clearance lies in the retiming band, which
would keep resetting its completion clock. Goal changes, unsafe clearance and
excessive tracking error are still checked first, and expiry is measured
against the absolute reference deadline so that it behaves exactly at callback
boundaries.

Retiming keeps the geometric route and rebuilds only its timing. Following
lines 14 and 15 of Algorithm 1, `rpo_retime_existing_plan` projects the chaser
onto the retained route, connects the measured position to the next route
vertex and times the whole connected path with the clearance and curvature
speed law; `HYPRRPOPlanner`'s `retime_rpo!` uses the same route rule. Time
samples on a straight stretch are dropped only when every removed point lies
within 1e-10 m of a monotone retained segment. Real bends, reversals and
excursions between coincident endpoints stay as vertices, so a small tracking
offset cannot create a sharp turn at the next time sample. A new plan keeps
the geometric polyline sampled from the planner's actual curve before time
integration changes the sample spacing, and later retimes use that polyline
rather than earlier measured-state connectors. Explicit polyline corners are
kept even when they are closer together than the sampling interval. A plan
without this record uses its delivered reference as the polyline; Bézier
control points are never inferred from timed samples. The retained geometry
must be finite 3-by-N points ending at the reference goal. The connector from
the measured position takes traversal time instead of overwriting the first
position, and a measured initial velocity passed to `rpo_retime_existing_plan`
becomes the reference's first velocity rather than being clipped to a newly
tightened cap.

For each run, inspect the installed reference, its sampled clearance checks
and the executed trajectory. A validated reference does not guarantee
collision-free tracking under arbitrary disturbances. A suddenly appearing
obstacle can leave the chaser in a state from which the clearance cannot be
kept with the available thrust; count such a case as a failed maneuver.
The controller, actuator limits and failure policy remain part of the
application. The [planner contract](../reference/rpo_planner_contract.md) and
[HYPR service reference](../reference/hypr_services.md) describe these
boundaries for developers extending the integration.

## Reconstructing the paper's experiment

The reconstruction of the final paper's experiment, with its pinned
environment, commands and guide, lives in the private SpaceAGORARPOStudies
repository. Lab members with access use it there; outside collaborators need
an invitation. It runs newly generated cases through the full SpaceAGORA
simulator with the paper's clearance, control rate, horizon and controller
weights. The original case coordinates, geometry draws and the source
revision that produced the paper were not recovered, so it is a
reconstruction, not a replay of the original runs. Planning times also depend
on the hardware and thread count used.

The public pieces it builds on are in core. `RetainedPathRPOPlanner(path_rtn,
config)` replays a saved geometric path without searching: it times the path
for the run's request with the HYPR retimer, and the run validates the timed
reference against the requested clearance before the controller uses it. A
reference the run refuses is not flown. Timing a fixed path is deterministic
and draws no random numbers. In `examples/Earth_RPO_CubeSat_MPC.jl`,
`build_rpo_cubesat_mpc_demo` accepts `precomputed_plan` for such a replay. Its
chaser has a 5 kg dry mass plus 0.2 kg of propellant by default, and
`chaser_dry_mass_kg=4.8` gives a 5 kg wet mass instead. Its
`minimum_throttle_fraction=0.2` suppresses commands at or below 20% of the
rated thrust and holds the remaining average thrust; this does not model
individual valve pulses or establish the paper's pulse-width modulation.

The comparison helpers' reduced two-body evaluator holds fixed inertial body
axes, as the full example does, and checks swept clearance against the
buffered station point cloud at every 0.05 s plant substep. That is a sampled
check, not a collision test against watertight CAD geometry, and the
evaluator's fuel figures stay diagnostic until they agree with the full
simulator for the experiment being reported.

`RPOGuidanceModel`, the guidance model that braked and held before each
replan, was removed in SpaceAGORA 0.5.0. On the planner path,
`on_rejection=:install_fallback` installs a hold instead, and the chaser
counts as settled on it within 0.02 m and 0.002 m/s by default
(`hold_position_tolerance_m`, `hold_speed_tolerance_mps`). These tolerances
are simulation settings, not values given in the paper.
[Moving from RPOGuidanceModel](rpo_planner_pilot.md#Moving-from-RPOGuidanceModel)
lists the policy differences. To keep an obstacle in the geometry for a while
after it was last observed, set `observation_window_s` on
`RPOReplanningSpheres` or `RPOReplanningConfig`.
