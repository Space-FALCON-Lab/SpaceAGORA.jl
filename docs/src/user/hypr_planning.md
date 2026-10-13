# Motion planning with HYPR

HYPR plans a route through a cluttered environment and turns that route into a
reference that a spacecraft controller can follow. In SpaceAGORA, its rendezvous
and proximity operations implementation is intended for motion near a target
spacecraft or station. The planner chooses the route and its timing; the
simulation then reveals how well the controlled vehicle follows it and how much
propellant the maneuver requires. These are separate questions. A low planning
cost is useful evidence about a candidate route, but it is not a measurement of
the fuel consumed during tracking.

The method is described by Jakob Gorisek-Gazze and Giusy Falcone in *HyPR:
Safety- and Fuel-Aware Motion Planning and Replanning for Autonomous SmallSat
Servicing*, accepted to the IEEE International Conference on Space Robotics
(iSpaRo), 2026. The [Space-FALCON Lab publication record](https://www.spacefalconlab.com/publications/)
provides the citation. The discussion here distinguishes that study from the
current software interfaces and from the smaller examples used to check an
installation.

## From a route to a controlled maneuver

The translational method first searches for a feasible route with RRT-Connect.
That route supplies an initial guess for a particle swarm, whose particles
represent candidate Bézier paths. The swarm improves the geometric path using
the configured obstacle and fuel-related objective. In manuscript mode, the
search effort and coefficient schedules adapt to the warm-start problem. The
returned path is the best evaluated candidate; the final particle positions
should not be interpreted as an additional solution.

Retiming assigns speeds and times to the selected geometry. Clearance, reaction
time and acceleration limits constrain how quickly the vehicle should traverse
each portion of the route. Acceleration-based speed caps are distinct from the
optional forward and backward passes that bound tangential acceleration;
those passes are disabled by default. SpaceAGORA's LQ-MPC controller then tracks the
reference using its own prediction horizon, weights and actuator model. The
planner's fuel proxy and the tracked maneuver's propellant consumption must
therefore be reported separately. Changing the controller can change the latter
even when the geometric path is unchanged.

Planning coordinates are expressed in the target's rotating radial, transverse
and normal frame, with distances in metres and time in seconds. Relative
velocity in that frame includes the rotation transport term. The station
geometry and the chaser's extent jointly determine the clearance calculation;
a distance from the chaser centre to a station sample is not, by itself, the
available collision margin. Keep the geometry, sampling resolution and safety
distance with the experiment settings when comparing results.

## Install and try the planner

The search implementation lives in the separate
[HYPR.jl repository](https://github.com/Space-FALCON-Lab/HYPR.jl). SpaceAGORA
provides the simulation, shared geometry and the versioned integration services.
The small `SpaceAGORAHYPR` package loads the two together. Loading SpaceAGORA alone
does not activate HYPR.

From a SpaceAGORA checkout, prepare a separate environment and run the public
planner example:

```sh
julia --startup-file=no scripts/setup_hypr.jl /path/to/hypr-env
JULIA_LOAD_PATH='@:@stdlib' julia --startup-file=no --project=/path/to/hypr-env \
  examples/rpo_planner_env/smoke.jl
```

The setup helper installs the exact HYPR revision in
[`packages/SpaceAGORAHYPR/HYPRSource.toml`](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/main/packages/SpaceAGORAHYPR/HYPRSource.toml). Keep that
pin and the SpaceAGORA revision in the run record. The short example executes
both the direct planner and HYPR through the public interface, checking that
references reach the controller and produce finite, bounded outputs. It covers
only a two-second portion of a maneuver. Its small swarm, simple geometry and
short controller preview are installation checks, not the paper experiment.

For a new application, load `SpaceAGORA` and `SpaceAGORAHYPR`, construct an
`HYPRRPOPlanner` with an explicit `RPOPSOConfig`, and pass it as `planner` to
`make_rpo_configuration`. Use `hypr_mode=:manuscript` when investigating the
manuscript formulation; the configuration default is `:legacy`. The
[public RPO guide](rpo_planner_pilot.md) explains the complete call sequence,
reference validation, seed mapping and the information returned by
`rpo_run_report`.

## Relating a run to the paper

The September 23 manuscript describes a 500-case static comparison and 900
dynamic cases. Its static comparison uses MersenneTwister seeds `740 + i`, a
0.5 m safety distance, a 0.1 s controller update and a 60-step MPC horizon. The
controller weights are `Q = diag(20,20,20,2,2,2)`, `R = 0.1I` and `Qf = 20Q`.
The study gives competing planners a per-case time budget based on HYPR's
measured planning time. Consequently, changing the processor or thread count
can change both the reported runtime and the amount of optimization performed
by a competing planner.

The repository's static and dynamic research entry points are
`examples/Earth_RPO_CubeSat_MPC_PlannerComparison.jl` and
`examples/Earth_RPO_CubeSat_MPC_Replanning.jl`. Their defaults are useful
configurations, but they are not an immutable transcription of Table I. In
particular, the current comparison defaults use a different MPC horizon and
weights, inherit the legacy planner mode, and leave acceleration-limited
retiming disabled. The batch driver also derives planner seeds from a master
random stream, rather than passing `740 + i` directly to every planner. Retain
the actual endpoints, sampled geometry and random states, along with every
effective setting; a seed number alone does not identify an experiment.

There is also a difference in the tracking model. The current comparison helper
propagates a discrete linear Hill-Clohessy-Wiltshire model and accounts for
continuous commanded impulse. The manuscript describes nonlinear two-body
propagation with a 20% pulse-width-modulation activation threshold. Matching
the weights in Table I does not reconcile that difference. A run of this helper
is useful for checking planning and tracking, but cannot by itself establish
agreement with the paper's fuel statistics.

The paper computes static fuel statistics over successful runs. Its matched
comparison uses the same successful cases for each competing planner and HYPR;
comparing unrelated successful subsets changes the question being answered.
Dynamic results use a different baseline, the HYPR static-map fuel consumption
for that experiment. Neither quantity should be substituted for the other.
Keep unsuccessful cases in the success-rate denominator and preserve their
failure reasons.

The package's manuscript tests check selected equations, seeded search behavior
and integration contracts. They do not establish agreement with all of the
paper's aggregate results. A reproduction also needs the paper's experiment
inputs, scenario definitions and complete outcome records. Until that comparison
is complete, cite the numerical results as manuscript-reported results rather
than as measurements of the current release.

## Replanning and the supported interface

The paper and the research replanning example distinguish keeping a reference,
retiming it and searching for a new route after a meaningful change in clearance,
tracking error or goal. The public planner pilot currently has a narrower
contract: it supports a static scene and explicit planning requests, and it
records why each request was accepted or rejected. Moving obstacles and the
paper's complete dynamic experiment are not established by that pilot.

For each run, inspect the installed reference, its sampled clearance checks and
the executed trajectory. A validated reference does not guarantee collision-free
tracking under arbitrary disturbances. The controller, actuator limits and
failure policy remain part of the application. The
[planner contract](../maintainer/rpo_planner_contract.md) and
[HYPR service reference](../maintainer/hypr_services.md) describe those boundaries
for developers extending the integration.
