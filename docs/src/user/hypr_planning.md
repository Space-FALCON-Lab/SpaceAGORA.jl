# Motion planning with HYPR

!!! note "Optional package"
    Requires the optional `HYPR` package; see
    [Install and try the planner](#Install-and-try-the-planner).

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
Load `SpaceAGORA` and `HYPR` directly. Loading SpaceAGORA alone
does not activate HYPR.

From a SpaceAGORA Git clone, prepare a separate environment and run the public
planner example:

```sh
julia --startup-file=no scripts/optional_packages.jl setup hypr /path/to/hypr-env
JULIA_LOAD_PATH='@:@stdlib' julia --startup-file=no --project=/path/to/hypr-env \
  examples/rpo_planner_env/smoke.jl
```

The setup helper uses the verified submodule at `packages/robotics/HYPR.jl`.
Its recorded revision agrees with `packages/HYPRSource.toml`.
The robotics folder is a source organization choice; HYPR supports spacecraft
RPO as well as robot arms. The [installation guide](installation_environment.md#Optional-package-commands)
describes source checks, and the [RPO guide](rpo_planner_pilot.md#Installing-optional-HYPR)
retains the separate `scripts/setup_hypr.jl` route. Keep the HYPR pin and
SpaceAGORA revision in the run record. The short example executes
both the direct planner and HYPR through the public interface, checking that
references reach the controller and produce finite, bounded outputs. It covers
only a two-second portion of a maneuver. Its small swarm, simple geometry and
short controller preview are installation checks, not the paper experiment.

For a new application, load `SpaceAGORA` and `HYPR`, construct an
`HYPRRPOPlanner` with an explicit `RPOPSOConfig`, and pass it as `planner` to
`make_rpo_configuration`. Use `hypr_mode=:manuscript` when investigating the
manuscript formulation; the configuration default is `:legacy`. The
[public RPO guide](rpo_planner_pilot.md) explains the complete call sequence,
reference validation, seed mapping and the information returned by
`rpo_run_report`.

## Relating a run to the paper

The selected final author-shared iSpaRo paper describes a 500-case static
comparison and 100 baseline paths under eight additional scenarios, or 800
dynamic cases per planner. The earlier September 23 draft described 900. Its static comparison uses MersenneTwister seeds `740 + i`, a
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
[planner contract](../reference/rpo_planner_contract.md) and
[HYPR service reference](../reference/hypr_services.md) describe those boundaries
for developers extending the integration.


## Final iSpaRo reconstruction candidate

The private SpaceAGORARPOStudies project owns
`examples/Earth_RPO_HyPR_FinalISpaRo_Reconstruction.jl`, which records a newly
generated experiment through the full SpaceAGORA simulator. It uses the final paper's
0.5 m clearance, 10 Hz control, 60-step horizon and controller weights, with
5 kg interpreted as initial wet mass and 0.2 kg of propellant. Exact original
case coordinates, geometry draws and the paper-producing source revision have
not been recovered. This example is a reconstruction, not an original-run replay.

The following reconstruction commands are for lab members with access to that
repository; outside collaborators need an invitation. Core users can use the
generic planner examples described above.

The example defaults to one HyPR case. Select other methods explicitly with
`HYPR_PLANNERS`; comparisons are optional, not an installation requirement.
The paper implementation walkthrough (SpaceAGORARPOStudies) gives runnable commands,
algorithm mapping and the retained dynamic failure. Write evidence outside the repository and any
cloud-synced vault, for example:

```sh
cd /path/to/SpaceAGORARPOStudies
JULIA_LOAD_PATH=@:@stdlib HYPR_CASES=1 HYPR_OUTPUT_ROOT="$HOME/spaceagora-hypr-final-isparo" \
  JULIA_NUM_THREADS=2 julia --startup-file=no --compiled-modules=existing \
  --project=/path/to/rpo-study-env \
  examples/Earth_RPO_HyPR_FinalISpaRo_Reconstruction.jl
```

Prepare `/path/to/rpo-study-env` using the research project's setup instructions.
It pins SpaceAGORA to `bfa32f6a44c5f242714b8cf182a451e863119a05` and HYPR to
`c434404dd3e4ac3ba351b2a5e0bb77d9be95fb94`; keep those dependency checkouts
at the recorded revisions. The command runs from the research checkout, not
from the core examples directory.
The example records cases, each planner's path/reference, successful and failed
outcomes, and actual delivered-force fuel consumption. Failed cases remain in
success/runtime accounting and must be excluded from fuel means. Comparisons
against HyPR must use the intersection of successful case IDs. Hardware and
thread count are part of the experiment; this local run does not reproduce the
paper's 16-thread AMD planning times.

Each planner's plan is flown on the planner contract: `RetainedPathRPOPlanner`
times it for the run's request, and the run validates it against the 0.5 m clearance
before the controller uses it. A plan the run refuses is not flown. Its
evaluation records the refusal (`refusal`) and fails, where the evaluation before
this change flew the plan and counted clearance violations. Replaying a plan does
not change its timing: the published case-4 plan replays to its recorded
reference exactly.

The full demo holds a fixed inertial attitude and samples clearance along
0.05 s swept chords against the buffered point cloud. This does not certify
watertight CAD collision geometry. The 20% activation option suppresses small
average-throttle commands; it does not establish resolved PWM valve behavior.
The reduced two-body evaluator remains diagnostic until its fuel agreement with
the full simulator is established for the experiment being reported.

Through SpaceAGORA 0.4.0, the HYPR guidance model `RPOGuidanceModel` braked and
held before replanning. It consumed measured planner time in simulation time,
checked the current map before release, and left failed replacements in hold. Its
position/speed release tolerances (0.02 m and 0.002 m/s) were explicit
reconstruction settings recovered from the later source, not numerical tolerances
specified by the final paper. SpaceAGORA 0.5.0 removed the model: on the planner
path a rejected request installs a hold (`on_rejection=:install_fallback`, with
the same default tolerances), with the policy differences listed in
[Moving from RPOGuidanceModel](rpo_planner_pilot.md#Moving-from-RPOGuidanceModel).
Retiming retains the remaining geometric route.
`RPOReplanningConfig.observation_window_s` optionally retains detections after
their last observation; a reconstruction can set it to 0.02 times the original
reference duration.

The replanning study `Earth_RPO_HyPR_FinalISpaRo_Replanning.jl`, removed with the
model, consumed a saved case record via `HYPR_BASELINE_RECORD`. It exercised a
baseline plus the eight reconstructed disturbance families in the full simulator,
and kept a separate obstructed-goal negative control. It used the callback
ordering of its release: the scenario input arrived before control; control ran
before guidance. Fixed endpoints that already violated clearance were rejected
before search, and the vehicle remained in hold. A failed hold is not a claim that
an already unsafe position has become safe. The original obstacle placements and
disturbance draws are missing, so these generated cases cannot be attributed to
the final table. The “Run the HyPR paper algorithm” guide in
SpaceAGORARPOStudies explains how to run the study with SpaceAGORA 0.4.0.

Remaining-route resampling subdivides each original segment so that a coarse
sampling interval cannot remove a corner and create a new shortcut. Replacement
release checks both the earlier planned route and the final retimed reference.
These checks preserve the 0.5 m limit; passing them does not qualify a complete
maneuver or its runtime. The reconstruction records unsuccessful and interrupted
runs separately from successful tracking.

Retiming follows Algorithm 1 lines 14-15: project onto the retained geometric
route, connect the measured position to the next geometric vertex, and time the
entire connected path with the existing clearance/curvature speed law. Collinear
time samples are removed only when all removed points lie within the existing
1e-10 m duplicate tolerance of a monotone retained segment. Genuine bends,
reversals and coincident-endpoint excursions remain vertices. This prevents a
small tracking offset from creating a sharp turn at the next time-grid sample.
Newly created plans retain the geometric sampling polyline produced from the
planner's actual curve configuration before time integration changes sample
spacing. Replacement release and subsequent retimes preserve that nominal
polyline separately from measured-state connectors. Explicit polyline corners
are retained even below the usual sampling interval. Legacy plans without this
metadata use their delivered reference as a polyline; Bezier control points are
never inferred from timed samples. Retained geometry must be finite 3-by-N
points ending at the reference goal, and release verifies the live goal again.
The connector consumes traversal time; the first position is never overwritten
after timing. The measured initial velocity is retained as boundary evidence,
not clipped to manufacture instantaneous compliance with a newly tightened cap.
Previous connectors are excluded from later geometric-progress calculations.
The future-time distance in line 7 remains the tracking-error criterion. Unsafe
stitched references trigger braking/replanning instead of being installed.

Replacement release checks the effective live goal, including the model goal
when no desired-goal override is set. A goal changed while planning, or an
override removed during that interval, cannot silently release the stale plan.

Dynamic appearance may create a state from which the clearance requirement is
physically infeasible under the available thrust. Such a case remains a failed
maneuver in unconditional success accounting. Braking and replanning do not
imply that every obstacle appearance admits a safe response. Operational
appearance assumptions and any stronger safety policy require their own
scientific disposition; neither the safety distance nor event timing is relaxed
by this reconstruction.

An expired reference, including a multipoint reference, is not repeatedly retimed solely
because its clearance lies in the retiming band. This preserves its completion
clock. Changed goals, unsafe clearance and excessive tracking error are still
checked before this no-traversal case. Expiry uses the absolute reference deadline
to preserve exact callback-boundary behavior without adding an arrival tolerance.
