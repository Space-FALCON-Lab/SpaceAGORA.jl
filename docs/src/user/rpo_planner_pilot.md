# Run an RPO case with interchangeable planners

The opt-in RPO pilot connects a planner to the existing simulation and LQ-MPC
controller. It supports truth observations, a passive target on a circular Earth
orbit, and a static station point cloud aligned with the target's rotating RTN
frame. The existing research examples remain available.

```julia
using SpaceAGORA
args = make_rpo_configuration(planner=DirectRPOPlanner(), seed=741)
solution = run_simulation(args; return_solution=true)
report = rpo_run_report(solution)
report.planners[1].records
report.commands[1].log
```

Use `HYPRRPOPlanner(RPOPSOConfig(...))` as the `planner` to choose HYPR. The separate
`examples/rpo_planner_env` project provides exact, deterministic settings and runs
both choices through this interface. Its README documents setup and offline
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
These are pilot settings, not a certification of tracking or actuator feasibility.

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
static scene; moving obstacles and replanning around a newly changed geometry are
outside this pilot.

Every failed initial, forced, background or retiming request stops this route with
`RPOPlanningError`, which carries copied diagnostics, spacecraft/request IDs, and
any original exception and backtrace. There is no implicit retain-old-plan or
emergency-controller fallback. An explicit `:retime` event requires the planner's
retiming capability; the two built-in adapters do not advertise it.

Immediately before control reads a reference, the lifecycle checks its lifetime
through the entire preview. Repeating the last reference column cannot bypass
expiry. `plan_validity_s` is a request lifetime, not the duration of the maneuver.
A short plan may repeat its final sample while still valid. A long plan must fit
inside its lifetime. Reference validity does not guarantee a stationary vehicle.

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

Checkpoint writing/resume and `isolate_state=false` are refused before propagation
or output creation. Faithful restart needs a future contract for planner/RNG,
controller, active-reference and trigger state. General restart and a core-only
installation with HYPR absent remain open work.

## Add your own planner or force

Subtype `AbstractRPOPlanner`. Implement `planner_capabilities`, optionally
`initialize_planner`, and `plan_rpo!(state, planner, request, rng)`. Return an
`RPOPlanningResult` containing an `RPOReference` for candidate status. Reference
times are uniform and relative to `origin_time_s`; vectors are SI target RTN,
including the rotating-frame transport term in velocities. The contract guide
lists validation rules. Add `retime_rpo!` only when you also declare retiming.

The separate-project smoke defines a tiny test planner using only these public
methods, including an explicit failed result. Its force example uses the existing
`AbstractForceTorqueModel` and `wrench` interface, passed through
`extra_effectors=(model,)`. No internal package module or source include is needed.
