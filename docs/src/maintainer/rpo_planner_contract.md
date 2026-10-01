# RPO planner contract development

The internal `SpaceAGORA.RPOPlannerInterfaces` module defines an initial
algorithm-independent request, timed reference, result and pure validation API.
It loads using Julia standard libraries without HYPR configuration. The module
is intentionally not exported as a stable package API while the integrated
planner pilot is being developed.

This change does not connect the contract to simulation callbacks, replace
`RPOGuidanceModel`, install an adapter, or make HYPR optional. Existing planner,
retiming and controller calculations continue through their existing route.

```@docs
SpaceAGORA.RPOPlannerInterfaces
```

## Data and ownership

`RPOPlanningRequest` carries stable spacecraft IDs, an epoch and simulation
time, truth state in target RTN coordinates, a position goal, copied geometry,
explicit constraints and reference timing. RTN velocities mean derivatives in
the rotating target frame, including its transport term. Geometry must be
static and aligned with that frame. The epoch and geometry representation are
opaque copied data; assembly is responsible for their meaning and provenance.

A request owns its geometry and epoch copies and stores state vectors as
immutable tuples. The lifecycle must retain its validation request and pass a
separate deep copy to the planner. A Julia struct does not make its contained
arrays immutable. `RPOReference` copies its arrays, and `RPOPlanningResult`
copies its reference and diagnostics. These boundaries prevent input-array
aliasing; they are not a sandbox for arbitrary plugin code. Integration must also
copy the result into lifecycle-owned storage before validation and install that
same validated copy. Later edits through a planner-held handle must not alter the
active reference. Mutation tests for both boundaries remain lifecycle acceptance
requirements.

References use uniformly spaced relative times beginning at zero, with 3 by N
position and velocity matrices. Their origin is simulation time. The existing
controller indexes columns using its control interval, so this contract rejects
nonuniform grids rather than resampling them. The sole initial terminal policy,
`:repeat_last_sample`, repeats both position and velocity. It does not imply a
stationary physical hold.

## Validation and explicit failures

`validate_rpo_result(request, result; clearance_at, time_s)` returns an
`RPOValidationResult` with acceptance, a reason symbol and available metrics.
It never installs a reference or advances a random stream. Checks include:

- request, spacecraft, frame and geometry-revision identity;
- finite values, dimensions, origin and uniform timing;
- position endpoints, declared and position-implied speed, and velocity-vector
  difference acceleration;
- reference lifetime, including every requested controller-preview time;
- clearance sampled along the reference polyline, including points between
  stored reference knots.

The trusted assembly supplies `clearance_at(point, geometry)` independently of
the planner. That query must return finite signed clearance with the declared
station/chaser approximation. Missing queries, nonfinite answers and excessive
validation work reject the candidate. Unexpected query exceptions propagate.
An algorithm's diagnostic path may be a control polygon; the validator never
uses it as an executable trajectory.

The default software tolerances are 1e-9 m for endpoints and 1e-10 s for time,
with clearance sampling at most every 0.05 m and at most 100,000 samples. These
are explicit validation settings, not newly accepted physical mission limits.
The time tolerance must remain smaller than half a sample interval. A caller
must prospectively select the settings appropriate to its fixture.

The speed limit checks both `max_declared_speed_mps` (velocity norms) and
`max_implied_speed_mps` (successive position displacements divided by the control
interval). `max_speed_mps` is their maximum. A reference that moves 100 m in
one second cannot pass a 1 m/s limit by declaring zero velocity. This does not
require the two velocity representations to be exact derivatives of each other.

Speed and acceleration comparisons allow only a bounded floating-point allowance:
`limit_roundoff_rtol` defaults to `128eps(Float64)` (about 2.84e-14), can be reduced
to zero for exact comparisons, and cannot exceed that ceiling. There is no
absolute allowance in SI units. An excess is rejected when it exceeds this
fraction of the enabled limit; measured values and the allowance are recorded.
This is a conservative software comparison rule, not an error bound for arbitrary
planner calculations or permission to loosen physical acceptance limits.

Common acceleration is the norm of successive velocity-vector differences per
control interval, including changes of direction. HYPR's tangential retiming
limit has a different meaning. The adapter must select and report any planning
headroom prospectively and still validate the returned reference against the
unchanged common limit. Increasing numerical tolerance to hide this distinction
is not permitted. The existing seed-742 fixture remains rejected at 0.1 m/s²;
its measured 0.10000004817916151 m/s² is beyond the roundoff allowance.

Only candidates terminated by completion or an iteration limit are eligible
by default. Wall-clock-budget candidates additionally require
`allow_time_budget_candidate=true` and all normal validation checks. This
permission does not make time-budgeted planning reproducible.

These are discrete reference checks, not continuous collision certification,
kinematic-consistency proof, actuator feasibility or physical tracking
acceptance. A reference can pass them and still require further validation in
its declared application.

## Extension methods

A planner can subtype `AbstractRPOPlanner` and implement `planner_capabilities`,
`initialize_planner`, `plan_rpo!` and optionally `retime_rpo!`. Planning methods
require an explicit `AbstractRNG`; there is no default global random stream.
The initial defaults report unsupported planning/retiming and consume no RNG.
Capabilities default to empty support. A retiming refusal does not silently
invoke full replanning.

Capability checks describe the algorithm only. They are not a claim that
simulation startup, isolation or checkpoint/restart has been implemented for
this interface. The later opt-in pilot must reject checkpoint/resume until its
complete run-owned state can be restored. Existing checkpoint behavior remains
unchanged in this contract-only change.

## Tests and next integration

`test/unit/gnc/rpo_planner_contract_tests.jl` loads only the contract source in a
standalone module and tests rejection cases and ownership. It can run with
`--project=@stdlib`. `rpo_planner_compatibility_tests.jl` captures repeated
seeded legacy/manuscript HYPR calls and both existing retiming policies using
the complete package, and validates each real reference against its physical
request limits. Two rounding-only cases pass; the tangential/total acceleration
mismatch remains an explicit rejection. Both are included in the unit driver and therefore in
the PR shard plan. Compatibility fingerprints are compared on the same Julia
and dependency environment; they are not portable golden values across
versions or platforms.

Next, add adapters preserving HYPR's returned effective settings, followed by
an opt-in lifecycle, supported external constructor and run-report accessor.
Their acceptance must include callback ordering, initialization, failure,
expiry, multiple controlled spacecraft and repeated/concurrent runs. Publish
top-level names only with the corresponding public API inventory, docstrings
and external example. Optional package extraction follows that pilot.
