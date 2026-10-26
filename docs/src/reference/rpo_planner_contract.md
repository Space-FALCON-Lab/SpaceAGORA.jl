# RPO planner contract

Use this reference when implementing a planner or interpreting an accepted
reference. The [RPO pilot guide](../user/rpo_planner_pilot.md) describes supported
simulation setup, and the [RPO Planner API](../generated/rpo_planner_api.md)
lists the public constructors and extension methods. Optional HYPR planning
uses the [versioned HYPR services](hypr_services.md).

## Data and ownership

`RPOPlanningRequest` carries stable spacecraft IDs, an epoch and simulation
time, truth state in target RTN coordinates, a position goal, copied geometry,
explicit constraints and reference timing. RTN velocities mean derivatives in
the rotating target frame, including its transport term. Geometry must be
static within a request and aligned with that frame. The epoch and geometry representation are
opaque copied data; assembly is responsible for their meaning and provenance.

A request owns its geometry and epoch copies and stores state vectors as
immutable tuples. The lifecycle must retain its validation request and pass a
separate deep copy to the planner. A Julia struct does not make its contained
arrays immutable. `RPOReference` copies its arrays, and `RPOPlanningResult`
copies its reference and diagnostics. These boundaries prevent input-array
aliasing; they are not a sandbox for arbitrary plugin code. Integration must also
copy the result into lifecycle-owned storage before validation and install that
same validated copy. Later edits through a planner-held handle must not alter the
active reference. The lifecycle tests cover mutation at both boundaries.

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

A planner can subtype `AbstractRPOPlanner` and implement `planner_capabilities`
and `plan_rpo!`, and optionally `initialize_planner` (the default returns
`nothing`, for stateless planners) and `retime_rpo!`. Planning methods
require an explicit `AbstractRNG`; there is no default global random stream.
The initial defaults report unsupported planning/retiming and consume no RNG.
Capabilities default to empty support. A retiming refusal does not silently
invoke full replanning.

## Geometry that changes during a run

A geometry provider, a subtype of `AbstractRPOGeometryProvider`, supplies the
geometry of each request. `rpo_geometry_revision(provider, revision, time_s)`
returns the revision at a time and `rpo_geometry_snapshot(provider, geometry,
time_s)` the geometry; both receive the configured base. Equal revisions must mean
equal snapshots, and providers draw no random numbers. The binding copies one
snapshot into each request, so the geometry is still static within a request and
validation checks the reference against the request's own copy and revision.

When the revision changes, the binding checks the installed reference's remaining
polyline, from its position at the current time to its end, against the new
snapshot. It uses the trusted `clearance_at` query, the request clearance and the
validation sample spacing and work cap; exceeding the cap counts as blocked. A
reference that still clears is kept, and the check is recorded. A blocked
reference, or a hold, is replaced by a `:replan` request. `RPOReplanningSpheres`
provides the spheres of the legacy replanning model.

Capability checks describe the algorithm only. Use the pilot's supported
lifecycle to install validated references. A checkpoint holds the binding's
run-owned state, including the algorithm's own state from `initialize_planner`;
a planner that keeps such state (anything but `nothing`) declares `restart` in
its capabilities to run with checkpoints, and the binding refuses it otherwise.

## The general algorithm contract

The RPO contract is the first family of SpaceAGORA's general algorithm contract,
the [GNC Algorithm Contract API](../generated/gnc_algorithm_contract_api.md).
`AbstractRPOPlanner` is an `AbstractGNCAlgorithm` of the `RPOReferenceFamily`,
and the general functions forward to the RPO names:

| General | RPO |
| --- | --- |
| `propose!(state, algorithm, request, rng)` | `plan_rpo!` |
| `propose_revision!(state, algorithm, request, active_product, rng)` | `retime_rpo!` |
| `algorithm_capabilities`, `initialize_algorithm`, `require_algorithm_support` | `planner_capabilities`, `initialize_planner`, `require_planner_support` |
| `validate_capabilities(family, ...)`, `validate_proposal(family, ...)` | `validate_rpo_capabilities`, `validate_rpo_result` |
| `is_current(family, ...)` | `rpo_reference_is_current` |

The family's rejection policy is `:fail` by default, and `RPOPlanningError` is a
`GNCAlgorithmError`. `RPOReferenceFamily(on_rejection=:install_fallback)` installs
the family's `fallback_product` instead: a hold of two samples at the chaser's
position with zero velocity, valid for the request's lifetime and carrying its
geometry revision. The hold is installed without validation, after the rejection
is recorded. A planner declares the default family; the binding configures the
policy through `bound_family`, which must return an instance of the planner's
family type. A planner that draws no random numbers declares
`consumes_rng(planner) = false`, as `DirectRPOPlanner` does; contract tests then
check that a proposal leaves the stream unchanged. The RPO names remain the
supported way to write an RPO planner.

The pilot lifecycle runs each request through the shared driver,
`run_algorithm_request!`. It builds and records the request, checks capabilities,
asks for a proposal with a separate copy, validates the copy against the retained
request, installs the validated reference into the plan buffer (`install!`) and
records the install. Records, reasons and reports are the same as before, with one
exception: an `RPOPlanningError` thrown by the planner's own code now fails the run
as `:planner_exception`, with that error as the cause, like any other planner
exception. The
planner's stream comes from `gnc_stream_words("RPO/v1", seed, chaser_id,
target_id)`, the same SHA-256 mapping as before.

## Built-in planners

`DirectRPOPlanner` builds a straight segment with `s(q)=3q²-2q³`, sampled
uniformly through the first grid time at or beyond its chosen duration. Its
analytic maxima are `1.5*distance/duration` for speed and
`6*distance/duration²` for acceleration; endpoints are at rest. The duration is
chosen from the reduced limits, rounded up to the request grid, and bounded by
`max_reference_samples` and the request lifetime before allocation. A
stationary request receives two identical samples. The baseline requires an
assembly-owned analytic `segment_clearance_at(start, goal, geometry)` query and
refuses missing queries, nonfinite clearance and blocked segments; it does not
search for an alternate route.

`HYPRRPOPlanner` copies its HYPR configuration, then maps request clearance and
sampling interval explicitly. Request limits with headroom cap configured
retiming limits. Conflicting configured minimum/initial speeds reject. A
requested acceleration limit requires the configured acceleration-limited
retimer; the adapter does not silently select another retiming algorithm.
`RPOPSOConfig()` leaves `retime_accel_limit_enable` off while the default
`make_rpo_configuration` constraints include an acceleration limit, so pair a
requested acceleration limit with `retime_accel_limit_enable=true`, as
`examples/rpo_planner_env/smoke.jl` does, or omit the limit.
`rrt_on_replan=true` explicitly enables the existing RRT warm start for a
`:replan` request.

`RetainedPathRPOPlanner(path_rtn, config)` replays a saved geometric path, such
as a published HYPR plan, instead of searching. For a request whose chaser
position and goal are within `endpoint_tolerance_m` (default 1e-6 m) of the
path's first and last points, it sets those points to the request's and times the
path with the HYPR retimer, under the request mapped onto `config` exactly as
`HYPRRPOPlanner` maps it for a new plan. Otherwise it refuses
(`:infeasible`, with `:start_not_on_path` or `:goal_not_on_path`). It draws no
random numbers, declares neither retiming nor restart, and needs the HYPR
companion. In a run, the timed reference passes the lifecycle's validation like
any other; a saved path is never installed unchecked. Replaying a HYPR plan's
path and effective configuration for the same request gives the same reference
as the plan.

`HYPRRPOPlanner` also retimes. Its `retime_rpo!` keeps the active reference's
validated polyline from the point nearest the chaser to its end, starts it at the
chaser's position (the route rule of the legacy replanning model) and times it as a
polyline with the HYPR retimer, under the request limits mapped as for a new plan.
It searches no new path and draws no random numbers. It refuses an active
reference that does not end at the request's goal (`:infeasible/:goal_not_on_reference`),
and applies the same lifetime, work and validator checks as `plan_rpo!`. With
`hypr_mode=:legacy` and the acceleration-limited retimer, a retime from off the
reference can turn the corner where it joins the route faster than the common
acceleration limit; validation then refuses it (`:failed/:reference_rejected`).

Both planners report `:infeasible/:insufficient_reference_lifetime` when the
absolute reference end exceeds the request lifetime, without adding a time
tolerance to this planning budget; equality is allowed. The HYPR adapter then
applies the unchanged validator and returns `:failed/:reference_rejected` with
the validation record for other output failures. Optimizer exceptions
propagate. A `:candidate` is still subject to lifecycle-owned validation before
installation, and every non-candidate outcome reaches the run as an
`RPOPlanningError` carrying these symbols.

Results report their outcome in the `status` and `termination` fields and
can carry diagnostics. A direct-baseline candidate records its planning budget, chosen
duration and analytic maxima. A HYPR adapter candidate distinguishes the
requested, mapped and optimizer-returned configurations and includes the
planning budget, path representation and deterministic cost and its history.
Refused results, whether refused before or after the optimizer runs, may carry
only the planning budget or no diagnostics. Because the adapter plans within the reduced limits described
below, its reference timing can differ from a legacy HYPR call that uses the
unreduced limits.

### Comparison planners

The seven modes of the RPO planner comparison are planners too:
`HYPRComparisonPlanner` (`:hypr`), `UnrefinedPSOPlanner` (`:pso_unrefined`),
`RRTConnectPlanner` (`:rrt_connect`), `RRTConnectBezierPlanner`
(`:rrt_connect_bezier`), `RRTStarPlanner` (`:rrt_star`), `CHOMPPlanner`
(`:chomp`) and `STOMPPlanner` (`:stomp`). Each holds a copy of its
`RPOPSOConfig`, which sets the cost and the retiming, and of its search or
optimizer settings. Its `plan_rpo!` runs the comparison's kernels under the
request's clearance, reference step, lifetime and limits, mapped with headroom as
for `HYPRRPOPlanner`; the comparison's fixed 0.5 m safety distance does not
apply. It plans without a wall-clock limit and times the path with the
configuration the kernel returned. A search that finds no path reports
`:failed/:path_not_found` and never offers its direct-segment diagnostic as a
candidate. Under an acceleration limit, the acceleration-limited retimer cannot
hold the limit at the corners of the polyline paths of `RRTConnectPlanner` and
`RRTStarPlanner`, so validation refuses them (`:failed/:reference_rejected`); the
curve-based modes can meet it. The comparison planners declare truth-state,
target-RTN planning without retiming or restart, draw from the planner's stream,
and need the optional HYPR companion, as the comparison does.

The comparison harness, `RPOPlannerComparisonConfig`, takes these objects in
`planners`; by default it compares all seven, built from its own settings
fields. Selecting a mode by name (`:rrt_star`, `"rrt-star"` and the other
aliases) was removed in SpaceAGORA 0.5.0: pass the object. The names remain the
identifiers and labels of the comparison's report rows.

## Adapter planning headroom

`RPOPlanningHeadroom(fraction=0.01)` prospectively reserves one percent of each
enabled request speed/acceleration limit for the new adapters. It changes the
planning target, never the physical request or the shared validator. The fraction
must be strictly between zero and one and is recorded with the result. A configured
HYPR cap that is already lower is preserved.

The adapters estimate output resolution using `32eps(Float64)*position_scale/dt`
for implied speed and `32eps(Float64)*velocity_scale/dt` for acceleration. Each
estimate must fit within one quarter of its reserved physical-limit budget.
They check request scales before planning and actual output scales before
returning a reference. Inadequate resolution yields
`:insufficient_reference_precision`; it never increases validation tolerance.
These conservative screening estimates are not a proof for arbitrary floating
point algorithms. Large coordinates or very small time steps may require a
better representation or different sampling, not permission to exceed a limit.

The fixed `128eps(Float64)` validator allowance does not cover arbitrary
finite-difference storage error. The earlier rounding-only fixtures pass at their
fixture scale. The adapter reserve covers a wider tested range of distances and
steps while retaining strict validation, including exact comparisons when the
caller selects zero allowance.

One percent is a declared pilot policy, not a guarantee of total-vector
acceleration. HYPR's tangential acceleration and curvature constraints can still
produce a larger vector difference. The adapter rejects a reference that fails
the unchanged physical checks; it does not retry with relaxed tolerances or choose
a larger reserve after seeing a failure. Replacing this policy requires a new
prospective setting and its own validation.

For a smooth path with orthogonal tangential and normal acceleration, saturating
its tangential part at `(1-f)*a_max` leaves at most
`a_max*sqrt(1-(1-f)^2)` for turning. At the default `f=0.01`, this is about
14.1 percent of the physical acceleration limit. This geometric estimate is
not a guarantee for sampled velocity differences; the unchanged discrete
validator still decides whether each candidate is admissible.
