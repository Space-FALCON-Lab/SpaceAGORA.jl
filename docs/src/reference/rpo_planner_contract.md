# RPO planner contract

This page is the reference for writing an RPO planner and for reading the
references it produces. The [RPO pilot guide](../user/rpo_planner_pilot.md)
describes the supported simulation setup, and the
[RPO Planner API](../generated/rpo_planner_api.md) lists the public constructors
and extension methods. Optional HYPR planning uses the
[versioned HYPR services](hypr_services.md).

## Data and ownership

An `RPOPlanningRequest` carries stable spacecraft IDs, an epoch and simulation
time, the truth state in target RTN coordinates, a position goal, copied
geometry, explicit constraints and reference timing. RTN velocities are
derivatives in the rotating target frame, including its transport term. Geometry
must be static within a request and aligned with that frame. The epoch and the
geometry are opaque copied data; the code that assembles the request is
responsible for their meaning and provenance.

A request owns its copies of the geometry and the epoch and stores state vectors
as immutable tuples. The lifecycle keeps its own request for validation and
passes a separate deep copy to the planner, because a Julia struct does not make
the arrays it contains immutable. `RPOReference` copies its arrays, and
`RPOPlanningResult` copies its reference and diagnostics. These copies prevent
aliasing of input arrays; they are not a sandbox for arbitrary plugin code. The
integration also copies the result into lifecycle-owned storage before
validation, and installs that same validated copy, so later edits through a
handle the planner kept cannot alter the active reference. The lifecycle tests
cover mutation at both boundaries.

A reference has uniformly spaced relative times beginning at zero, with 3-by-N
position and velocity matrices; its time origin is the simulation time. The
controller indexes the columns by its control interval, so this contract rejects
nonuniform grids rather than resampling them. The only terminal policy so far,
`:repeat_last_sample`, repeats both the last position and the last velocity. It
does not imply a stationary physical hold.

## Validation and explicit failures

`validate_rpo_result(request, result; clearance_at, time_s)` returns an
`RPOValidationResult` with the acceptance, a reason symbol and the metrics it
measured. It never installs a reference or advances a random stream. It checks:

- the identity of the request, the spacecraft, the frame and the geometry
  revision;
- finite values, dimensions, origin and uniform timing;
- the position endpoints, the declared and the position-implied speed, and the
  acceleration from velocity-vector differences;
- the reference lifetime, including every requested controller-preview time;
- clearance, sampled along the reference polyline, including points between the
  stored reference knots.

The trusted assembly code, not the planner, supplies the query
`clearance_at(point, geometry)`. It must return a finite signed clearance with
the declared station/chaser approximation. A missing query, a nonfinite answer or
validation work beyond the cap rejects the candidate; an unexpected exception from
the query propagates. An algorithm's diagnostic path may be a control polygon, but
the validator never uses it as an executable trajectory.

The default software tolerances (`RPOValidationSettings`) are 1e-9 m for
endpoints and 1e-10 s for time, with clearance sampled at most every 0.05 m and
at most 100,000 samples. These are validation settings, not physical mission
limits. The time tolerance must stay below half a sample interval. A caller
chooses the settings for its case before it plans.

The speed limit applies to both `max_declared_speed_mps` (the norms of the
velocity samples) and `max_implied_speed_mps` (successive position displacements
divided by the control interval); `max_speed_mps` is the larger of the two. A
reference that moves 100 m in one second cannot pass a 1 m/s limit by declaring
zero velocity. The two velocity representations need not be exact derivatives of
each other.

Speed and acceleration comparisons allow only a bounded floating-point excess.
`limit_roundoff_rtol` defaults to `128eps(Float64)` (about 2.84e-14), can be
reduced to zero for exact comparisons, and cannot be set higher. There is no
absolute allowance in SI units. A value is rejected when it exceeds the enabled
limit by more than this fraction of the limit, and the measured values and the
allowance are recorded. This is a conservative rule for comparing numbers, not an
error bound for arbitrary planner calculations or a way to loosen physical
acceptance limits.

The common acceleration is the norm of successive velocity-vector differences
per control interval, so it includes changes of direction. HYPR's tangential
retiming limit means something different. The adapter must choose and report any
planning headroom before it plans, and still validate the returned reference
against the unchanged common limit. Raising the numerical tolerance to hide the
difference is not permitted. At a 0.1 m/s² limit, for example, the contract tests
accept an acceleration `8eps` above the limit and reject 0.10000004817916151 m/s².

By default, only candidates that ended by completion or at an iteration limit
are eligible. A candidate that ended at its wall-clock budget also needs
`allow_time_budget_candidate=true`, and must still pass all the normal checks. This
permission does not make time-budgeted planning reproducible.

These are discrete checks of a reference. They are not continuous collision
certification, a proof of kinematic consistency, actuator feasibility or physical
tracking acceptance, and a reference that passes them may still need further
validation in its application.

## Extension methods

A planner subtypes `AbstractRPOPlanner` and implements `planner_capabilities` and
`plan_rpo!`, and optionally `initialize_planner` (the default returns `nothing`,
for stateless planners) and `retime_rpo!`. Planning methods require an explicit
`AbstractRNG`; there is no default global random stream. The defaults report
planning and retiming as unsupported and consume no random numbers, and the
default capabilities support nothing. A retiming refusal never falls back to full
replanning on its own.

## Geometry that changes during a run

A geometry provider, a subtype of `AbstractRPOGeometryProvider`, supplies the
geometry of each request. `rpo_geometry_revision(provider, revision, time_s)`
returns the revision at a time and `rpo_geometry_snapshot(provider, geometry,
time_s)` the geometry; both receive the configured base. Equal revisions must mean
equal snapshots, and providers draw no random numbers. The binding copies one
snapshot into each request, so the geometry is still static within a request, and
validation checks the reference against the request's own copy and revision.

When the revision changes, the binding checks the remaining polyline of the
installed reference, from its position at the current time to its end, against
the new snapshot. It uses the trusted `clearance_at` query, the request's
clearance and the validation sample spacing and work cap; exceeding the cap counts
as blocked. A reference that still clears is kept, and the check is recorded. A
blocked reference, or a hold, is replaced through a `:replan` request.
`RPOReplanningSpheres` provides the spheres of the older replanning model.

Capability checks describe the algorithm only; use the pilot's supported
lifecycle to install validated references. A checkpoint holds the binding's
run-owned state, including the algorithm's own state from `initialize_planner`.
A planner that keeps such state (anything but `nothing`) must declare `restart`
in its capabilities to run with checkpoints, and the binding refuses it
otherwise.

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
is recorded. A planner declares the default family; the binding sets the policy
through `bound_family`, which must return an instance of the planner's family
type. A planner that draws no random numbers declares
`consumes_rng(planner) = false`, as `DirectRPOPlanner` does, and the contract
tests then check that a proposal leaves the stream unchanged. The RPO names remain
the supported way to write an RPO planner.

The pilot lifecycle runs each request through the shared driver,
`run_algorithm_request!`. The driver builds and records the request, checks
capabilities, asks for a proposal with a separate copy, validates the copy
against the retained request, installs the validated reference into the plan
buffer (`install!`) and records the installation. An `RPOPlanningError` thrown by
the planner's own code fails the run as `:planner_exception`, with that error as
the cause, like any other planner exception. The planner's random stream is seeded
from `gnc_stream_words("RPO/v1", seed, chaser_id, target_id)`, a SHA-256 hash of
the seed and the two spacecraft IDs.

## Built-in planners

`DirectRPOPlanner` builds a straight segment with `s(q)=3q²-2q³`, sampled
uniformly up to the first grid time at or beyond its chosen duration. Its
analytic maxima are `1.5*distance/duration` for speed and
`6*distance/duration²` for acceleration, and both endpoints are at rest. The
duration is chosen from the reduced limits, rounded up to the request grid, and
checked against `max_reference_samples` and the request lifetime before anything
is allocated. A stationary request receives two identical samples. The planner
needs an analytic `segment_clearance_at(start, goal, geometry)` query, which the
assembly supplies, and refuses a missing query, a nonfinite clearance and a
blocked segment; it does not search for another route.

`HYPRRPOPlanner` copies its HYPR configuration, then maps the request's clearance
and sampling interval onto it explicitly. Request limits, less the headroom, cap
the configured retiming limits. Configured minimum or initial speeds that exceed
the planning limit are refused. A requested acceleration limit needs the
configured acceleration-limited retimer; the adapter does not silently choose
another retiming algorithm. `RPOPSOConfig()` leaves `retime_accel_limit_enable`
off, while the default `make_rpo_configuration` constraints include an
acceleration limit, so either pair a requested acceleration limit with
`retime_accel_limit_enable=true`, as `examples/rpo_planner_env/smoke.jl` does, or
leave the limit out. `rrt_on_replan=true` enables the RRT warm start for a
`:replan` request.

`RetainedPathRPOPlanner(path_rtn, config)` replays a saved geometric path, such as
a published HYPR plan, instead of searching. For a request whose chaser position
and goal are within `endpoint_tolerance_m` (default 1e-6 m) of the path's first
and last points, it moves those points onto the request's and times the path with
the HYPR retimer, under the request mapped onto `config` exactly as
`HYPRRPOPlanner` maps it for a new plan. Otherwise it refuses (`:infeasible`, with
`:start_not_on_path` or `:goal_not_on_path`). It draws no random numbers, declares
neither retiming nor restart, and needs the HYPR companion. In a run, the timed
reference passes the lifecycle's validation like any other, so a saved path is
never installed unchecked. Replaying a HYPR plan's path and effective
configuration for the same request gives the same reference as the plan.

`HYPRRPOPlanner` also retimes. Its `retime_rpo!` keeps the active reference's
validated polyline from the point nearest the chaser to its end, starts it at the
chaser's position (the route rule of the older replanning model), and times it as
a polyline with the HYPR retimer, under the request limits mapped as for a new
plan. It searches no new path and draws no random numbers. It refuses an active
reference that does not end at the request's goal
(`:infeasible/:goal_not_on_reference`), and applies the same lifetime, work and
validator checks as `plan_rpo!`. With `hypr_mode=:legacy` and the
acceleration-limited retimer, a retime from off the reference can turn the corner
where it joins the route faster than the common acceleration limit allows;
validation then refuses it (`:failed/:reference_rejected`).

Both planners report `:infeasible/:insufficient_reference_lifetime` when the end
of the reference, in absolute time, is later than the request lifetime; equality
is allowed, and no time tolerance is added to this planning budget. For other
output failures the HYPR adapter applies the unchanged validator and returns
`:failed/:reference_rejected` with the validation record. Optimizer exceptions
propagate. A `:candidate` still goes through lifecycle-owned validation before it
is installed, and every outcome other than a candidate reaches the run as an
`RPOPlanningError` that carries these symbols.

Results report their outcome in the `status` and `termination` fields and can
carry diagnostics. A direct-planner candidate records its planning budget, the
chosen duration and the analytic maxima. A HYPR adapter candidate distinguishes
the requested, mapped and optimizer-returned configurations, and includes the
planning budget, the path representation, and the deterministic cost and its
history. A refused result, whether refused before or after the optimizer runs,
may carry only the planning budget or no diagnostics. Because the adapter plans
within the reduced limits described below, its reference timing can differ from
that of a HYPR call outside the adapter, which uses the unreduced limits.

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
`planners`; by default it compares all seven, built from its own settings fields.
Selecting a mode by name (`:rrt_star`, `"rrt-star"` and the other aliases) was
removed in SpaceAGORA 0.5.0: pass the object. The names remain the identifiers
and labels of the comparison's report rows.

## Adapter planning headroom

`RPOPlanningHeadroom(fraction=0.01)` reserves, before planning, one percent of
each enabled speed and acceleration limit of the request for the adapters. It
changes the planning target, never the physical request or the shared validator.
The fraction must be strictly between zero and one, and it is recorded with the
result. A configured HYPR cap that is already lower is kept.

The adapters estimate the output resolution as `32eps(Float64)*position_scale/dt`
for the implied speed and `32eps(Float64)*velocity_scale/dt` for the
acceleration. Each estimate must fit within one quarter of its reserved budget.
They check the request's scales before planning and the actual output's scales
before returning a reference. Inadequate resolution yields
`:insufficient_reference_precision`; it never increases the validation tolerance.
These conservative screening estimates are not a proof for arbitrary
floating-point algorithms. Large coordinates or very small time steps may need a
better representation or different sampling, never permission to exceed a limit.

The fixed `128eps(Float64)` validator allowance does not cover arbitrary
finite-difference storage error; the adapter reserve does, over the range of
distances and steps the tests cover, while validation stays strict, exact if the
caller sets the allowance to zero.

One percent is a declared policy of the pilot, not a guarantee on the total
acceleration vector. HYPR's tangential acceleration and curvature constraints can
still produce a larger vector difference. The adapter rejects a reference that
fails the unchanged physical checks; it does not retry with relaxed tolerances or
choose a larger reserve after seeing a failure. Replacing this policy needs a new
setting chosen in advance and its own validation.

For a smooth path whose tangential and normal accelerations are orthogonal,
saturating the tangential part at `(1-f)*a_max` leaves at most
`a_max*sqrt(1-(1-f)^2)` for turning. At the default `f=0.01` this is about 14.1
percent of the physical acceleration limit. This geometric estimate is not a
guarantee for sampled velocity differences; the unchanged discrete validator
still decides whether each candidate is admissible.
