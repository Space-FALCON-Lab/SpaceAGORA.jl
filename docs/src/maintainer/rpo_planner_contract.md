# RPO planner contract development

The public opt-in [RPO pilot](../user/rpo_planner_pilot.md) connects neutral requests,
references and results to simulation initialization and accepted guidance updates.
The contract and direct baseline load with standard libraries without HYPR source.
This contract-level separation does not yet make the installed package independent
of HYPR. Existing `RPOGuidanceModel`, buffer layouts and controller mathematics
retain their compatibility route.

```@docs
SpaceAGORA.RPOPlannerInterfaces
SpaceAGORA.DirectRPOPlanning
SpaceAGORA.HYPRRPOPlanning
SpaceAGORA.RPOPlannerLifecycle
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

The adapters below preserve HYPR's returned effective settings. Next, add
an opt-in lifecycle, supported external constructor and run-report accessor.
Their acceptance must include callback ordering, initialization, failure,
expiry, multiple controlled spacecraft and repeated/concurrent runs. Publish
top-level names only with the corresponding public API inventory, docstrings
and external example. Optional package extraction follows that pilot.


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

## Internal planner implementations

`SpaceAGORA.DirectRPOPlanning.DirectRPOPlanner` builds a straight segment with
`s(q)=3q²-2q³`, sampled uniformly through the first grid time at or beyond its
chosen duration. Its analytic maxima are `1.5*distance/duration` for speed and
`6*distance/duration²` for acceleration; endpoints are at rest. The duration is
chosen from the reduced limits, rounded up to the request grid, and bounded by
`max_reference_samples` and request lifetime before allocation. A stationary
request receives two identical samples. This reference does not impose the
request's initial velocity as a boundary condition.

The baseline requires an assembly-owned analytic
`segment_clearance_at(start, goal, geometry)` query. It explicitly refuses missing
queries, nonfinite clearance and blocked segments. This checks a direct segment
against the declared geometry approximation, not an exact vehicle mesh or a
search for an alternate route. Independently validate the result with the trusted
point-clearance query before installing it.

`SpaceAGORA.HYPRRPOPlanning.HYPRRPOPlanner` copies its HYPR configuration, then
maps request clearance and sampling interval explicitly. Request limits with
headroom cap configured retiming limits. Conflicting configured minimum/initial
speeds reject. A requested acceleration limit requires the configured
acceleration-limited retimer; the adapter does not silently select another
retiming algorithm. `rrt_on_replan=true` explicitly enables the existing RRT warm
start for a `:replan` request.

The optimizer receives the caller's RNG and mapped configuration. Its returned
effective configuration, including adaptive changes, is retained and used for
retiming. For acceleration-limited retiming the adapter calls the same profile
construction and evaluation helpers, checking the profile duration before uniform
array allocation. Legacy retiming retains its configured step cap and checks
the produced duration before constructing a candidate. Both policies report
`:infeasible/:insufficient_reference_lifetime` when the absolute reference end
exceeds the request lifetime, without adding a time tolerance to this planning
budget. Equality is allowed. The adapter then applies the unchanged validator
and returns `:failed/:reference_rejected` with the validation record for other
output failures. Optimizer exceptions propagate.
A `:candidate` is still subject to lifecycle-owned validation before installation.

Diagnostics distinguish requested, mapped and optimizer-returned configurations;
include the planning budget, path representation, deterministic cost/history and
termination. Matched-input parity means comparing direct HYPR calls with the
adapter's mapped settings. The opt-in reserve can change reference timing from a
legacy call that uses unreduced limits. The old route and its fingerprints remain
unchanged.

Both planners support truth observations in target RTN and use the common
`initialize_planner`, `plan_rpo!` and capability interface. Neither advertises
retiming or restart in this packet. A separate test implementation loads through
the same extension methods. Contract and baseline tests run using standard
libraries with no HYPR source. The public pilot adds callback integration,
lifecycle ownership/expiry, restart refusal at run preparation and external
assembly. Optional-package installation remains a later acceptance gate.

## Lifecycle source ownership

`gnc/guidance/rpo/rpo_planner_module.jl` aggregates the opt-in lifecycle and
bounded public configuration builder. It owns no calculations. The include-chain
and source-completeness gates enforce this separation. `planner_lifecycle.jl`
owns run state, trusted validation/installation, events, failures and expiry;
`planner_configuration.jl` owns the pilot assembly. The engine and controller call
the neutral `SimulationLifecycle` hooks; neither dispatches on a HYPR type.

## Shared metric inputs

The shared path-normalization and finite-difference fuel kernels take required
numeric keyword inputs. Distances are metres, times and specific impulse are
seconds, mass is kilograms, and reference gravity is m/s². They do not own HYPR
defaults or require its configuration type. The existing configuration-based
methods remain in the same defining module through HYPR-owned forwarding in
metric_adapters.jl, preserving current callers and calculations.

Only these two kernels are independently exercised in the minimal shared-module
test. Other metrics retain their geometry/profile dependencies. Comparison
planners, RRT policy and configured retiming still need further separation before
a HYPR-free installation is demonstrated. This internal change adds no root
public API and makes no new physical fuel-model or numeric-type support claim.


### RRT search and HYPR policy ownership

RRT-Connect and RRT* have three-argument internal entry points in
`GuidanceHooks` that accept explicit `bounds`, `evaluate_components` and
`evaluate_cost` keywords. `evaluate_components(path)` must return a named result
with `total`. `refine_path` is optional and returns `(path, cost, improved)`;
`nothing` skips refinement. `edge_is_safe(a, b)` supplies the collision contract
for direct paths, extensions, connections, rewiring and shortcuts. It must be
symmetric because RRT-Connect reverses the goal-tree path when joining. The
legacy default uses the existing shared RPO geometry, clearance and sampling
settings. Adaptive samples depend on traversal direction, so the default is not
guaranteed to give the same answer in reverse. This pre-existing limitation is
preserved here; `path_found` alone does not certify collision clearance in
traversal order. Consumers needing that guarantee must supply a symmetric
predicate and validate the returned path against their collision policy. A
change to the default sampling or goal-tree validation requires separate numerical
review.

The search still uses geometric edge length for tree costs and RRT* rewiring.
The objective callback scores output paths and RRT* history. Callbacks must
agree on constraints and objective, preserve caller-owned input, and avoid
hidden random draws. Explicit `rng` controls search randomness. Failed search
retains the legacy direct-path diagnostic with `path_found=false`; callers must
check that flag before accepting a path. Direct safe paths bypass refinement,
as before. A runtime budget is best-effort wall-clock termination, so seeded
numerical comparisons use an infinite budget and fixed iteration counts.

HYPR's `ext/rpo/rrt_adapters.jl` retains the existing four-argument configuration
methods and result fields, including `config`. HYPR owns its search-box policy,
objective, optional refinement, Bezier fitting, and the decision to request an
RRT warm start. Existing public access and comparison callers remain intact.
The three-argument core omits `config`. Standalone tests load shared geometry,
sampling and tree operations plus RRT, without loading HYPR. The historical
`hypr_` names on shared helpers remain for compatibility with robot-arm users.
This separation does not establish optional HYPR package installation; that
also requires the configured retiming and package-loading work.


### Configured retiming ownership

`src/gnc/shared/rpo/path_retiming.jl` owns two internal calculations in the
existing `GuidanceHooks` module. `rpo_retime_samples(raw_samples, geometry; ...)`
advances along an already sampled path. The five-argument
`rpo_retime_profile(curve, samples, params, clearances, geometry; ...)` builds
the acceleration-limited profile. They require explicit numeric inputs and two
policies: `available_distance(clearance, distance, safe_distance)` and
`pointwise_speed(available_distance, curvature)`. Distances are metres, speed
is m/s, curvature is 1/m, time is seconds and acceleration is m/s². Production geometry is
`RPOReferenceGeometry`. `NavigationHooks` owns the clearance queries, whose
production methods require this type. The sampled-path kernel passes the signed
surface clearance and nearest-station-point distance to `available_distance`.
The profile kernel passes clearance and clearance plus the body margin, computed
from `geometry.station.keepout_radius_m` plus the maximum component of
`geometry.chaser.half_extents_body`. These are the same geometric distance in
exact arithmetic, with different floating-point constructions. The third argument
is the supplied safe distance. `pointwise_speed` receives the resulting available
distance and the sample curvature.

HYPR retains the legacy/manuscript policy distinction, reaction-time and speed
scaling rules, collision-sampling selection, and the choice between the two
retiming routes. Its existing configured methods forward those values and
policies to the shared calculations. The configured reference builder and the
sampling calls remain with HYPR. Existing qualified access, return fields and
module identities are preserved; these internal overloads add no root public API.

A supplied pointwise policy owns the physical speed cap and its order relative
to scaling. `max_speed_mps` in the shared calculation preserves the existing
fallback-speed handling; it does not impose an additional cap on policy output.
The shared calculation preserves minimum-speed floors, near-duplicate handling,
endpoint splitting, Bezier quadrature, forward/backward acceleration passes,
terminal rest and the legacy step-count limit. Warning levels, messages and
values are preserved. The step-cap and invalid-step warnings have different
source-derived tuple-field labels (`max_steps` and `dt_s` now name explicit
inputs instead of configuration expressions). Callbacks must agree with the caller's limits, preserve inputs and avoid hidden random draws.
Existing fallback paths are not a collision-free or feasibility certificate.

The shared geometry and profile-evaluation helpers remain shared and unchanged,
including the general helpers whose present consumers are HYPR-only. Independent
shared-module tests exercise retiming without HYPR definitions. Configured-call
comparisons and existing lifecycle/consumer tests cover compatibility. This
boundary alone does not prove an optional HYPR installation: the package load
chain and dependency declarations still need their separate acceptance work.


## Optional HYPR package boundary

The separate `HYPR` package owns the optimizer; its `HYPRSpaceAGORAExt`
extension owns configured objectives, refinement, RRT/retiming adapters and
robot-arm HYPR execution. HYPR's core has no SpaceAGORA dependency. The extension
uses SpaceAGORA's [versioned service contract](hypr_services.md).
`packages/SpaceAGORAHYPR` is the compatibility shim that loads both packages,
checks the extension and preserves historical module aliases. SpaceAGORA's core
has no required HYPR or shim dependency. OSQP remains a core dependency because
LQ-MPC uses it independently of HYPR.

Core retains configuration/result types, historical compatibility function
bindings and shared mathematical implementations in their original modules.
The HYPR extension imports these explicitly and supplies their configured methods.
The public `HYPRRPOPlanner` type remains in `HYPRRPOPlanning`; its execution delegates
to an extension-provided internal method after checking availability. This is a
coordinated extension of the package family's own contracts. No source is late
included or evaluated inside another package's module.

The HYPR extension registers availability in `__init__`; the shim verifies that
initialization succeeded. Precompilation alone does not activate HYPR.
`require_planner_support` is a baseline no-op that the built-in
HYPR adapter specializes. The public configuration builder and execution preflight
call it, keeping missing-support errors ahead of run output creation. Third-party
planners retain their existing behavior. Do not add a silent fallback planner.

Pure robot-arm path sampling, length, smoothness, segment distance and clearance
remain core-owned in `gnc/shared/robot_arm_geometry.jl`. The HYPR extension's cloth-wrench
bridge explicitly addresses `SpaceAGORA.SimulationModel`, preserving the dynamics
owner that the historical parent-module lookup reached. The mixed replanning file
stays in its original location. The seven general retiming helpers stay shared.

The configured comparison APIs intentionally still use HYPR policies. They require
the companion; the underlying explicit-policy RRT and shared scalar functions do
not. Configuration types can be constructed and inspected without loading execution.
Private compatibility bindings remain for existing qualified research consumers;
they are not a new general plugin API or a promise that every internal helper is
stable forever.

`test/optional_hypr/setup.jl <scratch-directory>` creates separate installation
projects. `core.jl` verifies the dependency graph, absence of the companion, baseline
simulation, shared planners and robot-arm use. `optional.jl core-first` and
`optional.jl hypr-first` run in fresh processes with `JULIA_LOAD_PATH=@:@stdlib`.
The full-feature test harness opts into the local companion explicitly; this
bootstrap is never used by the independent installation proof.

Architecture checks cover the core, compatibility shim and service boundary.
HYPR's coverage checks include the extracted core and extension implementation.
The required `tests` check also waits for the independent installation job.

The service compatibility and loading checks are distinct from the planner-result
contract.
