# RPO planner development

These notes say where the RPO planner interface lives in the source, which module
owns each part, and what the tests check. The interface itself is documented
elsewhere: the [RPO planner contract](../src/reference/rpo_planner_contract.md)
gives the request and reference types, the validation rules, the extension
methods, the built-in planners and the planning headroom, and the
[RPO pilot guide](../src/user/rpo_planner_pilot.md) covers running a case and
installing HYPR. When the code changes, keep those pages and this one consistent.

## Modules

`src/SpaceAGORA.jl` loads five modules for the planner interface:

- `SpaceAGORA.RPOPlannerInterfaces`
  (`src/models/gnc/planning/rpo/contract/rpo_planner.jl`, which includes
  `rpo_reference_validation.jl` and `rpo_planning_headroom.jl`): the request,
  reference, result and capability types, the extension methods,
  `validate_rpo_result` and `RPOPlanningHeadroom`. Apart from Julia standard
  libraries it needs only the general algorithm contract, `GNCAlgorithms`
  (`src/models/gnc/contracts/gnc_algorithms.jl`).
- `SpaceAGORA.DirectRPOPlanning` (`src/models/gnc/planning/rpo/direct_planner.jl`):
  `DirectRPOPlanner`, the direct-segment baseline. It needs only the contract.
- `SpaceAGORA.HYPRRPOPlanning`
  (`src/models/gnc/companions/hypr/rpo_planner_adapter.jl`): `HYPRRPOPlanner`
  and `RetainedPathRPOPlanner`, the adapters to HYPR planning and configured
  retiming.
- `SpaceAGORA.RPOComparisonPlanning`
  (`src/models/gnc/planning/rpo/comparison/planner_contract.jl`): the contract
  methods of the seven comparison planners.
- `SpaceAGORA.RPOPlannerLifecycle` (`src/models/gnc/planning/rpo/lifecycle.jl`):
  the run binding and `make_rpo_configuration`.

The contract loads before `SimulationModel`, because the comparison planner types
in `GuidanceHooks` subtype `AbstractRPOPlanner`;
`src/orchestration/simulation/composition/simulation_model.jl` includes the
contract itself when `SimulationModel` is loaded on its own. The other four
modules load after `SimulationModel`.

## Validation limits and planning headroom

The contract page lists the rules for
[data and ownership](../src/reference/rpo_planner_contract.md#data-and-ownership),
[validation](../src/reference/rpo_planner_contract.md#validation-and-explicit-failures)
and [extension methods](../src/reference/rpo_planner_contract.md#extension-methods).
One of them shapes most of the adapter code: `validate_rpo_result` checks every
reference against the request's physical limits, and the only allowance on its
speed and acceleration limits is `limit_roundoff_rtol`, whose default is also its
ceiling, `128eps(Float64)`. No
planner may loosen the validator to make its output pass. The adapters instead
reserve [planning headroom](../src/reference/rpo_planner_contract.md#adapter-planning-headroom)
before they plan, and reject a reference that still fails rather than retrying
with a larger reserve.

The reason is a difference of definitions. HYPR's acceleration-limited retimer
bounds the acceleration along the path and, through the curvature, the turning,
each on its own. The validator measures the norm of successive velocity-vector
differences per control interval, which combines the two. A HYPR reference timed
to the full limit can therefore exceed the common limit, and the compatibility
tests keep such a case.

## Tests

- `test/unit/gnc/rpo_planner_contract_tests.jl` and
  `test/unit/gnc/direct_rpo_planner_tests.jl` include only `gnc_algorithms.jl`,
  the contract source and, for the baseline, `direct_planner.jl` in standalone
  modules, so they run with `--project=@stdlib`. The contract tests cover the rejection cases, the
  copying at each ownership boundary and the roundoff rule: at a 0.1 m/s² limit,
  an acceleration `8eps` above the limit is accepted and 0.10000004817916151 m/s²
  is rejected.
- `test/unit/gnc/rpo_planner_compatibility_tests.jl` runs seeded HYPR planning
  through the full package, without the adapter, for three cases: seed 741
  (`:legacy` mode, legacy retimer), seed 742 (`:legacy`, acceleration-limited
  retimer) and seed 743 (`:manuscript`, acceleration-limited retimer). It repeats
  each call, compares fingerprints of the paths, costs and references, and
  validates each reference against the physical request limits. Seeds 741 and 743
  pass. Seed 742 is rejected with `:acceleration_limit`: its acceleration exceeds
  the 0.1 m/s² limit by more than a relative 1e-8, far beyond the roundoff
  allowance. The fingerprints are compared within one Julia and dependency
  environment; they are not golden values across versions or platforms.
- `rpo_planner_adapter_tests.jl`, `rpo_planner_lifecycle_tests.jl`,
  `rpo_planner_extensions_tests.jl` and `rpo_comparison_planner_object_tests.jl`,
  all in `test/unit/gnc/`, need the package and HYPR. The lifecycle tests check
  matched-input parity: each installed HYPR reference must equal a direct HYPR
  call made with the adapter's recorded `input_config`. The adapter plans within
  reduced limits, so compare it with direct HYPR calls only at those mapped
  settings.

All of these run in the unit driver, `test/unit/runtests.jl`.

## Planner implementation notes

`DirectRPOPlanner` sets the velocity to zero at both ends of its reference; it
does not use the request's initial velocity as a boundary condition. It checks
the direct segment with an analytic `segment_clearance_at` query, which
`make_rpo_configuration` supplies (`station_segment_clearance`) when the planner
has none. The lifecycle then validates its reference with the trusted point query
`clearance_at`, as it does for every planner.

`HYPRRPOPlanner`'s `plan_rpo!` checks that HYPR is loaded
(`HYPRSupport.require_hypr()`) and calls `_plan_hypr_rpo!`, which HYPR's
extension implements in `ext/rpo_planner_execution.jl`. That method maps the
request onto the planner's copy of the configuration, passes the caller's random
stream to the optimizer, and times the path with the configuration the optimizer
returns, adaptive changes included, never with the original one. With the
acceleration-limited retimer it builds the timing profile with the same sampling
and profile helpers as the configured retimer, and checks the profile's duration
against the sample budget and the request lifetime before it allocates the
uniform arrays. With the legacy retimer the step cap is the configured
`retime_max_steps`, lowered to the sample budget, and the end of the produced
reference is checked against the lifetime before a candidate is built. The
unchanged validator runs last. Optimizer exceptions propagate.

`HYPRRPOPlanner`'s `retime_rpo!` and `RetainedPathRPOPlanner` time a fixed path
and draw no random numbers. Both are in core and share `timed_candidate`, in
`rpo_planner_adapter.jl`, which applies the same request mapping and checks.
`_plan_hypr_rpo!` does not call it and has its own copy of that code, so keep
the two in step.

## Comparison planners

The seven comparison planner types (`HYPRComparisonPlanner` to `STOMPPlanner`)
are defined in `GuidanceHooks`, in
`src/models/gnc/planning/rpo/comparison/planner_objects.jl`, next to the
comparison harness in `planner_comparison.jl`, so that
`RPOPlannerComparisonConfig` can hold them in its `planners` field and build its
default list of all seven. Their contract methods are in `RPOComparisonPlanning`,
which loads after `SimulationModel`, as `HYPRRPOPlanning` does.

Each mode has one method of `_rpo_comparison_planner_path`, and both routes use
it: the harness with the comparison's fixed safety distance
(`RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M`, 0.5 m), and `plan_rpo!` with the
request's mapped settings. The two routes therefore share kernels, settings and
the order of random draws. Both take planner objects; selecting a mode by name
was removed in SpaceAGORA 0.5.0. `normalize_rpo_comparison_planner_type` and
`rpo_comparison_planner_label` still accept names and aliases, but only to give
report rows their canonical identifiers and labels.

### Shared result assembly

RRT-Connect and RRT* share their input conversion (`_rpo_rrt_search_inputs`) and
the result for a direct path that is already safe (`_rpo_rrt_direct_result`), in
`src/models/gnc/planning/rpo/rrt.jl`. Each algorithm keeps its own search tree,
random draws, stopping rule, failure result and refinement. In the direct result,
`path`, `raw_path` and the single `history` entry are the same array.

CHOMP and STOMP share the final refinement and result conversion through
`_rpo_refined_optimizer_result` in
`src/models/gnc/planning/rpo/comparison/trajectory_optimizers.jl`. It keeps the
raw path and cost, the component values, the histories, the iteration count and
the objective value that each search supplies. It is a private helper, not
another planner interface, and the two optimizers keep separate update rules and
termination checks.

`test/unit/gnc/rrt_boundary_tests.jl` covers direct, detoured and failed searches
without HYPR, and `rpo_comparison_planner_object_tests.jl` exercises each mode
with its refusals and metadata.

## Lifecycle source ownership

`src/models/gnc/planning/rpo/lifecycle.jl` defines `RPOPlannerLifecycle`, imports
what it needs and includes its two files, nothing more. The source completeness gate
(`test/gates/ci_src_completeness_contract_gate.jl`) requires its
`# Canonical aggregator: no behavior ownership.` header and refuses any function
or type definition in it.

`planner_lifecycle.jl` owns the binding, `PlannerGuidance`: its run state,
planning events and triggers, failures (`RPOPlanningError`), reference expiry,
checkpoint state, installation into the plan buffer, and `rpo_run_report`. Each
request goes through the shared driver `run_algorithm_request!`
(`src/models/gnc/contracts/gnc_algorithms.jl`), which validates the planner's
result before it is installed. `configuration.jl` owns the pilot assembly,
`make_rpo_configuration`, and the `RPOReplanningSpheres` geometry provider. The
engine calls the binding through neutral hooks (`SimulationLifecycle`,
`GNCRuntime` and `GuidanceHooks.update_guidance!`), and the LQ-MPC controller
(`RPOMPCControlModel`) reads the plan buffer; neither dispatches on a HYPR type.

## Shared calculations and HYPR policies

Core keeps the RPO path calculations that do not depend on HYPR's configuration
in `GuidanceHooks`, with explicit inputs. HYPR's extension adds the configured
methods that read an `RPOPSOConfig` and forward to them. The calculations keep
`GuidanceHooks` as their owner, so qualified names such as
`SimulationModel.GuidanceHooks.rpo_retime_samples` keep working
(`test/optional/hypr/optional.jl` checks these identities). None of them is
exported from `SpaceAGORA`.

### Metrics

`src/models/gnc/planning/rpo/path_metrics.jl` holds
`rpo_path_cost_normalization_refs(points; cost_ref_distance_m, sample_ds_m, tf_s, mass_kg, isp_s, g0_mps2)`
and the finite-difference fuel proxy
`rpo_fuel_proxy_from_samples(samples; tf_s, mass_kg, isp_s, g0_mps2)`. All the
keywords are required: distances in metres, times and specific impulse in
seconds, mass in kilograms and reference gravity in m/s². HYPR's
`ext/rpo/metric_adapters.jl` adds the methods that take an `RPOPSOConfig` and
pass its fields. `test/unit/gnc/shared_metrics_tests.jl` loads `path_metrics.jl`
alone and tests only these two kernels; the file's other metrics still depend on
geometry and profile code. The kernels compute in `Float64` and are not a
physical fuel model.

### RRT search

`rpo_rrt_connect_plan_path` and `rpo_rrt_star_plan_path` in
`src/models/gnc/planning/rpo/rrt.jl` take three positional arguments (start, goal
and geometry) and explicit keywords: `bounds`, the search box;
`evaluate_components(path)`, which returns a named result with `total`;
`evaluate_cost(path)`, the scalar score; an optional `refine_path(path)` returning
`(path, cost, improved)`, where `nothing` skips refinement; a required `rng`; and
`edge_is_safe(a, b)`, the collision test for the direct path, extensions,
connections, rewiring and shortcuts.

`edge_is_safe` must be symmetric, because RRT-Connect reverses the goal-tree path
when it joins the two trees. The default predicate uses the shared RPO geometry,
clearance and adaptive sampling settings, and its samples depend on the direction
of traversal, so it is not guaranteed to give the same answer in reverse.
`path_found` alone therefore does not certify clearance in traversal order; a
caller that needs that guarantee must supply a symmetric predicate and validate
the returned path against its own collision policy. A change to the default
sampling or to the goal-tree validation needs its own numerical check.

The trees grow, and RRT* rewires, by geometric edge length; the objective
callbacks score the output paths and the RRT* cost history. Callbacks must agree
on constraints and objective, leave caller-owned inputs unchanged and make no
hidden random draws; `rng` owns the search's draws. A failed search returns the
direct path as a diagnostic with `path_found=false`, and callers must check that
flag before using the path. A direct path that is already safe is returned
without refinement. The `max_runtime_s` budget is checked between iterations, so
it ends a search only approximately; seeded numerical comparisons use an infinite
budget and fixed iteration counts.

HYPR's `ext/rpo/rrt_adapters.jl` adds the four-argument methods that take an
`RPOPSOConfig` and return the result with its `config` field. HYPR supplies the
search box (`rpo_pso_bounds`), the objective, the refinement and the Bezier
fitting, and decides when to request an RRT warm start.
`test/unit/gnc/rrt_boundary_tests.jl` loads the shared geometry, sampling and
tree operations and `rrt.jl` without HYPR. The shared helpers that still carry a
`hypr_` prefix, in `src/models/gnc/planning/common/path_geometry.jl` and
`rrt_tree_operations.jl`, are used by RPO and robot-arm planning alike and keep
their names for qualified access through `HYPRUtils`.

### Retiming

`src/models/gnc/planning/rpo/path_retiming.jl` holds two retiming calculations.
`rpo_retime_samples(raw_samples, geometry; ...)` advances along an already
sampled path, and the five-argument
`rpo_retime_profile(curve, samples, params, clearances, geometry; ...)` builds
the acceleration-limited profile. Both take explicit numeric limits and two
policy callbacks, `available_distance(clearance, distance, safe_distance)` and
`pointwise_speed(available_distance, curvature)`. Distances are in metres, speeds
in m/s, curvature in 1/m, times in seconds and acceleration in m/s².

The production geometry is `RPOReferenceGeometry`. The clearance queries belong
to `NavigationHooks` (`src/models/gnc/planning/rpo/geometry/clearance.jl`), and
their methods are defined for that type. `rpo_retime_samples` passes `available_distance`
the signed surface clearance and the distance to the nearest station point.
`rpo_retime_profile` passes the clearance and the clearance plus a body margin,
`geometry.station.keepout_radius_m` plus the largest component of
`geometry.chaser.half_extents_body`. In exact arithmetic the two second arguments
are the same distance; in floating point they are computed differently. The third
argument is the supplied safe distance, and `pointwise_speed` receives the
resulting available distance and the sample's curvature.

The pointwise policy owns the physical speed cap and where it applies relative to
the speed scaling. The `max_speed_mps` keyword of the shared calculations caps
only the fallback speed and puts no further cap on the policy's output. The shared
calculations hold the minimum-speed floor, near-duplicate removal, the split of a
two-point path, Bezier quadrature, the forward and backward acceleration passes,
the terminal rest and, in `rpo_retime_samples`, the `max_steps` limit. Callbacks
must agree with the caller's limits, leave inputs unchanged and make no hidden
random draws. The fallback paths are not a certificate of clearance or
feasibility.

HYPR's `ext/rpo/path_retiming.jl` keeps the configured side: the `:legacy` and
`:manuscript` policies (`rpo_retime_available_distance`,
`rpo_retime_pointwise_speed`), the reaction-time and speed-scaling rules, and the
collision-sampling spacing (`rpo_retime_sampling_ds_m`). Its
`ext/rpo/rpo_reference_trajectory.jl` chooses between the two retiming routes
(`rpo_reference_from_path`). `test/unit/gnc/shared_retiming_tests.jl` loads the
shared files without any HYPR definitions and runs both calculations.

## Optional HYPR package boundary

HYPR is a separate package, pinned under `packages/robotics/HYPR.jl`, that owns
the swarm optimizer. Its `HYPRSpaceAGORAExt` extension, which Julia loads once
both `SpaceAGORA` and `HYPR` are loaded, holds the configured objectives,
refinement, sampling, RRT and retiming adapters, the planner execution
`_plan_hypr_rpo!` and robot-arm HYPR planning. SpaceAGORA is only a weak
dependency of HYPR, and SpaceAGORA's `Project.toml` does not list HYPR at all.
The extension reaches core code through the versioned
[HYPR services](../src/reference/hypr_services.md) facade,
`SpaceAGORA.HYPRServices`. It includes only its own files and evaluates nothing
inside SpaceAGORA's modules. OSQP stays a core dependency because the LQ-MPC
controller uses it, with or without HYPR.

Core keeps the configuration and result types (`RPOPSOConfig`,
`RobotArmHYPRConfig`, `RobotArmHYPRResult`), the function bindings that qualified
callers use, and the shared calculations, all in their original modules; the
extension adds the configured methods to those functions. A configuration can be
built and inspected without HYPR. `HYPRRPOPlanner` lives in `HYPRRPOPlanning`, and
its `plan_rpo!` checks availability before delegating to the extension.
Comparison planning configured with HYPR's policies needs the companion; the
explicit-policy RRT entry points and the shared scalar kernels do not. The
compatibility bindings kept for qualified research callers are not a general
plugin API, and internal helpers are not promised to stay stable.

The extension's `__init__` calls `HYPRServices.activate!`, which checks HYPR's
version and the service contract version before recording HYPR as available
(`hypr_available()`); precompiling HYPR does not activate it.
`require_planner_support` does nothing for a planner by default.
`HYPRRPOPlanner`, `RetainedPathRPOPlanner` and the comparison planners specialize
it to raise `HYPRUnavailableError` when HYPR is not loaded.
`make_rpo_configuration` and the binding's `preflight_guidance` call it, so the
error comes before the run creates any output. Other planners keep the default.
Do not add a silent fallback planner.

Robot-arm path sampling, length, smoothness, segment distance and clearance stay
in core, in `src/models/gnc/planning/robot_arm/path_geometry.jl`. The extension's
cloth-coupled base-wrench estimate reaches `SimulationModel.ClothRobotArmDynamics`
and `SimulationModel.ClothMultibody` through `HYPRServices.Cloth`, and falls back
to the rigid-base estimate when they are not defined.
`src/models/gnc/companions/hypr/replanning.jl`, which mixes the neutral
replanning geometry and decisions with HYPR-configured replanning, stays in core,
as do the retimer's helpers in `path_geometry.jl` and `profile_evaluation.jl`
under `src/models/gnc/planning/rpo/`.

The HYPR service contract has its own version (`HYPRServices.CONTRACT_VERSION`),
separate from the planner contract.

### Installation tests

`test/optional/hypr/setup.jl <scratch-directory>` builds two projects in that
directory: `core`, without HYPR, and `hypr`, with it. `test/optional/hypr/core.jl`,
run in `core`, checks that neither HYPR nor the retired `SpaceAGORAHYPR` loader is
installed or loaded, that selecting HYPR raises `HYPRUnavailableError` before a
run writes output, and that the baseline simulation, the shared planners and the
quintic robot-arm planner work. `test/optional/hypr/optional.jl core-first` and
`optional.jl hypr-first`, run in `hypr`, load the two packages in each order and
check that type and function identities are kept. The
`optional-hypr-installation` job in `.github/workflows/julia-ci.yml` runs each in
its own process with `JULIA_LOAD_PATH=@:@stdlib`, and the required `tests` job
waits for it.

The unit and integration suites load HYPR differently, through
`test/helpers/load_hypr.jl`, which the installation tests do not use.
`test/unit/gnc/hypr_services_tests.jl` checks the service facade, and HYPR's own
CI measures coverage of both its core and its SpaceAGORA extension.
