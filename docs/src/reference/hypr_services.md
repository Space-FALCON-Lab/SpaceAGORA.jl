# Versioned HYPR service contract

```@docs
SpaceAGORA.HYPRServices
```

`SpaceAGORA.HYPRServices` is the service module through which HYPR's optional
extension plans with SpaceAGORA. It lists the SpaceAGORA types and functions HYPR
may use and the functions HYPR adds methods to. It follows the shared
[optional adapter standard](adapter_contract.md).

This SpaceAGORA source provides contract **3.0.0** and works with HYPR **0.5.0**.
Load `SpaceAGORA` and `HYPR` directly. The setup helper pins the exact HYPR source
revision; the separate compatibility loader is retired
([Upgrading](../user/upgrading.md#Load-HYPR-directly)).

## Loading and failure behavior

Load SpaceAGORA and HYPR in either order. HYPR's extension checks the contract
version and the HYPR version before it defines any method, and checks them again
when it initializes. HYPR becomes available only when initialization succeeds.
Loading the same HYPR version again changes nothing. A different contract or HYPR
version is refused. This contract does not support the no-argument activation of
the old companion package that carried its own implementation, and HYPR's
dependency metadata also excludes companion 0.1.

If you tried an incompatible version, quit that Julia process and start a new one
with a supported pair. Julia can define methods before a package's initialization
fails, and the version check does not remove them. After such a conflict, the
availability check and the planning preflight report HYPR as unavailable. Do not
catch the loading error and keep calling planner methods you imported before.

HYPR's own extension performs these checks and activates the provider.
`SpaceAGORA.hypr_available()` reports successful activation, and planner
preflight refuses an unavailable or conflicted provider. The service module
exposes `CompatibilityError` as an alias of the existing error type.

Availability belongs to one Julia process: each worker loads its own supported
packages.

## Ownership and version policy

The service module gives existing SpaceAGORA names a supported address. It moves
no type, copies no shared calculation and does not change which package owns a
method. Three lists say what each name is for:

- `IMPLEMENTED_FUNCTIONS`: the functions HYPR adds methods to, its extension points;
- `CONSUMED_SERVICES`: the SpaceAGORA types and functions HYPR may use;
- `CLOTH_SERVICES`: three compliant-multibody services under the
  [`Cloth`](#Cloth) namespace.

The types HYPR uses keep their fields and results. HYPR adds methods only for its
own domain to the listed extension points; shared methods stay in SpaceAGORA.

Some listed names begin with an underscore, because the code that already called
them used that spelling. Through the service module they are part of this
versioned contract. Use them only through `HYPRServices`: their original paths in
SpaceAGORA's internal modules are not supported, and a name that is not listed is
not part of the contract.

A change to a listed name's signature, meaning, required fields or ownership
needs a new contract version, reviewed together with HYPR. Contract versions are
matched exactly: while this interface settles, even a change that only adds a
name needs a new reviewed version pair.

## Changes from contract 2.0.0

Contract 3.0.0 narrows the service module to the names that SpaceAGORA's own code
calls or exports and the names HYPR uses. No name that remains changes its
signature, meaning or owner.

- **HYPR's own helpers leave the contract.** 75 RPO and robot-arm helpers that
  only HYPR's extension called are now HYPR's own functions, defined in that
  extension: path cost, path retiming and the acceleration-limited reference
  (`rpo_path_cost`, `rpo_retime_path`, `rpo_retimed_reference`), the PSO,
  refinement and warm-start helpers, and the `_robot_arm_*` swarm, refinement,
  retiming and RRT helpers. SpaceAGORA no longer declares them. Code that reached
  one through `HYPRServices` or SpaceAGORA's internal modules now uses the
  extension's own namespaces, `Base.get_extension(HYPR, :HYPRSpaceAGORAExt).RPO.<name>`
  or `Base.get_extension(HYPR, :HYPRSpaceAGORAExt).RobotArm.<name>`. The extension
  is available after `using SpaceAGORA, HYPR` has activated HYPR; before that, or
  after a refused activation, `Base.get_extension` returns `nothing`.
- **The swarm policy is HYPR's.** The `SwarmPolicy` namespace and the policy
  entries under `RPO` and `RobotArm` are removed. `hypr_iteration_weights`,
  `hypr_material_improvement` and `hypr_protected_particle_mask` are HYPR's
  exported functions, and the extension's `SwarmPolicy` namespace still names them.
- **Three unused services leave the list.** HYPR did not use
  `_robot_arm_segment_distance`, `plan_robot_arm_motion` or
  `rpo_sample_path_bezier`. The functions stay in SpaceAGORA, and
  `plan_robot_arm_motion` stays a public name.
- **What remains.** 20 extension points (HYPR's configured planners, samplers,
  cost components, retiming profile and the RPO planner entry), 50 services
  SpaceAGORA supplies, and the three `Cloth` services.
- **HYPR versions.** This source accepts HYPR 0.5.0 only. HYPR 0.4.0 and earlier
  add methods to helpers that SpaceAGORA no longer declares, so they are refused
  before their methods are defined. HYPR 0.5.0 refuses contract 2.0.0.

`RPOStationGeometry` and `RPOCubeSatGeometry`, which build the station and chaser
of an `RPOReferenceGeometry`, are public SpaceAGORA names from 0.6.0. They are not
part of this contract.

## Changes from contract 1.0.0

Contract 2.0.0 changes one rule for the methods HYPR supplies. No binding is
added, removed or renamed, and SpaceAGORA's own calls are unchanged.

- **The random stream is required.** `rpo_pso_plan_path`, the HYPR-configured
  methods of `rpo_rrt_connect_plan_path`, `rpo_rrt_star_plan_path` and
  `rpo_rrt_connect_bezier_plan_path`, and `plan_robot_arm_motion_hypr` take
  `rng::AbstractRNG` as a required keyword. A call without it fails instead of
  drawing from the task-local global stream. SpaceAGORA's GNC planning kernels
  already require it, and every SpaceAGORA caller passes one.
- **HYPR versions.** HYPR 0.2.0 was the first version to implement the rule
  above. SpaceAGORA 0.5 accepted HYPR 0.4.0 under this contract and refused HYPR
  0.1.0, 0.1.1, 0.2.0 and 0.3.0 before their methods were defined.
- **Argument order.** `_plan_hypr_rpo!(state, planner, request, rng)` already puts
  the state it changes first and the stream last, and is unchanged. The names in
  the lists, underscores included, are unchanged.

This contract covers configured RPO and robot-arm planning. The geometry and
dynamics it uses stay part of SpaceAGORA's simulator.

## Domains and callers

HYPR's names fall into four groups, one for each planning domain. The table says
what each group is for; the [binding inventory](#Explicit-binding-inventory)
below lists every name. The names themselves are defined in
`src/models/gnc/companions/hypr/services.jl`.

| Domain | What HYPR supplies | What SpaceAGORA supplies, and who uses the result |
|---|---|---|
| `Planner` | The HYPR RPO planner that the GNC adapter runs | Request and result types, validation and budgets; SpaceAGORA's configured planner calls HYPR's implementation. |
| `RPO` | Configured planning, sampling, post-refinement, cost components and the retiming profile | Geometry, shared numerical methods and configuration; SpaceAGORA's planner and comparison workflows use the references HYPR returns. |
| `RobotArm` | HYPR's arm planner and its cost components | Arm types, geometry and shared kernels; the arm planning examples and comparison workflows use the plans. |
| `Cloth` | Calls into the compliant-arm services during retiming | Three narrowly exposed dynamics and multibody services in `CLOTH_SERVICES`. |

A name in a group's **HYPR supplies methods** list is an extension point: HYPR
adds methods to that function, but need not write every method of it. For
example, the extension's metric and sampling methods convert HYPR settings and
then call SpaceAGORA's own kernels. Because both packages add methods to the same
function, moving or replacing it changes both. Retiring the loader does not
move those functions or narrow the service contract.

SpaceAGORA's RPO planner adapter, comparison planners and examples call these
names, and so does HYPR's extension. Some HYPR tests still use private
SpaceAGORA paths. These paths are not supported; maintainers find them listed,
with the checks for a boundary change, in
`docs/contributing/companion_boundaries.md` in the
[contributor documentation](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/tree/main/docs/contributing).

## Explicit binding inventory

Each group below is a namespace inside `HYPRServices`, not a copy of a SpaceAGORA
source module. The extension imports each name explicitly; it may not import a
whole internal module.

### RobotArm

**HYPR supplies methods:** `plan_robot_arm_motion_hypr`, `robot_arm_hypr_path_cost_components`

**SpaceAGORA supplies services and types:** `ClothArmBasePose`, `ClothArmModel`, `RobotArmHYPRConfig`, `RobotArmHYPRResult`, `RobotArmPlan`, `RobotArmPlannerConfig`, `RobotArmRRTConnectTree`, `RobotArmSphereObstacle`, `_reference_times`, `_robot_arm_path_length`, `_robot_arm_path_smoothness`, `_validate_robot_arm_hypr_config`, `cloth_fk`, `cloth_ik`, `hypr_rrt_join_paths`, `hypr_rrt_nearest_index`, `hypr_rrt_steer`, `robot_arm_clearance_stats_from_samples`, `robot_arm_sample_hypr_path`

### Planner

**HYPR supplies methods:** `_plan_hypr_rpo!`

**SpaceAGORA supplies services and types:** `HYPRRPOPlanner`, `RPOPlanningRequest`, `RPOPlanningResult`, `RPOReference`, `rpo_planning_budget`, `validate_rpo_result`, `AbstractRNG`

### RPO

**HYPR supplies methods:** `rpo_adaptive_sampling_min_ds_m`, `rpo_fuel_proxy_from_samples`, `rpo_normalized_path_cost_components`, `rpo_path_cost_normalization_refs`, `rpo_post_refine_path`, `rpo_pso_plan_path`, `rpo_reference_from_path`, `rpo_retime_profile`, `rpo_retime_sampling_ds_m`, `rpo_rrt_connect_bezier_plan_path`, `rpo_rrt_connect_plan_path`, `rpo_rrt_star_plan_path`, `rpo_sample_path`, `rpo_sample_path_bezier_adaptive`, `rpo_sample_path_bezier_adaptive_with_params`, `rpo_sample_path_polyline_adaptive`, `rpo_sample_path_with_params`

**SpaceAGORA supplies services and types:** `RPOAdaptiveSamplingSettings`, `RPOPSOConfig`, `RPORRTConnectSettings`, `RPORRTStarSettings`, `RPOReferenceGeometry`, `RPORetimeCurve`, `rpo_adaptive_segment_samples`, `rpo_arc_length_params`, `rpo_clearance_distance_to_station`, `rpo_clearance_stats_from_samples`, `rpo_hcw_fuel_proxy`, `rpo_hypr_refinement_sampling_density_m`, `rpo_hypr_sampling_density_m`, `rpo_obstacle_sigmoid_penalty`, `rpo_obstacle_sigmoid_threshold`, `rpo_path_clearance_stats`, `rpo_path_length`, `rpo_profile_hcw_fuel_proxy`, `rpo_pso_config`, `rpo_resample_polyline_points`, `rpo_retime_samples`, `rpo_retimed_reference_from_profile`, `rpo_sample_path_polyline`, `validate_rpo_pso_config`

### Cloth

`Cloth.ClothRobotArmDynamics` supplies `simulate_cloth_robot_arm_plan` and
`assign_coupled_cloth_robot_arm_rhs!`. `Cloth.ClothMultibody` supplies
`compliant_state_parts`. Each of these small namespaces gives access to those
original functions and to nothing else in the dynamics code.

These are the names in contracts 1.0.0, 2.0.0 and 3.0.0, spelled as they are there. They
still identify the compliant multibody services, although ordinary SpaceAGORA
code now uses `RobotArmModel`, `robot_arm_fk` and the `CompliantRobotArm` names.
The `Cloth` namespace, `CLOTH_SERVICES` and the binding inventory above are
unchanged. In new application code, use the
[preferred root API names](../public_api_policy.md#Compliant-multibody-and-robot-arm-names);
in an adapter that uses this contract, keep the names listed here.
