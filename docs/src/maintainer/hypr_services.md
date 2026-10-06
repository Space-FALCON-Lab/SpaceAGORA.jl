# Versioned HYPR service contract

```@docs
SpaceAGORA.HYPRServices
```

`SpaceAGORA.HYPRServices` is the supported boundary for HYPR's optional extension.
SpaceAGORA **0.2.0** provides contract **1.0.0**, supporting HYPR **0.1.0** and **0.1.1** and the
`SpaceAGORAHYPR` **0.2.0** compatibility package. The setup helper pins the exact
HYPR source revision. A release requires both hosted checks and installation
from the published sources for that version pair.

## Loading and failure behavior

The extension checks the contract and provider version before defining its methods
and checks again during runtime initialization. Only successful initialization
activates HYPR. Loading the same provider/version again is idempotent. Contract or
provider-version mismatch is refused. The old implementation-bearing companion's
no-argument activation protocol is unsupported by this SpaceAGORA contract.
HYPR's dependency metadata also excludes companion 0.1.

If an incompatible implementation was attempted, discard that Julia process and
start a fresh one with the supported pair. Julia can define methods before a
package's initialization fails; this check does not roll them back. Availability
and supported planning preflight fail closed after a provider conflict. Do not
catch a loading failure and continue calling previously imported planner methods.

The compatibility shim checks required extension bindings before constructing aliases
and checks successful initialization and availability again from its `__init__`.
The facade exposes `CompatibilityError` as an alias of the existing error type.

Both core-first and HYPR-first loading are supported. Availability belongs to the
current process; a worker must load its own supported packages. CI precompiles the core, extension and shim, then tests fresh processes with
compiled modules required, including either load order and failed activation.
Local package-image evidence and published-source installation are separate
checks; both are required alongside hosted CI before release acceptance.

## Ownership and version policy

The facade aliases existing SpaceAGORA bindings. It moves no type, duplicates no
shared calculation and changes no method owner. Its explicit
`IMPLEMENTED_FUNCTIONS`, `CONSUMED_SERVICES` and `CLOTH_SERVICES` tuples
distinguish extension points from services. Types in the consumed list retain their existing field and result
contracts. The extension may add only its domain-specific methods to the listed
extension points; shared methods stay in SpaceAGORA.

Some names begin with an underscore because their historical callers use that
spelling. Those listed below are now part of this versioned integration contract
through the facade. Their original internal-module paths are not a supported
integration API. Unlisted internals are outside the contract. A breaking change
to signatures, semantics, required fields or ownership requires a new contract
version and coordinated adapter review. Version 1 is checked exactly; additive
changes also require an explicitly reviewed version pair while this API settles.

This boundary covers configured RPO and robot-arm planning. It does not make
SpaceAGORA geometry or dynamics independent of the simulator.

## Explicit binding inventory

The groups below are namespaces within `HYPRServices`, not copies of the source
modules. The extension imports explicit names only; whole internal-module imports
are not permitted.

### RobotArm

**HYPR supplies methods:** `_robot_arm_control_points`, `_robot_arm_empty_rrt_warmstart_diagnostics`, `_robot_arm_flatten_internal_points`, `_robot_arm_hypr_base_wrench_ratios`, `_robot_arm_hypr_cloth_base_wrench_ratios`, `_robot_arm_hypr_cloth_state_for_reaction`, `_robot_arm_hypr_cull_swarm!`, `_robot_arm_hypr_early_stopping_feasible`, `_robot_arm_hypr_iteration_weights`, `_robot_arm_hypr_link_com_history`, `_robot_arm_hypr_material_improvement`, `_robot_arm_hypr_post_refine_points`, `_robot_arm_hypr_reaction_scale`, `_robot_arm_hypr_reference_times_from_scales`, `_robot_arm_hypr_refinement_better`, `_robot_arm_hypr_retime_reference`, `_robot_arm_hypr_rigid_base_wrench_ratios`, `_robot_arm_plan_from_q_reference`, `_robot_arm_resample_polyline_points`, `_robot_arm_rrt_connect!`, `_robot_arm_rrt_connect_warmstart_path`, `_robot_arm_rrt_extend!`, `_robot_arm_rrt_join_paths`, `_robot_arm_rrt_nearest_index`, `_robot_arm_rrt_path_score`, `_robot_arm_rrt_random_state`, `_robot_arm_rrt_segment_is_safe`, `_robot_arm_rrt_segment_samples`, `_robot_arm_rrt_shortcut_path`, `_robot_arm_rrt_steer`, `_robot_arm_rrt_warmstart_fields`, `_robot_arm_seed_control_points`, `hypr_iteration_weights`, `hypr_material_improvement`, `plan_robot_arm_motion_hypr`, `robot_arm_hypr_path_cost_components`

**SpaceAGORA supplies services and types:** `ClothArmBasePose`, `ClothArmModel`, `RobotArmHYPRConfig`, `RobotArmHYPRResult`, `RobotArmPlan`, `RobotArmPlannerConfig`, `RobotArmRRTConnectTree`, `RobotArmSphereObstacle`, `_reference_times`, `_robot_arm_path_length`, `_robot_arm_path_smoothness`, `_robot_arm_segment_distance`, `_validate_robot_arm_hypr_config`, `cloth_fk`, `cloth_ik`, `hypr_rrt_join_paths`, `hypr_rrt_nearest_index`, `hypr_rrt_steer`, `plan_robot_arm_motion`, `robot_arm_clearance_stats_from_samples`, `robot_arm_sample_hypr_path`

### SwarmPolicy

**HYPR supplies methods:** `hypr_iteration_weights`, `hypr_material_improvement`, `hypr_protected_particle_mask`

**SpaceAGORA supplies services and types:** None.

### Planner

**HYPR supplies methods:** `_plan_hypr_rpo!`

**SpaceAGORA supplies services and types:** `HYPRRPOPlanner`, `RPOPlanningRequest`, `RPOPlanningResult`, `RPOReference`, `rpo_planning_budget`, `validate_rpo_result`, `AbstractRNG`

### RPO

**HYPR supplies methods:** `_rpo_rrt_configured_result`, `_rpo_sampling_settings`, `hypr_iteration_weights`, `hypr_material_improvement`, `hypr_protected_particle_mask`, `rpo_adaptive_pso_config`, `rpo_adaptive_sampling_min_ds_m`, `rpo_estimate_geometry_complexity`, `rpo_fit_bezier_fixed_endpoints`, `rpo_fuel_proxy_dt_s`, `rpo_fuel_proxy_from_samples`, `rpo_manuscript_adaptive_pso_config`, `rpo_manuscript_exploration_score`, `rpo_manuscript_path_cost_components`, `rpo_normalized_path_cost_components`, `rpo_path_cost`, `rpo_path_cost_normalization_refs`, `rpo_position_to_path`, `rpo_post_refine_path`, `rpo_probe_geometry_metrics`, `rpo_pso_bounds`, `rpo_pso_cull_swarm!`, `rpo_pso_early_stopping_feasible`, `rpo_pso_effective_safe_distance`, `rpo_pso_empty_warmstart_diagnostics`, `rpo_pso_iteration_weights`, `rpo_pso_material_improvement`, `rpo_pso_plan_path`, `rpo_pso_project_to_segment`, `rpo_pso_protected_particle_mask`, `rpo_pso_rrt_warmstart_path`, `rpo_pso_stagnation_count_after_learning`, `rpo_pso_station_bounds`, `rpo_pso_tapered_noise_scale`, `rpo_pso_warmstart_bounds`, `rpo_reference_from_path`, `rpo_refine_lower_degree`, `rpo_refine_shortcut_refit`, `rpo_refine_tighten_handles`, `rpo_refinement_bernstein`, `rpo_refinement_better`, `rpo_refinement_clamp_path`, `rpo_refinement_config`, `rpo_refinement_project_to_segment`, `rpo_refinement_sample_params`, `rpo_refinement_segment_is_safe`, `rpo_refinement_segment_samples`, `rpo_refinement_shortcut_samples`, `rpo_retime_available_distance`, `rpo_retime_path`, `rpo_retime_pointwise_speed`, `rpo_retime_profile`, `rpo_retime_sampling_ds_m`, `rpo_retimed_reference`, `rpo_rrt_connect_bezier_plan_path`, `rpo_rrt_connect_plan_path`, `rpo_rrt_star_plan_path`, `rpo_sample_path`, `rpo_sample_path_bezier_adaptive`, `rpo_sample_path_bezier_adaptive_with_params`, `rpo_sample_path_polyline_adaptive`, `rpo_sample_path_with_params`, `rpo_try_accept_refinement`

**SpaceAGORA supplies services and types:** `RPOAdaptiveSamplingSettings`, `RPOPSOConfig`, `RPORRTConnectSettings`, `RPORRTStarSettings`, `RPOReferenceGeometry`, `RPORetimeCurve`, `rpo_adaptive_segment_samples`, `rpo_arc_length_params`, `rpo_clearance_distance_to_station`, `rpo_clearance_stats_from_samples`, `rpo_hcw_fuel_proxy`, `rpo_hypr_refinement_sampling_density_m`, `rpo_hypr_sampling_density_m`, `rpo_obstacle_sigmoid_penalty`, `rpo_obstacle_sigmoid_threshold`, `rpo_path_clearance_stats`, `rpo_path_length`, `rpo_profile_hcw_fuel_proxy`, `rpo_pso_config`, `rpo_resample_polyline_points`, `rpo_retime_samples`, `rpo_retimed_reference_from_profile`, `rpo_sample_path_bezier`, `rpo_sample_path_polyline`, `validate_rpo_pso_config`

### Cloth

`Cloth.ClothRobotArmDynamics` supplies `simulate_cloth_robot_arm_plan` and
`assign_coupled_cloth_robot_arm_rhs!`. `Cloth.ClothMultibody` supplies
`compliant_state_parts`. These narrow namespaces preserve the original function
identities; they do not expose the entire dynamics implementation.
