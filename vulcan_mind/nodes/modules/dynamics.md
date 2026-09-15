---
id: module.dynamics
label: DynamicEffectors
kind: module
source:
  file: src/dynamics/coupled/force_torque_models.jl
  symbol: DynamicEffectors
inputs: []
outputs:
- id: api
  type: Module
  units: n/a
  description: The exported effector surface — `calcForceTorque`, `wrench`, `wrench_caching!`,
    `environment_requirements`, `solver_partition`, and the model structs `ConstantGravityModel`,
    `InverseSquaredGravityModel`, `InverseSquaredJ2GravityModel`, `NBodyGravityModel`,
    `GravitationalHarmonicsModel`, `SolarRadiationPressureModel`, `MagneticTorqueRodModel`,
    `EddyCurrentDampingModel`, `AerodynamicCoefficientConstant`, `AerodynamicCoefficientfM`,
    `AerodynamicCoefficientNoBallisticFlight`, `BaseThrusterModel`, `RobotArmReactionEffector`.
tags:
- module
charts:
- master
origin: agent
---

# DynamicEffectors

## Purpose
`DynamicEffectors` is the wrapper module that owns every force and torque acting on a
SpaceAGORA spacecraft, plus the rigid-body equations of motion those wrenches drive. It
declares the generic function `calcForceTorque` at `src/dynamics/coupled/force_torque_models.jl:5`
before any submodule is included, so `GravityEffectors`, `AerodynamicEffectors`,
`PerturbationEffectors`, `ThrusterModels`, `GuidanceModels` and `RobotArmReactionEffectors`
can each add methods to one shared symbol. The companion submodules `DynamicsRotational`
(`src/dynamics/rotational/rotational_models.jl`) and `DynamicsTranslational`
(`src/dynamics/translational/translational_models.jl`) turn the summed wrench into state
derivatives consumed by `src/simulation/engine/dynamics_rhs.jl`.

## Theory & Math
Translational motion of a point mass with variable mass in the inertial (J2000) frame:

$$
\dot{\mathbf{r}}_{ii} = \mathbf{v}_{ii}, \qquad
\dot{\mathbf{v}}_{ii} = \frac{1}{m}\,\mathbf{F}_{ii}, \qquad
\dot{m} = \dot{m}_{\text{rate}}
$$

where $\mathbf{r}_{ii}$ is inertial position [m], $\mathbf{v}_{ii}$ inertial velocity [m/s],
$m$ spacecraft mass [kg], $\mathbf{F}_{ii}$ the summed inertial force [N] and
$\dot{m}_{\text{rate}}$ the propellant mass flow [kg/s], negative while thrusting.

Attitude is propagated with Euler's rotational equation including a reaction-wheel
momentum bias, and scalar-last quaternion kinematics on the body rate:

$$
\boldsymbol{\dot\omega}_b = \mathbf{J}^{-1}\Big(\boldsymbol{\tau}_b - \boldsymbol{\omega}_b \times
\big(\mathbf{J}\boldsymbol{\omega}_b + \mathbf{h}_w\big)\Big), \qquad
\dot{\mathbf{q}} = \tfrac{1}{2}\,\mathbf{q}\otimes\begin{bmatrix}\boldsymbol{\omega}_b\\0\end{bmatrix}
$$

where $\boldsymbol{\omega}_b$ is the body angular velocity [rad/s], $\mathbf{J}$ the body
inertia tensor [kg·m²], $\boldsymbol{\tau}_b$ the summed body torque [N·m], $\mathbf{h}_w$
the wheel angular momentum in body axes [N·m·s], and $\mathbf{q} = (q_1,q_2,q_3,q_4)$ the
scalar-last body-to-inertial quaternion [dimensionless].

## Model & Assumptions
- The spacecraft hub is a rigid body: `angular_acceleration` uses a constant `inertia_tensor`
  taken from `spacecraft[i].inertia_tensor`, so flexible-mode and propellant-slosh coupling
  are not represented.
- Quaternion kinematics and Euler's equation are both expressed on the *body* angular rate.
  The docstring in `src/dynamics/rotational/attitude_kinematics.jl` records that the pre-2026-07
  inertial-rate form drifted inertial angular momentum by 68% over 2000 s for a torque-free
  asymmetric tumble with `I = diag(1.5, 1.0, 2.0)` and `|omega| ~ 0.04 rad/s`.
- Wrenches superpose linearly. Each effector returns an `(force_ii, torque_body)` tuple and
  `dynamics_rhs.jl` accumulates them, so aerodynamic-gravity and thruster-plume interactions
  are excluded by construction.
- Translational force is inertial-frame; torque is body-frame. No effector may mix the two.

## Design & Implementation
`force_torque_models.jl` declares `function calcForceTorque end` first, then includes the six
effector files and re-exports their public names. `GravityEffectors` is given a
`@eval GravityEffectors using ..PerturbationEffectors` at line 16 so the legacy
`aerobraking_gravity_force_ii` path can reach `gravity_n_bodies` and `acc_gravity_pines!`
without a circular `include`. Each effector implements two calling conventions: the legacy
`calcForceTorque(model, x::ComponentVector, param::ODEParams, i)` and the sampled
`wrench(model, x::StateSample, env::EnvironmentSample, t)`; `environment_requirements` lets a
model declare that it needs `planet_frame=true` so the RHS builds the planet-fixed rotation
once per stage. `DynamicsTranslational` exposes four RHS writers —
`assign_full_translational_rhs!`, `assign_slow_translational_rhs!`,
`assign_control_only_translational_rhs!` and `assign_force_only_translational_rhs!` — matching
the multirate solver partitions. `DynamicsRotational` exposes `body_torque`,
`body_angular_velocity`, `combine_torques`, `quaternion_derivative` and `angular_acceleration`,
all `@inline` and returning `SVector` values so the RHS allocates nothing per step.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `api` | Module | n/a | — | The exported effector surface — `calcForceTorque`, `wrench`, `wrench_caching!`, `environment_requirements`, `solver_partition`, and the model structs `ConstantGravityModel`, `InverseSquaredGravityModel`, `InverseSquaredJ2GravityModel`, `NBodyGravityModel`, `GravitationalHarmonicsModel`, `SolarRadiationPressureModel`, `MagneticTorqueRodModel`, `EddyCurrentDampingModel`, `AerodynamicCoefficientConstant`, `AerodynamicCoefficientfM`, `AerodynamicCoefficientNoBallisticFlight`, `BaseThrusterModel`, `RobotArmReactionEffector`. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `api` → [[dynamics.aerodynamic_wrench_models__aero_link_atmosphere_query|_aero_link_atmosphere_query]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__aero_workspace_for_sat_bang|_aero_workspace_for_sat!]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__attitude_link_alpha|_attitude_link_alpha]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__ensure_aero_workspace_capacity_bang|_ensure_aero_workspace_capacity!]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__make_aero_scratch_workspace|_make_aero_scratch_workspace]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__multibody_max_threads|_multibody_max_threads]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__multibody_outer_parallel_hint|_multibody_outer_parallel_hint]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__multibody_parallel_mode|_multibody_parallel_mode]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__multibody_thread_decision|_multibody_thread_decision]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__multibody_thread_threshold|_multibody_thread_threshold]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__multibody_use_threads|_multibody_use_threads]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__quaternion_link_alpha|_quaternion_link_alpha]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__simulation_model_module_for_aero|_simulation_model_module_for_aero]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__store_aero_caches_bang|_store_aero_caches!]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__store_vector_cache_bang|_store_vector_cache!]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models__threadid_capacity|_threadid_capacity]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models_aerodynamiccoefficientconstant|AerodynamicCoefficientConstant]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models_aerodynamiccoefficientfm|AerodynamicCoefficientfM]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models_aerodynamiccoefficientnoballisticflight|AerodynamicCoefficientNoBallisticFlight]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models_collect_and_reset_link_wrenches_bang|collect_and_reset_link_wrenches!]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models_set_per_link_atmosphere_bang|set_per_link_atmosphere!]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models_solver_partition|solver_partition]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.aerodynamic_wrench_models_wrench_caching_bang|wrench_caching!]] · `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[dynamics.calc_force_torque|calcForceTorque]] · `module_api` · call · `src/dynamics/coupled/force_torque_models.jl`
- `api` → [[dynamics.cloth_multibody__actuator_torque_child_world|_actuator_torque_child_world]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody__body_offset_velocity|_body_offset_velocity]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody__grid_value|_grid_value]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody__joint_rest|_joint_rest]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody__normalize_state_quaternions__normalize_state_quaternions_bang|_normalize_state_quaternions!]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody__rest_child_parent_quat|_rest_child_parent_quat]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_build_compliant_topology|build_compliant_topology]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_build_rectangular_compliant_grid|build_rectangular_compliant_grid]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_clothmultibody|ClothMultibody]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_compliant_state_vector|compliant_state_vector]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_compliantbody|CompliantBody]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_compliantjoint|CompliantJoint]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_compliantjointactuator|CompliantJointActuator]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_compliantjointload|CompliantJointLoad]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_compliantmultibodymodel|CompliantMultibodyModel]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_compliantmultibodytrajectory|CompliantMultibodyTrajectory]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_complianttopologybuild|CompliantTopologyBuild]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_complianttopologyedge|CompliantTopologyEdge]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_complianttopologynode|CompliantTopologyNode]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_rectangular_prism_inertia|rectangular_prism_inertia]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_simulate_compliant_multibody|simulate_compliant_multibody]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_step_compliant_multibody_implicit_midpoint|step_compliant_multibody_implicit_midpoint]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_multibody_step_compliant_multibody_rk4|step_compliant_multibody_rk4]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics__body_offset_velocity_world|_body_offset_velocity_world]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics__compliance_matrix|_compliance_matrix]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics__rest_child_parent_quat|_rest_child_parent_quat]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics_cloth_reference_state|cloth_reference_state]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_actuators|cloth_robot_arm_actuators]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_end_effector|cloth_robot_arm_end_effector]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_initial_state|cloth_robot_arm_initial_state]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_multibody|cloth_robot_arm_multibody]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics_clothrobotarmdynamics|ClothRobotArmDynamics]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics_clothrobotarmreferencestate|ClothRobotArmReferenceState]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics_clothrobotarmsimulation|ClothRobotArmSimulation]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics_coupled_cloth_robot_arm_state_shape|coupled_cloth_robot_arm_state_shape]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics_initialize_coupled_cloth_robot_arm_state_bang|initialize_coupled_cloth_robot_arm_state!]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.cloth_robot_arm_dynamics_simulate_cloth_robot_arm_plan|simulate_cloth_robot_arm_plan]] · `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`
- `api` → [[dynamics.force_torque_models_dynamiceffectors|DynamicEffectors]] · `module_api` · call · `src/dynamics/coupled/force_torque_models.jl`
- `api` → [[dynamics.perturbations__cannonball_radiation_accel|_cannonball_radiation_accel]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__canonical_harmonics_normalization|_canonical_harmonics_normalization]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__clamp_unit|_clamp_unit]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__convert_harmonics_coefficients_to_full_bang|_convert_harmonics_coefficients_to_full!]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__ensure_nbody_workspace_capacity_bang|_ensure_nbody_workspace_capacity!]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__fully_normalized_legendre_scale|_fully_normalized_legendre_scale]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__get_harmonics_batch_pool|_get_harmonics_batch_pool]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__get_harmonics_batch_pool_cached_bang|_get_harmonics_batch_pool_cached!]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__harmonics_calcforcetorque_with_lpi|_harmonics_calcforcetorque_with_lpi]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__harmonics_flat_batch_kernel_bang|_harmonics_flat_batch_kernel!]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__harmonics_lpi_at_bang|_harmonics_lpi_at!]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__harmonics_lpi_cache_key|_harmonics_lpi_cache_key]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__harmonics_model_cache_key|_harmonics_model_cache_key]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__harmonics_workspace_for_sat_bang|_harmonics_workspace_for_sat!]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__infer_harmonics_reference_radius_m|_infer_harmonics_reference_radius_m]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__interp_vec3_catmull_rom|_interp_vec3_catmull_rom]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__interp_vec3_linear|_interp_vec3_linear]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__lambert_phase_function|_lambert_phase_function]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__lvlh_cascade_torque|_lvlh_cascade_torque]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__magnetic_field_inertial|_magnetic_field_inertial]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__make_harmonics_batch_workspace|_make_harmonics_batch_workspace]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__make_harmonics_scratch_workspace|_make_harmonics_scratch_workspace]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__make_nbody_scratch_workspace|_make_nbody_scratch_workspace]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__nbody_acceleration_ii|_nbody_acceleration_ii]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__nbody_body_force_ii|_nbody_body_force_ii]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__nbody_body_position_from_cache_j2000_m|_nbody_body_position_from_cache_j2000_m]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__nbody_body_position_from_spice_direct_j2000_m|_nbody_body_position_from_spice_direct_j2000_m]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__nbody_body_position_from_spice_j2000_m|_nbody_body_position_from_spice_j2000_m]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__nbody_body_position_from_spice_memoized_j2000_m|_nbody_body_position_from_spice_memoized_j2000_m]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__nbody_workspace_for_sat_bang|_nbody_workspace_for_sat!]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__resolve_third_body_mu|_resolve_third_body_mu]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__spice_lock|_spice_lock]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__spice_rhs_memo_reset_if_stale_bang|_spice_rhs_memo_reset_if_stale!]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__srp_sun_position_from_cache_j2000_m|_srp_sun_position_from_cache_j2000_m]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__srp_sun_position_from_spice_direct_j2000_m|_srp_sun_position_from_spice_direct_j2000_m]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__srp_sun_position_from_spice_j2000_m|_srp_sun_position_from_spice_j2000_m]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__srp_sun_position_from_spice_memoized_j2000_m|_srp_sun_position_from_spice_memoized_j2000_m]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__srp_total_acceleration_ii|_srp_total_acceleration_ii]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_eclipse_area_calc|eclipse_area_calc]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_eddycurrentdampingmodel|EddyCurrentDampingModel]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_get_magnetic_field_dipole|get_magnetic_field_dipole]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_gravitationalharmonicsmodel|GravitationalHarmonicsModel]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_gravity_backbone_kick_acceleration_ii|gravity_backbone_kick_acceleration_ii]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_gravity_backbone_kick_structure|gravity_backbone_kick_structure]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_gravity_backbone_structure|gravity_backbone_structure]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_harmonicsbatchworkspace|HarmonicsBatchWorkspace]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_lvlhcascadeattitudecontrolmodel|LVLHCascadeAttitudeControlModel]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_magnetictorquerodmodel|MagneticTorqueRodModel]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_nbodygravitymodel|NBodyGravityModel]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_planetary_albedo_accel|planetary_albedo_accel]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_planetary_ir_accel|planetary_ir_accel]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.perturbations_solarradiationpressuremodel|SolarRadiationPressureModel]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynamics.point_mass_dynamics_assign_control_only_translational_rhs_bang|assign_control_only_translational_rhs!]] · `module_api` · call · `src/dynamics/translational/point_mass_dynamics.jl`
- `api` → [[dynamics.point_mass_dynamics_assign_force_only_translational_rhs_bang|assign_force_only_translational_rhs!]] · `module_api` · call · `src/dynamics/translational/point_mass_dynamics.jl`
- `api` → [[dynamics.point_mass_dynamics_assign_full_translational_rhs_bang|assign_full_translational_rhs!]] · `module_api` · call · `src/dynamics/translational/point_mass_dynamics.jl`
- `api` → [[dynamics.point_mass_dynamics_assign_slow_translational_rhs_bang|assign_slow_translational_rhs!]] · `module_api` · call · `src/dynamics/translational/point_mass_dynamics.jl`
- `api` → [[dynamics.point_mass_dynamics_mass_derivative|mass_derivative]] · `module_api` · call · `src/dynamics/translational/point_mass_dynamics.jl`
- `api` → [[dynamics.position_kinematics_zero_position_derivative|zero_position_derivative]] · `module_api` · call · `src/dynamics/translational/position_kinematics.jl`
- `api` → [[dynamics.robot_arm_reaction_effector_calcforcetorque|calcForceTorque]] · `module_api` · call · `src/dynamics/coupled/force_torque_models/robot_arm_reaction_effector.jl`
- `api` → [[dynamics.robot_arm_reaction_effector_robotarmreactioneffectors|RobotArmReactionEffectors]] · `module_api` · call · `src/dynamics/coupled/force_torque_models/robot_arm_reaction_effector.jl`
- `api` → [[dynamics.torque_models_body_angular_velocity|body_angular_velocity]] · `module_api` · call · `src/dynamics/rotational/torque_models.jl`
- `api` → [[dynamics.torque_models_body_torque|body_torque]] · `module_api` · call · `src/dynamics/rotational/torque_models.jl`
- `api` → [[dynx.coupled_perturbations_calculate_magnetic_torque|calculate_magnetic_torque]] · `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- `api` → [[dynx.ftm_aerodynamic_effectors_aerodynamiceffectors|AerodynamicEffectors]] · `module_api` · call · `src/dynamics/coupled/force_torque_models/aerodynamic_effectors.jl`
- `api` → [[dynx.ftm_gravity_effectors_gravityeffectors|GravityEffectors]] · `module_api` · call · `src/dynamics/coupled/force_torque_models/gravity_effectors.jl`
- `api` → [[dynx.ftm_guidance_models_guidancemodels|GuidanceModels]] · `module_api` · call · `src/dynamics/coupled/force_torque_models/guidance_models.jl`
- `api` → [[dynx.ftm_perturbation_effectors_perturbationeffectors|PerturbationEffectors]] · `module_api` · call · `src/dynamics/coupled/force_torque_models/perturbation_effectors.jl`
- `api` → [[dynx.ftm_robot_arm_reaction_effector_robotarmreactioneffector|RobotArmReactionEffector]] · `module_api` · call · `src/dynamics/coupled/force_torque_models/robot_arm_reaction_effector.jl`
- `api` → [[dynx.ftm_thruster_models_thrustermodels|ThrusterModels]] · `module_api` · call · `src/dynamics/coupled/force_torque_models/thruster_models.jl`
- `api` → [[dynx.rotational_models_dynamicsrotational|DynamicsRotational]] · `module_api` · call · `src/dynamics/rotational/rotational_models.jl`
- `api` → [[dynx.rotational_rigid_body_dynamics_angular_acceleration|angular_acceleration]] · `module_api` · call · `src/dynamics/rotational/rigid_body_dynamics.jl`
- `api` → [[dynx.rotational_torque_models_combine_torques|combine_torques]] · `module_api` · call · `src/dynamics/rotational/torque_models.jl`
- `api` → [[dynx.translational_models_dynamicstranslational|DynamicsTranslational]] · `module_api` · call · `src/dynamics/translational/translational_models.jl`
- `api` → [[grp.src_dynamics_coupled|dynamics/coupled/]] · `members_in` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- `api` → [[grp.src_dynamics_multibody_cloth|dynamics/multibody_cloth/]] · `members_in` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`
- `api` → [[grp.src_dynamics_rotational|dynamics/rotational/]] · `members_in` · call · `src/dynamics/rotational/torque_models.jl`
- `api` → [[grp.src_dynamics_translational|dynamics/translational/]] · `members_in` · call · `src/dynamics/translational/point_mass_dynamics.jl`
- `api` → [[module.core|SimulationModel]] · `dynamics` · call · `src/core/simulation_model.jl:38-40`
<!-- vulcan:connections:end -->

## Limitations
- `acceleration_from_force` returns the zero vector when mass is non-finite or
  `abs(mass) <= eps(Float64)`, silently masking a depleted-mass bug rather than throwing.
- `mass_derivative` clamps a non-finite `mass_rate` to `0.0`, so a NaN thruster flow rate
  becomes a frozen mass instead of an integrator failure.
- `angular_acceleration` solves `inertia_tensor \ rhs` every call; a singular or
  near-singular inertia tensor produces `Inf`/`NaN` rates with no guard.
- The quaternion is not renormalised inside `quaternion_derivative`; norm drift is bounded
  only by the integrator tolerance and by whatever renormalisation the callback layer applies.
- `ThrusterModels` in this tree is a re-export shim over `...ThrusterModels: BaseThrusterModel`
  and contains no thrust physics of its own.

## Provenance
Mapped from `src/dynamics/coupled/force_torque_models.jl`, with the equations of motion read
from `src/dynamics/rotational/attitude_kinematics.jl`, `src/dynamics/rotational/rigid_body_dynamics.jl`,
`src/dynamics/rotational/torque_models.jl`, `src/dynamics/translational/point_mass_dynamics.jl`
and `src/dynamics/translational/position_kinematics.jl`. Call sites verified in
`src/simulation/engine/dynamics_rhs.jl`.
