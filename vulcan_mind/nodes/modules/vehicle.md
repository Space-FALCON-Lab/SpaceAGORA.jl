---
id: module.vehicle
label: Robotics
kind: module
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: Robotics
inputs: []
outputs:
- id: api
  type: Module
  units: n/a
  description: Robotics model records, forward/inverse kinematics, surface targeting,
    and the vehicle spacecraft, actuator, structure, thermal, and kinematics APIs
    aggregated by SimulationModel.
tags:
- module
charts:
- master
origin: agent
---

# Robotics

## Purpose
The vehicle region supplies the physical configuration carried by a SpaceAGORA spacecraft. Its robotics module defines cloth-arm links, joints, poses, state, forward kinematics, inverse kinematics, and surface-target helpers. Neighboring vehicle modules own spacecraft components and assembly, actuator and thruster models, coordinate transforms, mass and inertia, geometry, thermal behavior, and resource interfaces. `SimulationModel` re-exports these APIs in dependency order.

## Theory & Math
Forward kinematics composes link transforms from the base pose and joint coordinates. Inverse kinematics uses a damped least-squares update, `Δq = (JᵀJ + λ²I)⁻¹Jᵀe`, where `J` is the end-effector Jacobian, `e` is Cartesian error, and `λ` is the damping factor. Structural helpers combine component masses and offsets to obtain center of mass and inertia properties.

## Model & Assumptions
The cloth-arm representation treats links and joints as rigid geometric elements with configured axes, lengths, radii, and masses. Kinematic states use consistent frame conventions and joint-vector length. Vehicle assembly assumes component definitions are complete enough to compute mass properties before dynamics are initialized.

## Design & Implementation
`robotics.jl` declares the arm records and implements `default_cloth_arm_model`, `cloth_fk`, `cloth_fk_state`, `cloth_end_effector_pose`, `cloth_ik`, and `closest_surface_target`. Other vehicle files provide spacecraft containers, assembly graph traversal, actuator hooks, rotations, structure calculations, and thermal models. The dynamics and guidance regions consume the resulting records rather than rebuilding geometry.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `api` | Module | n/a | — | Robotics model records, forward/inverse kinematics, surface targeting, and the vehicle spacecraft, actuator, structure, thermal, and kinematics APIs aggregated by SimulationModel. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `api` → [[grp.src_vehicle_actuators|vehicle/actuators/]] · `members_in` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl`
- `api` → [[grp.src_vehicle_kinematics|vehicle/kinematics/]] · `members_in` · call · `src/vehicle/kinematics/kinematics.jl`
- `api` → [[grp.src_vehicle_robotics|vehicle/robotics/]] · `members_in` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[grp.src_vehicle_spacecraft|vehicle/spacecraft/]] · `members_in` · call · `src/vehicle/spacecraft/assembly.jl`
- `api` → [[grp.src_vehicle_structure|vehicle/structure/]] · `members_in` · call · `src/vehicle/structure/geometry_properties.jl`
- `api` → [[grp.src_vehicle_thermal|vehicle/thermal/]] · `members_in` · call · `src/vehicle/thermal/thermal_models.jl`
- `api` → [[module.core|SimulationModel]] · `vehicle` · call · `src/core/simulation_model.jl:34-54`
- `api` → [[vehicle.assembly_add_facet_bang|add_facet!]] · `module_api` · call · `src/vehicle/spacecraft/assembly.jl`
- `api` → [[vehicle.assembly_add_joint_bang|add_joint!]] · `module_api` · call · `src/vehicle/spacecraft/assembly.jl`
- `api` → [[vehicle.assembly_add_magnet_bang|add_magnet!]] · `module_api` · call · `src/vehicle/spacecraft/assembly.jl`
- `api` → [[vehicle.assembly_add_thruster_bang|add_thruster!]] · `module_api` · call · `src/vehicle/spacecraft/assembly.jl`
- `api` → [[vehicle.assembly_assembly|Assembly]] · `module_api` · call · `src/vehicle/spacecraft/assembly.jl`
- `api` → [[vehicle.components_components|Components]] · `module_api` · call · `src/vehicle/spacecraft/components.jl`
- `api` → [[vehicle.components_create_facet_list|create_facet_list]] · `module_api` · call · `src/vehicle/spacecraft/components.jl`
- `api` → [[vehicle.components_facet|Facet]] · `module_api` · call · `src/vehicle/spacecraft/components.jl`
- `api` → [[vehicle.components_magnet|Magnet]] · `module_api` · call · `src/vehicle/spacecraft/components.jl`
- `api` → [[vehicle.components_reactionwheelassembly|ReactionWheelAssembly]] · `module_api` · call · `src/vehicle/spacecraft/components.jl`
- `api` → [[vehicle.geometry_properties_get_normal_vector|get_normal_vector]] · `module_api` · call · `src/vehicle/structure/geometry_properties.jl`
- `api` → [[vehicle.geometry_properties_get_sa_area|get_SA_area]] · `module_api` · call · `src/vehicle/structure/geometry_properties.jl`
- `api` → [[vehicle.geometry_properties_get_sc_area|get_SC_area]] · `module_api` · call · `src/vehicle/structure/geometry_properties.jl`
- `api` → [[vehicle.geometry_properties_get_spacecraft_length|get_spacecraft_length]] · `module_api` · call · `src/vehicle/structure/geometry_properties.jl`
- `api` → [[vehicle.geometry_properties_get_tangent_vector|get_tangent_vector]] · `module_api` · call · `src/vehicle/structure/geometry_properties.jl`
- `api` → [[vehicle.kinematics_kinematics|Kinematics]] · `module_api` · call · `src/vehicle/kinematics/kinematics.jl`
- `api` → [[vehicle.kinematics_rotate_link|rotate_link]] · `module_api` · call · `src/vehicle/kinematics/kinematics.jl`
- `api` → [[vehicle.mass_properties_get_com|get_COM]] · `module_api` · call · `src/vehicle/structure/mass_properties.jl`
- `api` → [[vehicle.mass_properties_get_inertia_tensor|get_inertia_tensor]] · `module_api` · call · `src/vehicle/structure/mass_properties.jl`
- `api` → [[vehicle.mass_properties_get_spacecraft_mass|get_spacecraft_mass]] · `module_api` · call · `src/vehicle/structure/mass_properties.jl`
- `api` → [[vehicle.mass_properties_set_inertia_tensor_bang|set_inertia_tensor!]] · `module_api` · call · `src/vehicle/structure/mass_properties.jl`
- `api` → [[vehicle.model__initial_condition_apsis_direction_ii|_initial_condition_apsis_direction_ii]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model__initial_condition_lpi|_initial_condition_lpi]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model__initial_condition_oblate_altitude|_initial_condition_oblate_altitude]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model__initial_condition_oblate_surface_radius|_initial_condition_oblate_surface_radius]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model__initial_condition_radius_for_oblate_altitude|_initial_condition_radius_for_oblate_altitude]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model_abstractinitialcondition|AbstractInitialCondition]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model_cartesianinitialcondition|CartesianInitialCondition]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model_controlmodel|ControlModel]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model_dynamicsmodel|DynamicsModel]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model_guidancemodel|GuidanceModel]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model_initialcondition|InitialCondition]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model_joint|Joint]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model_link|Link]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model_navigationmodel|NavigationModel]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.model_spacecraftmodels|SpacecraftModels]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehicle.robotics__ee_position_jacobian|_ee_position_jacobian]] · `module_api` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[vehicle.robotics__normalize_axis|_normalize_axis]] · `module_api` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[vehicle.robotics_closest_surface_target|closest_surface_target]] · `module_api` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[vehicle.robotics_cloth_end_effector_pose|cloth_end_effector_pose]] · `module_api` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[vehicle.robotics_cloth_fk_state|cloth_fk_state]] · `module_api` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[vehicle.robotics_cloth_total_reach|cloth_total_reach]] · `module_api` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[vehicle.robotics_clotharmbasepose|ClothArmBasePose]] · `module_api` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[vehicle.robotics_clotharmjoint|ClothArmJoint]] · `module_api` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[vehicle.robotics_clotharmlink|ClothArmLink]] · `module_api` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[vehicle.robotics_clotharmmodel|ClothArmModel]] · `module_api` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[vehicle.robotics_clotharmstate|ClothArmState]] · `module_api` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[vehicle.robotics_default_cloth_arm_model|default_cloth_arm_model]] · `module_api` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[vehicle.robotics_robotics|Robotics]] · `module_api` · call · `src/vehicle/robotics/robotics.jl`
- `api` → [[vehicle.thermal_models_heatrate_convective|heatrate_convective]] · `module_api` · call · `src/vehicle/thermal/thermal_models.jl`
- `api` → [[vehicle.thermal_models_heatrate_radiative|heatrate_radiative]] · `module_api` · call · `src/vehicle/thermal/thermal_models.jl`
- `api` → [[vehicle.thruster_hooks_integrate_impulse_bang|integrate_impulse!]] · `module_api` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl`
- `api` → [[vehicle.thruster_hooks_schmitt_trigger|schmitt_trigger]] · `module_api` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl`
- `api` → [[vehicle.thruster_hooks_thruster_debug_enabled|thruster_debug_enabled]] · `module_api` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl`
- `api` → [[vehicle.thruster_hooks_thrusterhooks|ThrusterHooks]] · `module_api` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl`
- `api` → [[vehx.actuators_thruster_models_basethrustermodel|BaseThrusterModel]] · `module_api` · call · `src/vehicle/actuators/thruster/thruster_models.jl`
- `api` → [[vehx.actuators_thruster_models_module_thrustermodels|ThrusterModels]] · `module_api` · call · `src/vehicle/actuators/thruster/thruster_models_module.jl`
- `api` → [[vehx.spacecraft_components_thruster|Thruster]] · `module_api` · call · `src/vehicle/spacecraft/components.jl`
- `api` → [[vehx.spacecraft_model_spacecraftmodel|SpacecraftModel]] · `module_api` · call · `src/vehicle/spacecraft/model.jl`
- `api` → [[vehx.structure_models_structure|Structure]] · `module_api` · call · `src/vehicle/structure/structure_models.jl`
- `api` → [[vehx.thermal_models_getheatrate|getHeatRate]] · `module_api` · call · `src/vehicle/thermal/thermal_models.jl`
- `api` → [[vehx.thermal_models_module_vehiclethermalmodels|VehicleThermalModels]] · `module_api` · call · `src/vehicle/thermal/thermal_models_module.jl`
<!-- vulcan:connections:end -->

## Limitations
The inverse-kinematics solver is local and can converge to a configuration that is valid numerically but undesirable mechanically. Singular or poorly conditioned Jacobians depend on damping and stopping tolerances. Geometry and mass-property helpers do not certify structural strength, collision freedom, actuator saturation, or thermal feasibility; those concerns remain downstream validation responsibilities.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` and the vehicle files included through `src/core/simulation_model.jl`.
