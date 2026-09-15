---
id: module.core
label: SimulationModel
kind: module
source:
  file: src/core/simulation_model.jl
  symbol: SimulationModel
inputs:
- id: dynamics
  type: Module
  units: n/a
  required: true
  description: ClothMultibody, ClothRobotArmDynamics, DynamicsRotational, DynamicsTranslational
    and DynamicEffectors, included at lines 38, 40, 95, 99 and 103.
- id: environment
  type: Module
  units: n/a
  required: true
  description: Planets, EphemeridesModels and EnvironmentModels plus the gravity and
    aerodynamic effector files, included at lines 17, 48, 50, 82, 105 and 106.
- id: gnc
  type: Module
  units: n/a
  required: true
  description: CommandTypes, RobotArmPlanning, GuidanceModels, NavigationHooks, GuidanceHooks
    and ControlHooks, included at lines 29, 36, 45, 121, 124 and 127.
- id: io
  type: Module
  units: n/a
  required: true
  description: IOConfig, IOSerialization and IOOutputs, included at lines 109, 111
    and 113, owning configuration parsing, state serialization and result writing.
- id: mission
  type: Module
  units: n/a
  required: true
  description: AerobrakingPolicy, included at line 117, owning the mission-operations
    policy types consumed by guidance and control hooks.
- id: vehicle
  type: Module
  units: n/a
  required: true
  description: Robotics, ThrusterModels, Components, SpacecraftModels, Assembly, Kinematics,
    ThrusterHooks, Structure and VehicleThermalModels, included at lines 34 to 131.
outputs:
- id: api
  type: Module
  units: n/a
  description: The union namespace produced by the twenty-eight @reexport using statements,
    exposing spacecraft, environment, effector, sampling and configuration symbols
    under one module.
tags:
- module
charts:
- master
origin: agent
---

# SimulationModel

## Purpose
`SimulationModel` is the canonical aggregator described by its own first line: it owns no
behavior. It fixes the load order for every physics, vehicle, guidance and IO sub-module, and
republishes their exported symbols through `Reexport.@reexport` so that downstream code holds a
single dependency instead of twenty-eight.

## Model & Assumptions
- Load order encodes the real dependency graph. Abstract types come first (line 24), then
  sampling types (26), commands (29), robotics and dynamics (34 to 41), actuators and guidance
  (43 to 46), planets and ephemerides (48 to 51), vehicle components and assembly (54 to 75),
  configuration types (78), physical models (82), dynamics (95 to 106), IO (109 to 113), mission
  policy (117) and finally the GNC hook surfaces (121 to 128) and callbacks (139).
- Lines 11 and 12 assume `RuntimeServices` may not yet exist in the parent module and
  conditionally `Base.include` it into `parentmodule(@__MODULE__)`, so this file is loadable
  both from `SpaceAGORA.jl` and standalone.
- Three files are included without a matching `@reexport`: `hypr_utils.jl` (32),
  `parallel_policy.jl` (92), and the gravity and aerodynamic effector files (105, 106). Those
  define methods on already-exported generic functions rather than new public names.
- `@reexport` assumes no two sub-modules export the same identifier.

## Design & Implementation
The body is a straight sequence of `include(joinpath(@__DIR__, "..", ...))` calls paired with
`@reexport using .Submodule`. Utility includes at lines 15 to 17 splice
`quaternion_utils.jl`, `reference_system_config.jl` and `planet_shapes.jl` directly into
`SimulationModel`'s own scope rather than into a nested module, so `quat_mult`, `rot`, `hat`,
`error_quaternion`, `project_unit_quaternion`, `qToEulerAngles` and `dcm_to_quaternion` become
`SimulationModel` methods and are visible to every sub-module that resolves names upward through
`..`. `NoGramPresets` at line 135 is the clearest example of that upward resolution: it declares
`using ..AbstractTypes`, `..Planets`, `..EnvironmentModels`, `..EphemeridesModels`,
`..VehicleThermalModels` and `..SimConfig`, all of which exist only because earlier lines already
included them.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dynamics` | Module | n/a | yes | ClothMultibody, ClothRobotArmDynamics, DynamicsRotational, DynamicsTranslational and DynamicEffectors, included at lines 38, 40, 95, 99 and 103. |
| in | `environment` | Module | n/a | yes | Planets, EphemeridesModels and EnvironmentModels plus the gravity and aerodynamic effector files, included at lines 17, 48, 50, 82, 105 and 106. |
| in | `gnc` | Module | n/a | yes | CommandTypes, RobotArmPlanning, GuidanceModels, NavigationHooks, GuidanceHooks and ControlHooks, included at lines 29, 36, 45, 121, 124 and 127. |
| in | `io` | Module | n/a | yes | IOConfig, IOSerialization and IOOutputs, included at lines 109, 111 and 113, owning configuration parsing, state serialization and result writing. |
| in | `mission` | Module | n/a | yes | AerobrakingPolicy, included at line 117, owning the mission-operations policy types consumed by guidance and control hooks. |
| in | `vehicle` | Module | n/a | yes | Robotics, ThrusterModels, Components, SpacecraftModels, Assembly, Kinematics, ThrusterHooks, Structure and VehicleThermalModels, included at lines 34 to 131. |
| out | `api` | Module | n/a | — | The union namespace produced by the twenty-eight @reexport using statements, exposing spacecraft, environment, effector, sampling and configuration symbols under one module. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `dynamics` · call · `src/core/simulation_model.jl:38-40`
- [[module.environment|EnvironmentModels]] · `api` → `environment` · call · `src/core/simulation_model.jl:17-50`
- [[module.gnc|CommandTypes]] · `api` → `gnc` · call · `src/core/simulation_model.jl:29-36`
- [[module.io|IOConfig]] · `api` → `io` · call · `src/core/simulation_model.jl:109-113`
- [[module.mission|AerobrakingPolicy]] · `api` → `mission` · call · `src/core/simulation_model.jl:117-117`
- [[module.vehicle|Robotics]] · `api` → `vehicle` · call · `src/core/simulation_model.jl:34-54`

**Downstream**

- `api` → [[core.abstract_types_abstractcontroleffectormodel|AbstractControlEffectorModel]] · `module_api` · call · `src/core/types/abstract_types.jl`
- `api` → [[core.abstract_types_abstractdensitymodel|AbstractDensityModel]] · `module_api` · call · `src/core/types/abstract_types.jl`
- `api` → [[core.abstract_types_abstractephemeridesmodel|AbstractEphemeridesModel]] · `module_api` · call · `src/core/types/abstract_types.jl`
- `api` → [[core.abstract_types_abstractforcetorquemodel|AbstractForceTorqueModel]] · `module_api` · call · `src/core/types/abstract_types.jl`
- `api` → [[core.abstract_types_abstractguidancemodel|AbstractGuidanceModel]] · `module_api` · call · `src/core/types/abstract_types.jl`
- `api` → [[core.abstract_types_abstractplanet|AbstractPlanet]] · `module_api` · call · `src/core/types/abstract_types.jl`
- `api` → [[core.abstract_types_abstractthermalmodel|AbstractThermalModel]] · `module_api` · call · `src/core/types/abstract_types.jl`
- `api` → [[core.abstract_types_abstractthrustermodel|AbstractThrusterModel]] · `module_api` · call · `src/core/types/abstract_types.jl`
- `api` → [[core.compat_model_codes__compat_enum_parse|_compat_enum_parse]] · `module_api` · call · `src/core/types/compat_model_codes.jl`
- `api` → [[core.compat_model_codes_legacymodelcodes|LegacyModelCodes]] · `module_api` · call · `src/core/types/compat_model_codes.jl`
- `api` → [[core.effector_sampling_atmospheresample|AtmosphereSample]] · `module_api` · call · `src/core/types/effector_sampling.jl`
- `api` → [[core.effector_sampling_effectorsampling|EffectorSampling]] · `module_api` · call · `src/core/types/effector_sampling.jl`
- `api` → [[core.effector_sampling_gravity_backbone_kick_structure|gravity_backbone_kick_structure]] · `module_api` · call · `src/core/types/effector_sampling.jl`
- `api` → [[core.effector_sampling_gravity_backbone_structure|gravity_backbone_structure]] · `module_api` · call · `src/core/types/effector_sampling.jl`
- `api` → [[core.effector_sampling_planetframesample|PlanetFrameSample]] · `module_api` · call · `src/core/types/effector_sampling.jl`
- `api` → [[core.effector_sampling_solarephemerissample|SolarEphemerisSample]] · `module_api` · call · `src/core/types/effector_sampling.jl`
- `api` → [[core.effector_sampling_solver_partition|solver_partition]] · `module_api` · call · `src/core/types/effector_sampling.jl`
- `api` → [[core.effector_sampling_statesample|StateSample]] · `module_api` · call · `src/core/types/effector_sampling.jl`
- `api` → [[core.effector_sampling_thirdbodyephemerissample|ThirdBodyEphemerisSample]] · `module_api` · call · `src/core/types/effector_sampling.jl`
- `api` → [[core.no_gram_presets_nogrampresets|NoGramPresets]] · `module_api` · call · `src/core/state/no_gram_presets.jl`
- `api` → [[core.quaternion_utils_dcm_to_quaternion|dcm_to_quaternion]] · `module_api` · call · `src/core/numerics/quaternion_utils.jl`
- `api` → [[core.quaternion_utils_error_quaternion|error_quaternion]] · `module_api` · call · `src/core/numerics/quaternion_utils.jl`
- `api` → [[core.quaternion_utils_hat|hat]] · `module_api` · call · `src/core/numerics/quaternion_utils.jl`
- `api` → [[core.quaternion_utils_qtoeulerangles|qToEulerAngles]] · `module_api` · call · `src/core/numerics/quaternion_utils.jl`
- `api` → [[core.quaternion_utils_quat_mult|quat_mult]] · `module_api` · call · `src/core/numerics/quaternion_utils.jl`
- `api` → [[core.reference_system__body_fixed_state_xform|_body_fixed_state_xform]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system__body_fixed_to_j2000_state|_body_fixed_to_j2000_state]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system__j2000_to_body_fixed_state|_j2000_to_body_fixed_state]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system__planet_flattening|_planet_flattening]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system__rtn_rate_rad_s|_rtn_rate_rad_s]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system__safe_acos|_safe_acos]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system__spice_body_fixed_frame|_spice_body_fixed_frame]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system__spice_lock|_spice_lock]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system_alfadeltartor|alfadeltartor]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system_config_cartesian|cartesian]] · `module_api` · call · `src/core/state/reference_system_config.jl`
- `api` → [[core.reference_system_config_clock|clock]] · `module_api` · call · `src/core/state/reference_system_config.jl`
- `api` → [[core.reference_system_config_h_lan_lon|H_LAN_LON]] · `module_api` · call · `src/core/state/reference_system_config.jl`
- `api` → [[core.reference_system_config_r_ra_dec|R_RA_DEC]] · `module_api` · call · `src/core/state/reference_system_config.jl`
- `api` → [[core.reference_system_config_referencesystems|ReferenceSystems]] · `module_api` · call · `src/core/state/reference_system_config.jl`
- `api` → [[core.reference_system_config_udunue|uDuNuE]] · `module_api` · call · `src/core/state/reference_system_config.jl`
- `api` → [[core.reference_system_latlongtooe|latlongtoOE]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system_latlongtor|latlongtor]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system_orbital_elements_to_lvlh_quaternion|orbital_elements_to_lvlh_quaternion]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system_r_pintor_i|r_pintor_i]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system_rotate_vector_by_quaternion|rotate_vector_by_quaternion]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system_rtn_dcm_from_inertial|rtn_dcm_from_inertial]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system_rtn_to_inertial_relative_state|rtn_to_inertial_relative_state]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system_rtoalfadeltar|rtoalfadeltar]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.reference_system_rtolatlongrad|rtolatlongrad]] · `module_api` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[core.runtime_types__typed_nothing_vector|_typed_nothing_vector]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_aerodynamics|Aerodynamics]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_aeroscratchworkspace|AeroScratchWorkspace]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_closed_form|Closed_form]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_cnf|Cnf]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_configtypes|ConfigTypes]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_controller|Controller]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_engines|Engines]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_forces|Forces]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_gramtrackcache|GramTrackCache]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_harmonicsscratchworkspace|HarmonicsScratchWorkspace]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_initial_condition|Initial_condition]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_initialparameters|InitialParameters]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_intermediatesolution|IntermediateSolution]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_mission|Mission]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_model|Model]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_nbodyephemeriscache|NBodyEphemerisCache]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_nbodyscratchworkspace|NBodyScratchWorkspace]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_orientation|Orientation]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_performance|Performance]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_physical_properties|Physical_properties]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_planetframeephemeriscache|PlanetFrameEphemerisCache]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_rhsplanenvconfig|RhsPlanEnvConfig]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_savecache|SaveCache]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_sharedbuffers|SharedBuffers]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_simulation|Simulation]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_solution|Solution]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_spicerhsmemo|SpiceRhsMemo]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_spiceruntimecounters|SpiceRuntimeCounters]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_srpsunephemeriscache|SRPSunEphemerisCache]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.runtime_types_vacuumpredictedgramcache|VacuumPredictedGRAMCache]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[core.simulation_configuration__parse_mission_type|_parse_mission_type]] · `module_api` · call · `src/core/state/simulation_configuration.jl`
- `api` → [[core.simulation_configuration__warn_deprecated_config_enabled|_warn_deprecated_config_enabled]] · `module_api` · call · `src/core/state/simulation_configuration.jl`
- `api` → [[core.simulation_configuration__warn_deprecated_mission_type_input_bang|_warn_deprecated_mission_type_input!]] · `module_api` · call · `src/core/state/simulation_configuration.jl`
- `api` → [[core.simulation_configuration_initialtime|InitialTime]] · `module_api` · call · `src/core/state/simulation_configuration.jl`
- `api` → [[core.simulation_configuration_simconfig|SimConfig]] · `module_api` · call · `src/core/state/simulation_configuration.jl`
- `api` → [[grp.src_core_interfaces|core/interfaces/]] · `members_in` · call · `src/core/interfaces/reference_system.jl`
- `api` → [[grp.src_core_numerics|core/numerics/]] · `members_in` · call · `src/core/numerics/quaternion_utils.jl`
- `api` → [[grp.src_core_simulation_model_jl|core/simulation_model.jl]] · `members_in` · call · `src/core/simulation_model.jl`
- `api` → [[grp.src_core_state|core/state/]] · `members_in` · call · `src/core/state/no_gram_presets.jl`
- `api` → [[grp.src_core_types|core/types/]] · `members_in` · call · `src/core/types/abstract_types.jl`
- `api` → [[module.spaceagora|SpaceAGORA]] · `core` · call · `src/SpaceAGORA.jl:9-9`
- `api` → [[parcore.abstract_types_abstracttypes|AbstractTypes]] · `module_api` · call · `src/core/types/abstract_types.jl`
- `api` → [[parcore.compat_model_codes_legacydensitymodelcode|LegacyDensityModelCode]] · `module_api` · call · `src/core/types/compat_model_codes.jl`
- `api` → [[parcore.reference_system_config_oe|OE]] · `module_api` · call · `src/core/state/reference_system_config.jl`
- `api` → [[parcore.runtime_types_odeparams|ODEParams]] · `module_api` · call · `src/core/types/runtime_types.jl`
- `api` → [[parcore.simulation_model_simulationmodel|SimulationModel]] · `module_api` · call · `src/core/simulation_model.jl`
<!-- vulcan:connections:end -->

## Limitations
The aggregator is a single point of serialization for compile time: touching any included file
invalidates the whole `SimulationModel` precompile unit. Because utility files are spliced into
module scope rather than namespaced, a name defined in `quaternion_utils.jl` can collide with a
sub-module export without an obvious diagnostic. The `..`-relative `using` statements in leaf
modules hard-code this exact parent structure, so a sub-module cannot be loaded in isolation for
unit testing without reproducing the include order. There is no conditional inclusion: a run
that needs neither robotics nor cloth multibody dynamics still pays their load cost.

## Provenance
Mapped from `src/core/simulation_model.jl`, with socket evidence read from
`src/core/types/effector_sampling.jl` and `src/core/state/no_gram_presets.jl`.
