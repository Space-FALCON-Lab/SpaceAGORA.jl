# Public API and compatibility

Applications and extensions should depend only on the
root `SpaceAGORA` module, `using SpaceAGORA`, and on the versioned facades
described below.
The generated public API pages list the supported names. Anything else is
internal and can move between releases.

## Stable surface

The generated [Public API](generated/public_api.md) reference and the pages it
links to, among them the [Simulation Setup API](generated/simulation_setup_api.md),
[Terrain API](generated/terrain_api.md),
[Mesh Aerodynamics API](generated/mesh_aerodynamics_api.md) and
[Visualization API](generated/visualization_api.md), document the supported
exported names. They include:

- the simulation entry points, such as `run_simulation`;
- the typed run configuration objects;
- the setup types and helpers (`SimulationConfiguration`, `SpacecraftModel`,
  `Link`, `InitialTime`, `make_example_config` and others), so `using SpaceAGORA`
  alone is enough to configure a run;
- the no-GRAM baseline constructors and builders;
- the abstract interfaces and hook functions for extensions (see
  [Extension hooks](#Extension-hooks) below);
- the verification, CLI and asset-report entry points.

Three versioned facades rooted at `SpaceAGORA` are also supported boundaries:
`HYPRServices`, `ExternalPropagationServices` and `EDGServices`. Use them by
qualified name, for example `SpaceAGORA.EDGServices`. Their documented
metadata, lifecycle operations and explicit binding inventories make up the
supported boundary, and their [adapter contract](reference/adapter_contract.md)
versions them separately from the package version. Treat any name that is
neither in the documented root API nor in a facade inventory as internal.

## Compliant multibody and robot-arm names

Use the robot-arm names for kinematics and the compliant multibody names for
spring-connected dynamics. The [Robotics API](generated/robotics_api.md) and
[Compliant Multibody API](generated/compliant_multibody_api.md) document both
the preferred names and the historical ones. The model represents flexibility
with discrete bodies and compliant joints; it is not a continuum cloth solver.
The preferred names are aliases of the same types and generic functions, so
changing how a call is spelled changes neither its dispatch, nor its state
layout, nor its physics.

| Historical name | Preferred name |
| --- | --- |
| `ClothArmBasePose`, `ClothArmLink`, `ClothArmJoint`, `ClothArmModel`, `ClothArmPose`, `ClothArmState` | `RobotArmBasePose`, `RobotArmLink`, `RobotArmJoint`, `RobotArmModel`, `RobotArmPose`, `RobotArmState` |
| `default_cloth_arm_model` | `default_robot_arm_model` |
| `cloth_fk`, `cloth_fk_state`, `cloth_ik`, `cloth_total_reach`, `cloth_end_effector_pose` | the corresponding `robot_arm_` names |
| `ClothRobotArmSimulation`, `ClothRobotArmReferenceState` | `CompliantRobotArmSimulation`, `CompliantRobotArmReferenceState` |
| `cloth_reference_state` | `compliant_robot_arm_reference_state` |
| functions containing `cloth_robot_arm`, such as `simulate_cloth_robot_arm_plan` | the same names with `compliant_robot_arm`, such as `simulate_compliant_robot_arm_plan` |

For example, write `arm = default_robot_arm_model()` and
`pose = robot_arm_fk(arm, RobotArmBasePose([0.0, 0.0, 0.0]), q)`, where `q`
holds the arm's joint angles. Calls that use the historical names keep working.
They are still exported, except `ClothRobotArmReferenceState` and
`cloth_reference_state`, which you reach as `SpaceAGORA.ClothRobotArmReferenceState`
and `SpaceAGORA.cloth_reference_state`. The defining modules and types keep
their identities, so printed types and serialized owner names can still contain
`Cloth`. The root also provides the module aliases
`SpaceAGORA.CompliantMultibody` and `SpaceAGORA.CompliantRobotArmDynamics`;
they don't make every member of those modules public, and the supported
operations remain the ones in the public reference.

Some historical spellings are configuration or data identifiers rather than
display names, and they stay: `planner=:cloth_quintic`, the `retime_cloth_*`
configuration fields, and the existing example file names, output directories
and CSV column names. The `HYPRServices.Cloth` namespace of the versioned HYPR
facade and its inventory keep their exact spellings; see the
[HYPR service contract](reference/hypr_services.md#Cloth). None of this adds
contact physics or couples the core compliant model to MuJoCo; the optional
[MuJoCo companion](user/mujoco_scenes.md) has its own propagation ownership and
coupling limits.

## Internal modules

These modules are internal implementation namespaces:

- `RuntimeServices`
- `SimulationModel`
- `SimulationEngine`
- `SimulationCampaigns`
- `ParallelProfiles`
- `ParallelProcess`
- `TelemetryVerification`
- `SpaceAGORACLI`
- the subsystem modules under `src/models/gnc`, `src/models/dynamics`,
  `src/orchestration/simulation` and `src/orchestration/parallel`

The package and its tests use them heavily, but they are not the user contract.
They can be reorganized without a deprecation window, as long as the root
exported API stays stable.

## Extension hooks

A model author depends on two things: the documented abstract interfaces
exported by `SpaceAGORA`, which a new model subtypes, and the exported hook
functions it implements, which follow the
[Extension Conventions](reference/extension_conventions.md). By family:

- force and torque models: `wrench`, `environment_requirements`,
  `solver_partition` and the gravity-backbone hooks;
- atmospheres: `getDensity` and `getDensityBatch!`;
- GNC ticks: `update_control!`, `update_guidance!` and `update_navigation!`;
- control outputs: `control_wrench`, `control_propellant_rate` and
  `control_wheel_torque`;
- actuators: `actuator_wrench`, `actuator_limits` and the optional actuator
  methods `actuation_mode` and `actuator_levels`;
- sensors: `sensor_measure`, `sensor_dimension` and the optional sensor methods;
- planners: `plan_rpo!`, `retime_rpo!` and the GNC algorithm contract.

The [Extension Hook Reference](generated/extension_hook_reference.md) lists
every hook with its shape, signature, fallback and status. Extension code can
therefore stay on `SpaceAGORA.*` names instead of reaching into internal
modules.

Two of these hooks concern the solver. `solver_partition` puts a dynamic
effector on the implicit or the explicit side of `solver_mode=:split_imex`; it
does not select the mode, and it does not move a custom effector into the fast
part of `:multirate_mri`. `gravity_backbone_structure` and
`gravity_backbone_acceleration_ii` declare an effector's part in the symplectic
gravity core of `:gravity_backbone_split`, and `gravity_backbone_kick_structure`
and `gravity_backbone_kick_acceleration_ii` its part in the explicit velocity
kicks, such as solar radiation pressure and n-body gravity.
[Solver Configuration](user/solver_configuration.md) describes the modes.

## Deprecations

A public name that is renamed or retired keeps working for at least one tagged
release:

- The old name still works, and its docstring opens with "Deprecated, use ...".
  The [Extension Hook Reference](generated/extension_hook_reference.md) gives a
  deprecated hook's replacement.
- Using the old name logs a warning once per owner and name.
- `julia --depwarn=error` turns these warnings into errors, so a test run can
  refuse deprecated names.
- The name is removed in a later release, once a tagged release has carried the
  warning and the companion packages have moved. SpaceAGORA 0.5.0, for example,
  removed the hook names that 0.3.0 deprecated.
- Release numbers follow `version` in `Project.toml`.

The [deprecation rules](reference/extension_conventions.md#Deprecations) give
the details, and [Migrating Extension Hooks](reference/migrating_extension_hooks.md)
maps every old hook name to its replacement, with the schedule.
[Upgrading](user/upgrading.md) lists what each release asks you to change.

Internal modules are not covered by this window; they can change at any time, as
long as the root exported API stays stable.

## Opt-in RPO planner surface

The root exports the neutral planner contracts and extension methods,
`DirectRPOPlanner`, `HYPRRPOPlanner`, `RetainedPathRPOPlanner`, `RPOPSOConfig`,
`rpo_pso_config`, `RPOPlanningEvent`, `RPOPlanningError`,
`make_rpo_configuration` and `rpo_run_report`, and it also exports the output
settings type `SimulationSettings`. The [RPO Planner API](generated/rpo_planner_api.md)
lists them, and the [pilot guide](user/rpo_planner_pilot.md) defines the
supported truth and static-RTN scope, the failure policy and the restart
refusal. The internal lifecycle modules and the full `SimulationModel` are not
public interfaces. The root also exports the seven comparison planners,
`HYPRComparisonPlanner`, `UnrefinedPSOPlanner`, `RRTConnectPlanner`,
`RRTConnectBezierPlanner`, `RRTStarPlanner`, `CHOMPPlanner` and `STOMPPlanner`;
the comparison harness itself is internal.
