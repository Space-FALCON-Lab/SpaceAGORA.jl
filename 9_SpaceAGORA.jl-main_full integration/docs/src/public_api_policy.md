# Public API Policy

Use this page when you need to know which package surface is stable and which
parts of the repository can still move.

This page is for maintainers and extension authors. If you only want to run
SpaceAGORA, use the User Guide and the generated Public API page instead.

Shortest stable package example:

```julia
using SpaceAGORA
```

What to read next:

- [Public API](generated/public_api.md)
- [Extensibility](extensibility.md)
- [Maintainer Overview](maintainer/index.md)

SpaceAGORA now treats the root `SpaceAGORA` module as the only stable package
surface.

## Stable surface

The generated [Public API](generated/public_api.md) reference and its linked
[Simulation Setup API](generated/simulation_setup_api.md),
[Terrain API](generated/terrain_api.md),
[Mesh Aerodynamics API](generated/mesh_aerodynamics_api.md) and
[Visualization API](generated/visualization_api.md) pages document the supported
exported symbols. That surface currently includes:

- simulation entrypoints such as `run_simulation`
- typed runtime configuration objects
- simulation setup types and helpers (`SimulationConfiguration`, `SpacecraftModel`,
  `Link`, `InitialTime`, `make_example_config`, ...), so `using SpaceAGORA` alone
  is enough to configure a run
- no-GRAM baseline constructors and builders
- selected abstract extension interfaces
- selected extension hook functions
- verification entrypoints
- CLI and asset-report entrypoints

Documented versioned facades rooted at `SpaceAGORA` are also supported integration
boundaries: `HYPRServices` and `ExternalPropagationServices`. Their documented
metadata, lifecycle operations and explicit binding inventories define the
supported boundary; use them by qualified name. Their
[adapter contract](maintainer/adapter_contract.md) controls compatibility
separately from package versions. Other symbols absent from the documented root
API or a facade inventory should be treated as internal.

## Compliant multibody and robot-arm names

Use the robot-arm names for kinematics and the compliant multibody names for
spring-connected dynamics. The [Robotics API](generated/robotics_api.md) and
[Compliant Multibody API](generated/compliant_multibody_api.md) document both the
preferred names and the compatibility names. The model represents flexibility through discrete
bodies and compliant joints; it is not a continuum cloth solver. The preferred
names are aliases of the existing types and generic functions, so changing a
call's spelling does not change its dispatch, state layout or physics.

| Historical name | Preferred name |
| --- | --- |
| `ClothArmBasePose`, `ClothArmLink`, `ClothArmJoint`, `ClothArmModel`, `ClothArmPose`, `ClothArmState` | `RobotArmBasePose`, `RobotArmLink`, `RobotArmJoint`, `RobotArmModel`, `RobotArmPose`, `RobotArmState` |
| `default_cloth_arm_model` | `default_robot_arm_model` |
| `cloth_fk`, `cloth_fk_state`, `cloth_ik`, `cloth_total_reach`, `cloth_end_effector_pose` | Corresponding `robot_arm_` names |
| `ClothRobotArmSimulation`, `ClothRobotArmReferenceState` | `CompliantRobotArmSimulation`, `CompliantRobotArmReferenceState` |
| `cloth_reference_state` | `compliant_robot_arm_reference_state` |
| Functions containing `cloth_robot_arm`, including `simulate_cloth_robot_arm_plan` | The same names with `compliant_robot_arm`, including `simulate_compliant_robot_arm_plan` |

For example, write `arm = default_robot_arm_model()` and
`pose = robot_arm_fk(arm, RobotArmBasePose([0.0, 0.0, 0.0]), q)` when `q` contains
the arm's joint angles. Existing calls using the historical names remain valid.
The defining modules and types retain their identities, so printed types and
serialized owner names may still contain `Cloth`.
The root also provides the qualified module aliases
`SpaceAGORA.CompliantMultibody` and `SpaceAGORA.CompliantRobotArmDynamics`.
These aliases do not make every module member public; supported operations
remain those listed in the public reference.

Some historical strings are configuration or data identifiers, rather than
display names. Keep `planner=:cloth_quintic`, `mode=:cloth`, the `retime_cloth_*`
configuration fields, and existing example filenames, output directories and
CSV column names. The versioned `HYPRServices.Cloth` namespaces and their
inventory retain their exact spellings; see the [HYPR service contract](maintainer/hypr_services.md#Cloth).
This compatibility change neither adds contact physics nor couples the core
compliant model to MuJoCo. The optional [MuJoCo companion](user/mujoco_scenes.md)
has its own propagation ownership and coupling limits.

## Internal modules

The following modules remain internal implementation detail namespaces:

- `RuntimeServices`
- `SimulationModel`
- `SimulationEngine`
- `ParallelProfiles`
- `TelemetryVerification`
- `SpaceAGORACLI`
- subsystem modules under `src/models/gnc`, `src/models/dynamics`, `src/orchestration/simulation`, and
  `src/orchestration/parallel`

These modules are still used heavily inside the package and tests, but they are
not the stable user contract. They may be reorganized without a deprecation
window as long as the root exported API remains stable.

## Extension hooks

For model authors, the stable extension contract is:

- subtype the documented abstract interfaces exported by `SpaceAGORA`
- implement the matching exported hook functions such as:
  - `wrench`
  - `environment_requirements`
  - `solver_partition`
  - `gravity_backbone_structure`
  - `gravity_backbone_acceleration_ii`
  - `calcForceTorque`
  - `getDensity`
  - `getDensityBatch!`
  - `calcControlEffect!`
  - `calcControlForceTorque`
  - `calcControlMassFlowRate`

This lets extension code stay rooted in `SpaceAGORA.*` rather than reaching into
internal modules.

Solver-side policy note:

- `split_imex` now means the atmosphere-implicit IMEX partition
- `multirate` still means the control-focused split
- `solver_partition` is the public hook that places dynamic effectors on the implicit or explicit IMEX side
- `gravity_backbone_split` is the public fixed-step gravity-backbone split mode with a symplectic gravity core plus explicit SRP / N-body velocity kicks
- `gravity_backbone_structure` / `gravity_backbone_acceleration_ii` declare gravity-core participation
- `gravity_backbone_kick_structure` / `gravity_backbone_kick_acceleration_ii` declare explicit translational kick participation
- the current `heat_loads` state remains explicit because it is an accumulated heat-rate integral, not a thermal-capacitance model

## Documentation enforcement

The docs build checks the root exported public API list and fails if any symbol
listed in the generated public API page is missing documentation.

The generated undocumented-export report is written to:

```text
docs/build/undocumented_public_exports.txt
```

That report is intentionally scoped to the stable root API, not to internal
modules.

## Opt-in RPO planner surface

The root exports the neutral planner contracts and extension methods,
`DirectRPOPlanner`, `HYPRRPOPlanner`, `RPOPSOConfig`, `rpo_pso_config`,
`RPOPlanningEvent`, `RPOPlanningError`, `make_rpo_configuration` and
`rpo_run_report`. `SimulationSettings` remains owned by the existing configuration
module and has a root alias for output settings. The [RPO Planner API](generated/rpo_planner_api.md) lists these symbols. The [pilot guide](user/rpo_planner_pilot.md)
defines supported truth/static-RTN scope, failure policy and restart refusal.
Internal lifecycle modules and the full `SimulationModel` are not public interfaces.
