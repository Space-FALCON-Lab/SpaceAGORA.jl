# SpaceAGORA Canonical Topology Contract

This contract defines canonical ownership for the topology cleanup that answers the six architecture placement questions.

## Locked Ownership Decisions
1. Benchmarks and studies are owned only by top-level `benchmarks/` (never under `src/`).
2. Domain dynamics ownership is:
   - `src/dynamics/translational/*`
   - `src/dynamics/rotational/*`
   - `src/dynamics/coupled/*`
   and `src/dynamics/models/` is retired.
3. Mission definition owners are at mission root:
   - `src/mission/define_mission.jl`
   - `src/mission/mission_model.jl`
   - `src/mission/initial_conditions.jl`
4. Simulation execution dispatch belongs to `src/simulation/execution/*`; runtime internals stay in `src/simulation/engine/*`.
5. Actuator ownership is split between:
   - generic hook surface: `src/vehicle/actuators/actuator_hooks.jl`
   - thruster-specific hooks: `src/vehicle/actuators/thruster/thruster_hooks.jl`
6. Vehicle boundary is:
   - `src/vehicle/spacecraft/*` for composition/integration
   - `src/vehicle/structure/*` for mass/inertia/geometry/structure analysis

## Required Canonical Files
1. `src/simulation/execution/run.jl`
2. `src/simulation/execution/dispatch.jl`
3. `src/dynamics/coupled/force_torque_models.jl`
4. `src/vehicle/actuators/actuator_hooks.jl`
5. `src/vehicle/actuators/thruster/thruster_hooks.jl`
6. `src/vehicle/structure/structure_models.jl`
7. `src/vehicle/structure/assembly_graph.jl`
8. `src/vehicle/structure/mass_properties.jl`
9. `src/vehicle/structure/geometry_properties.jl`

## Retired Paths
1. `src/benchmarks/`
2. `src/dynamics/models/`
3. `src/mission/campaigns/run.jl`
4. `src/mission/campaigns/set_and_run.jl`
5. `src/mission/campaigns/define_mission.jl`
6. `src/mission/campaigns/mission_model.jl`
7. `src/mission/campaigns/initial_cond_calc.jl`
8. `src/vehicle/actuators/thruster_effectors.jl`
9. `src/vehicle/spacecraft/spacecraft_analysis.jl`

## Enforcement Gates
1. `test/ci_no_src_benchmarks_root_gate.jl`
2. `test/ci_no_dynamics_models_gate.jl`
3. `test/ci_vehicle_structure_boundary_gate.jl`
4. `test/ci_canonical_path_contract_gate.jl`
5. `test/ci_architecture_contract_gate.jl`
