# SpaceAGORA Canonical Topology Contract

This contract defines canonical ownership for the topology cleanup that answers the six architecture placement questions.

## Locked Ownership Decisions
1. Benchmarks and studies are owned only by top-level `benchmarks/` (never under `src/`).
2. Runnable examples are owned only by top-level `examples/` (never under `src/`).
3. Operational plotting scripts are owned only by top-level `scripts/plotting/` (never under `src/analysis/plotting/`).
4. Domain dynamics ownership is:
   - `src/dynamics/translational/*`
   - `src/dynamics/rotational/*`
   - `src/dynamics/coupled/*`
   and `src/dynamics/models/` is retired.
5. Mission configuration ownership is split between:
   - `src/mission/initial_conditions.jl` for canonical initial-state definitions
   - `src/mission/operations/*` for plans, schedules, and aerobraking policy ownership
6. Simulation execution dispatch belongs to `src/simulation/execution/*`; runtime internals stay in `src/simulation/engine/*`.
7. Actuator ownership is split between:
   - generic hook surface: `src/vehicle/actuators/actuator_hooks.jl`
   - thruster-specific hooks: `src/vehicle/actuators/thruster/thruster_hooks.jl`
8. Vehicle boundary is:
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
2. `src/examples/`
3. `src/analysis/plotting/`
4. `src/analysis/reports/`
5. `src/parallel/state/`
6. `src/parallel/tuning/`
7. `src/vehicle/sensors/`
8. `src/core/utils/typed_example_utils.jl`
9. `src/dynamics/models/`
10. `src/mission/campaigns/run.jl`
11. `src/mission/campaigns/set_and_run.jl`
12. `src/mission/campaigns/define_mission.jl`
13. `src/mission/campaigns/mission_model.jl`
14. `src/mission/campaigns/initial_cond_calc.jl`
15. `src/vehicle/actuators/thruster_effectors.jl`
16. `src/vehicle/spacecraft/spacecraft_analysis.jl`

## Enforcement Gates
1. `test/ci_no_src_benchmarks_root_gate.jl`
2. `test/ci_no_dynamics_models_gate.jl`
3. `test/ci_vehicle_structure_boundary_gate.jl`
4. `test/ci_canonical_path_contract_gate.jl`
5. `test/ci_architecture_contract_gate.jl`
