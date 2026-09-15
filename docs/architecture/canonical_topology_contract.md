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
5. Mission configuration ownership: initial-state definitions live with the
   spacecraft model (`src/vehicle/spacecraft/model.jl`); plans, schedules, and
   aerobraking policy live under `src/mission/operations/*`.
6. Simulation execution ownership belongs to `src/simulation/engine/*`; legacy `src/simulation/events/*`, `src/simulation/execution/*`, and `src/simulation/solver_orchestration/*` are retired.
7. Shared runtime serialization services are owned by `src/simulation/runtime_services.jl`, not by `SimulationModel`.
8. Actuator ownership: thruster-specific hooks live in
   `src/vehicle/actuators/thruster/thruster_hooks.jl`.
9. Vehicle boundary is:
   - `src/vehicle/spacecraft/*` for composition/integration
   - `src/vehicle/structure/*` for mass/inertia/geometry/structure analysis,
     including the CAD/mesh readers (`mesh_geometry.jl`) shared by the
     visualization scene and the aerodynamic panel method
10. Visualization boundary is:
   - `src/analysis/visualization/scene/*` owns the renderer-independent scene
     layer (scene types, spacecraft geometry conversion, planet rotation
     sampling, the JSON sidecar) as a `SimulationModel` submodule
   - `src/analysis/visualization/rpo/*` owns the PlotlyJS RPO plots
   - the browser viewer itself is owned only by top-level `viewer/` (never
     under `src/viewer/`)

## Required Canonical Files
1. `src/simulation/engine/public_api.jl`
2. `src/simulation/engine/execution.jl`
3. `src/simulation/runtime_services.jl`
4. `src/dynamics/coupled/force_torque_models.jl`
5. `src/vehicle/actuators/thruster/thruster_hooks.jl`
6. `src/vehicle/structure/structure_models.jl`
7. `src/vehicle/structure/assembly_graph.jl`
8. `src/vehicle/structure/mass_properties.jl`
8a. `src/vehicle/structure/mesh_geometry.jl`
8b. `src/dynamics/coupled/aerodynamic_mesh_surrogate.jl`
9. `src/vehicle/structure/geometry_properties.jl`

## Path Policy
1. This contract documents canonical owners only.
2. Removed or forbidden legacy paths are intentionally not listed here as a
   tombstone inventory.
3. Canonical-path enforcement belongs to CI gates, not to user-facing topology
   documentation.

## Canonical GNC and Runtime IO Owners
1. GNC hook owners: `src/gnc/control/control_hooks.jl`, `src/gnc/guidance/guidance_hooks.jl`, `src/gnc/navigation/navigation_hooks.jl`.
2. GNC bridge helper owner: `src/gnc/internal/bridge_helpers.jl`.
3. Runtime IO owners: `src/io/config/io_config.jl`, `src/io/serialization/io_serialization.jl`, `src/io/outputs/io_outputs.jl`; `src/simulation/engine/*` orchestrates IO but does not own serialization or output implementations.
4. No `legacy_` identifiers in `src/`.

## Enforcement Gates
1. `test/gates/ci_path_policy_gate.jl` (retired paths, required owner files, forbidden path tokens)
2. `test/gates/ci_vehicle_structure_boundary_gate.jl`
3. `test/gates/ci_architecture_contract_gate.jl`
4. `test/gates/ci_thin_entry_files_gate.jl` (benchmark launchers stay thin forwarders)
