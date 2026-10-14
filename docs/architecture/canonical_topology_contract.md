# SpaceAGORA Canonical Topology Contract

This contract defines where SpaceAGORA implementation responsibilities live.

The source tree has four groups: `src/interfaces/` contains `cli`, `io` and
`analysis`; `src/orchestration/` contains `simulation`, `mission` and `parallel`;
`src/models/` contains `vehicle`, `environment`, `dynamics` and `gnc`; and
`src/foundation/` contains `core` and `assets`. These groups organize files while
preserving the existing Julia modules and their include order. They do not impose
a new dependency hierarchy: `orchestration/simulation/composition/simulation_model.jl` remains the
composition loader for definitions owned across the groups.

`src/SpaceAGORA.jl` and `src/precompile_workload.jl` remain at the package source
root. Julia discovers package extensions in repository-root `ext/`. Optional
companions retain their own repositories and source trees outside core `src/`.

## Locked Ownership Decisions
1. Benchmarks and studies are owned only by top-level `benchmarks/` (never under `src/`).
2. Runnable examples are owned only by top-level `examples/` (never under `src/`).
3. Operational plotting scripts are owned only by top-level `scripts/plotting/` (never under `src/interfaces/analysis/plotting/`).
4. Domain dynamics ownership is:
   - `src/models/dynamics/translational/*`
   - `src/models/dynamics/rotational/*`
   - `src/models/dynamics/coupled/*`
   and `src/models/dynamics/models/` is retired.
5. Mission configuration ownership: initial-state definitions live with the
   spacecraft model (`src/models/vehicle/spacecraft/model.jl`); plans, schedules, and
   aerobraking policy live under `src/orchestration/mission/operations/*`.
   One-run scenario configuration lives under `src/orchestration/simulation/config/`; the existing
   spacecraft collection and selected effectors are defined in
   `constellation_configuration.jl`, retaining the `SpacecraftModels.DynamicsModel`
   identity. Physical spacecraft assembly remains under vehicle.
6. Simulation execution ownership belongs to `src/orchestration/simulation/engine/*`; legacy `src/orchestration/simulation/events/*`, `src/orchestration/simulation/execution/*`, and `src/orchestration/simulation/solver_orchestration/*` are retired.
7. Shared runtime serialization services are owned by `src/orchestration/simulation/runtime_services.jl`, not by `SimulationModel`.
8. Actuator ownership: thruster-specific hooks live in
   `src/models/vehicle/actuators/thruster/thruster_hooks.jl`.
9. Vehicle boundary is:
   - `src/models/vehicle/spacecraft/*` for composition/integration
   - `src/models/vehicle/structure/*` for mass/inertia/geometry/structure analysis,
     including the CAD/mesh readers (`mesh_geometry.jl`) shared by the
     visualization scene
10. Visualization boundary is:
   - `src/interfaces/analysis/visualization/scene/*` owns the renderer-independent scene
     layer (scene types, spacecraft geometry conversion, planet rotation
     sampling, the JSON sidecar) as a `SimulationModel` submodule
   - `src/interfaces/analysis/visualization/rpo/*` owns the PlotlyJS RPO plots
   - the browser viewer and HTML assembly are owned by the independent
     `SpaceAGORAVisualizer.jl` package; core keeps portable payload export and
     `ext/SpaceAGORAVisualizerExt.jl` connects the optional renderer. No renderer
     source is duplicated under `src/viewer/` or top-level `viewer/`.

11. Externally propagated spacecraft (shadow entries held in the state vector while an
    optional package integrates them) are owned in four places, none of which
    knows any concrete integrator:
    - `src/models/dynamics/coupled/external_propagation.jl`: the `SimulationModel.ExternalPropagation`
      hook implementation identities. Optional adapters use the supported facade below;
      default hooks throw `Not implemented`, except for the accepting preflight default.
    - `src/orchestration/simulation/engine/external_propagation.jl`: setup-time refusals, per-run runtimes
      and the shadow right-hand side.
    - `src/orchestration/simulation/callbacks/external_propagation_callbacks.jl`: the sync callback that
      overwrites shadow entries after every accepted step, first among the discrete callbacks.
    - `src/orchestration/simulation/external_propagation_services.jl`: the versioned
      `SpaceAGORA.ExternalPropagationServices` adapter boundary, aliasing existing hook
      identities and supplying narrow configuration, state and frame services.
    `SimulationConfiguration.external_propagators` carries the owners. Concrete owners live
    in their domain repositories, including the independent
    [SpaceAGORAMuJoCo.jl companion](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl).
    The root package never depends on them. Native bindings, adapter tests and the
    Basilisk/MuJoCo comparison harness belong to that companion; generic lifecycle,
    conformance tests and shared physics remain in core.

## Required Canonical Files
1. `src/orchestration/simulation/engine/public_api.jl`
2. `src/orchestration/simulation/engine/execution.jl`
3. `src/orchestration/simulation/runtime_services.jl`
4. `src/models/dynamics/coupled/force_torque_models.jl`
5. `src/models/vehicle/actuators/thruster/thruster_hooks.jl`
6. `src/models/vehicle/structure/structure_models.jl`
7. `src/models/vehicle/structure/assembly_graph.jl`
8. `src/models/vehicle/structure/mass_properties.jl`
8a. `src/models/vehicle/structure/mesh_geometry.jl`
9. `src/models/vehicle/structure/geometry_properties.jl`

## Path Policy
1. This contract documents canonical owners only.
2. Removed or forbidden legacy paths are intentionally not listed here as a
   tombstone inventory.
3. Canonical-path enforcement belongs to CI gates, not to user-facing topology
   documentation.

## Canonical GNC and Runtime IO Owners
1. GNC hook owners: `src/models/gnc/control/control_hooks.jl`, `src/models/gnc/guidance/guidance_hooks.jl`, `src/models/gnc/navigation/navigation_hooks.jl`.
2. GNC bridge helper owner: `src/models/gnc/internal/bridge_helpers.jl`.
3. Runtime IO owners: `src/interfaces/io/paths/io_paths.jl`, `src/interfaces/io/checkpoints/checkpoint_storage.jl`, `src/interfaces/io/outputs/io_outputs.jl`; `src/orchestration/simulation/engine/*` orchestrates IO but does not own serialization or output implementations.
4. No `legacy_` identifiers in `src/`.

## Enforcement Gates
1. `test/gates/ci_path_policy_gate.jl` (retired paths, required owner files, forbidden path tokens)
2. `test/gates/ci_vehicle_structure_boundary_gate.jl`
3. `test/gates/ci_architecture_contract_gate.jl`
4. `test/gates/ci_thin_entry_files_gate.jl` (benchmark launchers stay thin forwarders)

## Shared Configuration and Mathematics

- `SimulationModel.SimConfig._with_configuration` copies an existing simulation
  configuration and applies named overrides. Unchanged fields retain their
  references, and model types are inferred again by the constructor. Callers
  retain responsibility for nested settings and intentional resets. The telemetry
  runner explicitly clears the typed solver so its environment-based retry
  policy remains authoritative. Engine state isolation remains separate.
- `SimulationModel.QuaternionMath` owns the quaternion functions in
  `src/foundation/core/numerics/quaternion_utils.jl`. Consumers import this module rather
  than including another copy of the implementation. Quaternions are scalar-last;
  `rot(q)` maps inertial components to body components.
- `SimulationModel.FrameTransforms` owns the coordinate functions in
  `src/foundation/core/interfaces/reference_system.jl`. Engine, callbacks, GNC, aerodynamics
  and telemetry import these same bindings. `rotate_vector_by_quaternion` is an
  active vector rotation, so it is not interchangeable with `rot(q)`.
- These are internal owners, not new root-level public APIs. Existing qualified
  consumer bindings remain available. `ReferenceSystems` continues to own the
  reference-system types.
