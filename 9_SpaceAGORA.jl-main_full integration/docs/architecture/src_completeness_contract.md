# SpaceAGORA `src/` Completeness Contract

This contract defines how canonical source files are classified and what is forbidden in canonical code.

## File Categories
1. `Canonical Owner`
   - Owns runtime behavior.
   - Must contain executable logic, typed data structures, or concrete API methods.
2. `Canonical Aggregator`
   - Include/re-export composition only.
   - No behavior ownership.
   - Must start with header comment:
     - `# Canonical aggregator: no behavior ownership.`
   - Allowlisted paths:
     - `src/orchestration/simulation/engine/simulation_engine.jl`
     - `src/orchestration/simulation/callbacks/callbacks.jl`
     - `src/orchestration/parallel/routing/parallel_profiles.jl`
     - `src/orchestration/parallel/process/parallel_process.jl`
     - `src/orchestration/simulation/composition/simulation_model.jl`
3. `Experimental Interface`
   - Typed contracts for planned subsystems.
   - Lives under top-level `experimental/*`, never under `src/*`.
   - Unimplemented behavior fails explicitly with:
     - `throw(ErrorException("Not implemented: <method>"))`
4. `Forbidden`
   - Silent stubs such as `*_stub() = nothing`
   - Canonical files that only forward to legacy source ownership
   - Computed-but-unused force terms in active dynamics paths

## Translational Boundary
1. `src/models/dynamics/translational/*` is reserved for pure translational terms only.
2. Gravity field kernels live under `src/models/environment/gravity/`; spacecraft/runtime
   adapters live under `src/models/dynamics/coupled/gravity/`. The compatibility include
   manifest preserves the existing `DynamicEffectors.GravityEffectors` module ownership
   and declaration order. Neither responsibility belongs under translational dynamics.
3. Translational files must not own attitude-coupled torque logic or actuator scheduling/commands.

## Topology Boundary
1. Benchmark ownership is top-level `benchmarks/*`; `src/benchmarks/*` is forbidden.
2. Runnable example ownership is top-level `examples/*`; `src/examples/*` is forbidden.
3. Plotting/report tooling ownership is top-level `scripts/*`; `src/interfaces/analysis/plotting/*` is forbidden.
4. Unfinished subsystem scaffolds belong under top-level `experimental/*`; mirrored `src/*` copies are forbidden.
5. Only the canonical roots described in this contract are valid ownership
   locations for plotting, reporting, and parallel support surfaces.
6. Coupled dynamics owners live under `src/models/dynamics/coupled/*`.
7. Simulation execution ownership is `src/orchestration/simulation/engine/public_api.jl` and `src/orchestration/simulation/engine/execution.jl`; `src/orchestration/simulation/events/*`, `src/orchestration/simulation/execution/*`, and `src/orchestration/simulation/solver_orchestration/*` are forbidden.
8. Mission configuration ownership: initial-state definitions live with the spacecraft model (`src/models/vehicle/spacecraft/model.jl`); mission plans and policies live under `src/orchestration/mission/operations/*`.
9. Structural analysis ownership is `src/models/vehicle/structure/*`, not `src/models/vehicle/spacecraft/*`; CAD/mesh readers live in `src/models/vehicle/structure/mesh_geometry.jl`.
10. Visualization scene ownership is `src/interfaces/analysis/visualization/scene/*` (a canonical owner inside `SimulationModel`); portable payload export stays with that scene owner. The browser viewer lives in independent `SpaceAGORAVisualizer.jl`; `ext/SpaceAGORAVisualizerExt.jl` connects it, and duplicated renderer source under `src/viewer/*` or top-level `viewer/*` is forbidden.

## GNC Aerobraking Boundary
1. `E-EDG` and `T-EDG` strategy ownership is `src/models/gnc/guidance/aerobraking/*`.
2. `src/models/gnc/control/aerobraking/*` is tracking/execution only and must not define strategy solvers.
3. Guidance files must not include control source files directly.
4. Strategy selection contracts are mission-owned in `src/orchestration/mission/operations/aerobraking_policy/*`.
5. Propulsive maneuver guidance/control coupling must use typed command buffers; guidance must not mutate thruster runtime state directly.
6. The typed EDG numerical kernels live in the optional `SpaceAGORAEDG.jl` package.
   Core `typed_edg/algorithms.jl` is a compatibility dispatcher, not an include
   aggregator: it must contain exactly the declared package forwarders and no
   numerical implementations. The source-completeness and EDG ownership gates
   enforce that exact inventory. Core physical services include only their
   reviewed `adapter_queries.jl` capability wrapper; the companion receives no
   raw integrator state or mutable vehicle geometry.

## Experimental Interface Set
1. `experimental/vehicle/resources/resources.jl`
2. `experimental/vehicle/resources/battery/battery_model.jl`
3. `experimental/vehicle/resources/solar_array/solar_array_model.jl`
4. `experimental/vehicle/resources/power_bus/power_bus_model.jl`
5. `experimental/vehicle/resources/loads/load_model.jl`
6. `experimental/mission/constellation/network/network.jl`
7. `experimental/gnc/estimation/estimation.jl`
8. `experimental/vehicle/laser_terminal/laser_terminal.jl`

## Enforced Gates
1. `test/gates/ci_src_completeness_contract_gate.jl`
2. `test/gates/ci_no_legacy_include_chains_gate.jl`
3. `test/gates/ci_no_inert_force_terms_gate.jl`
4. `test/gates/ci_scaffold_interface_gate.jl`
5. `test/gates/ci_translational_ownership_gate.jl`
6. `test/gates/ci_path_policy_gate.jl`
7. `test/gates/ci_vehicle_structure_boundary_gate.jl`
8. `test/gates/ci_gnc_aerobraking_boundary_gate.jl`
9. `test/gates/ci_no_guidance_control_cross_include_gate.jl`
10. `test/gates/ci_gnc_typed_command_boundary_gate.jl`
11. `test/gates/ci_aerobraking_strategy_contract_gate.jl`
