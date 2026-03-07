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
     - `src/simulation/engine/simulation_engine.jl`
     - `src/simulation/callbacks/callbacks.jl`
     - `src/parallel/routing/parallel_profiles.jl`
     - `src/core/simulation_model.jl`
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
1. `src/dynamics/translational/*` is reserved for pure translational terms only.
2. Gravity ownership is canonical under `src/environment/gravity/gravity_models.jl` (not under translational).
3. Translational files must not own attitude-coupled torque logic or actuator scheduling/commands.

## Topology Boundary
1. Benchmark ownership is top-level `benchmarks/*`; `src/benchmarks/*` is forbidden.
2. Runnable example ownership is top-level `examples/*`; `src/examples/*` is forbidden.
3. Plotting/report tooling ownership is top-level `scripts/*`; `src/analysis/plotting/*` is forbidden.
4. Unfinished subsystem scaffolds belong under top-level `experimental/*`; mirrored `src/*` copies are forbidden.
5. Only the canonical roots described in this contract are valid ownership
   locations for plotting, reporting, and parallel support surfaces.
6. Coupled dynamics owners live under `src/dynamics/coupled/*`.
7. Simulation execution ownership is `src/simulation/engine/public_api.jl` and `src/simulation/engine/execution.jl`; `src/simulation/events/*`, `src/simulation/execution/*`, and `src/simulation/solver_orchestration/*` are forbidden.
8. Mission configuration ownership is `src/mission/initial_conditions.jl`; mission plans and policies live under `src/mission/operations/*`.
9. Structural analysis ownership is `src/vehicle/structure/*`, not `src/vehicle/spacecraft/*`.

## GNC Aerobraking Boundary
1. `E-EDG` and `T-EDG` strategy ownership is `src/gnc/guidance/aerobraking/*`.
2. `src/gnc/control/aerobraking/*` is tracking/execution only and must not define strategy solvers.
3. Guidance files must not include control source files directly.
4. Strategy selection contracts are mission-owned in `src/mission/operations/aerobraking_policy/*`.
5. Propulsive maneuver guidance/control coupling must use typed command buffers; guidance must not mutate thruster runtime state directly.

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
1. `test/ci_src_completeness_contract_gate.jl`
2. `test/ci_no_legacy_include_chains_gate.jl`
3. `test/ci_no_inert_force_terms_gate.jl`
4. `test/ci_scaffold_interface_gate.jl`
5. `test/ci_translational_ownership_gate.jl`
6. `test/ci_no_src_benchmarks_root_gate.jl`
7. `test/ci_no_dynamics_models_gate.jl`
8. `test/ci_vehicle_structure_boundary_gate.jl`
9. `test/ci_gnc_aerobraking_boundary_gate.jl`
10. `test/ci_no_guidance_control_cross_include_gate.jl`
11. `test/ci_gnc_typed_command_boundary_gate.jl`
12. `test/ci_aerobraking_strategy_contract_gate.jl`
