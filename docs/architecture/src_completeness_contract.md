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
3. `Scaffold Interface`
   - Typed contracts for planned subsystems.
   - Unimplemented behavior must fail explicitly with:
     - `throw(ErrorException("Not implemented: <method>"))`
4. `Forbidden`
   - Silent stubs such as `*_stub() = nothing`
   - Canonical files that only forward to legacy source ownership
   - Computed-but-unused force terms in active dynamics paths

## Translational Boundary
1. `src/dynamics/translational/*` is reserved for pure translational terms only.
2. Gravity ownership is canonical under `src/environment/gravity/gravity_models.jl` (not under translational).
3. Translational files must not own attitude-coupled torque logic or actuator scheduling/commands.

## Scaffold Interface Set
1. `src/vehicle/resources/resources.jl`
2. `src/vehicle/resources/battery/battery_model.jl`
3. `src/vehicle/resources/solar_array/solar_array_model.jl`
4. `src/vehicle/resources/power_bus/power_bus_model.jl`
5. `src/vehicle/resources/loads/load_model.jl`
6. `src/mission/constellation/network/network.jl`
7. `src/gnc/estimation/estimation.jl`
8. `src/vehicle/actuators/laser_terminal/laser_terminal.jl`

## Enforced Gates
1. `test/ci_src_completeness_contract_gate.jl`
2. `test/ci_no_legacy_include_chains_gate.jl`
3. `test/ci_no_inert_force_terms_gate.jl`
4. `test/ci_scaffold_interface_gate.jl`
5. `test/ci_translational_ownership_gate.jl`
