# Source completeness contract

This contract classifies the files under `src/` and lists what canonical code
must not contain. Each rule names the gate in `test/gates/` that enforces it; a
rule without a gate is applied in review.

## File categories

A canonical owner owns runtime behavior: it contains executable logic, typed
data structures or concrete API methods.

A canonical aggregator only includes and re-exports. It starts with the comment
`# Canonical aggregator: no behavior ownership.` and defines no function or
type. `ci_src_completeness_contract_gate.jl` checks the header and parses each
file for definitions; its `CANONICAL_AGGREGATOR_FILES` list is:

- `src/models/gnc/planning/rpo/lifecycle.jl`
- `src/orchestration/simulation/engine/simulation_engine.jl`
- `src/orchestration/simulation/callbacks/callbacks.jl`
- `src/orchestration/parallel/routing/parallel_profiles.jl`
- `src/orchestration/parallel/process/parallel_process.jl`
- `src/orchestration/simulation/composition/simulation_model.jl`

An experimental interface is a typed contract for a planned subsystem. It lives
under top-level `experimental/`, never under `src/`, and stays out of the
package load graph until it owns real behavior. Each unimplemented method throws
`ErrorException("Not implemented: <method>")`; `ci_scaffold_interface_gate.jl`
loads the scaffolds and checks every entry point for that message prefix.

`ci_src_completeness_contract_gate.jl` refuses these in any `src/` file:

- silent stubs of the form `*_stub() = nothing`;
- the header `Compatibility wrapper: canonical path forwarding to legacy implementation.`,
  which marked a file that only forwarded to a legacy owner;
- the token `No-op `.

A force term computed in an active dynamics path but never applied is also
forbidden; no gate detects it.

## Translational boundary

1. `src/models/dynamics/translational/` holds pure translational terms.
   `ci_translational_ownership_gate.jl` refuses inline position, velocity and
   mass derivative assignments (`du_view.pos .= ...`, `du_view.vel .= ...`,
   `du_view.mass = ...`) anywhere else in `src/`. `ci_architecture_contract_gate.jl`
   requires `translational_models.jl`, `position_kinematics.jl` and
   `point_mass_dynamics.jl` there and checks that the engine right-hand side
   calls the `DynamicsTranslational.assign_*_translational_rhs!` owners.
2. Gravity field kernels live under `src/models/environment/gravity/`, and the
   spacecraft and runtime adapters under `src/models/dynamics/effectors/gravity/`.
   The ordered include manifest `src/models/environment/gravity/gravity_models.jl`
   keeps the existing `DynamicEffectors.GravityEffectors` module ownership and
   definition order. Neither belongs under translational dynamics: the
   architecture gate requires the kernel and adapter files and fails if
   `src/models/dynamics/translational/gravity_models.jl` exists.
3. Translational files do not own attitude-coupled torque logic or actuator
   scheduling and commands (review rule).

## Topology boundary

`ci_path_policy_gate.jl` refuses each retired path below (its `RETIRED_PATHS`
list holds the full set).

1. Benchmarks live in top-level `benchmarks/`; `src/benchmarks/` is retired.
2. Runnable examples live in top-level `examples/`; `src/examples/` is retired.
3. Plotting and report tooling lives in top-level `scripts/`;
   `src/interfaces/analysis/plotting/` and `src/interfaces/analysis/reports/`
   are retired.
4. Unfinished subsystem scaffolds live in top-level `experimental/`, and copies
   under `src/` are retired (`src/models/vehicle/resources/`,
   `src/models/gnc/estimation/`, `src/orchestration/mission/constellation/`).
5. Plotting, reporting and parallel support code has no owner outside the roots
   this contract names.
6. Force and torque effectors live under `src/models/dynamics/effectors/`;
   multibody dynamics lives under `src/models/dynamics/multibody/`.
7. Simulation execution is owned by
   `src/orchestration/simulation/engine/public_api.jl` and
   `src/orchestration/simulation/engine/execution.jl`;
   `src/orchestration/simulation/events/`,
   `src/orchestration/simulation/execution/` and
   `src/orchestration/simulation/solver_orchestration/` are retired.
8. Initial-state definitions live with the spacecraft model
   (`src/models/vehicle/spacecraft/model.jl`); mission plans and policies live
   under `src/orchestration/mission/operations/`.
9. Structural analysis lives in `src/models/vehicle/structure/`, not
   `src/models/vehicle/spacecraft/`, and the CAD and mesh readers in
   `src/models/vehicle/structure/mesh_geometry.jl`.
   `ci_vehicle_structure_boundary_gate.jl` fails if a structure function such as
   `get_COM` or `get_inertia_tensor` is defined under `spacecraft/`.
10. The visualization scene is owned by
    `src/interfaces/analysis/visualization/scene/`, a canonical owner inside
    `SimulationModel`, and portable payload export stays with it. The browser
    viewer lives in the independent `SpaceAGORAVisualizer.jl`, which
    `ext/SpaceAGORAVisualizerExt.jl` connects; `src/viewer/` is retired.

## Include chains

`ci_no_legacy_include_chains_gate.jl` scans the core, environment, GNC,
simulation and parallel source trees. It allows a raw `include(...)` only in the
files on its `ALLOWED_RAW_INCLUDE_FILES` list (ordered model manifests, module
aggregators and composition files), refuses the `__legacy_` token, and refuses
includes that climb three directories
(`include(joinpath(@__DIR__, "..", "..", "..", ...`) outside the two
composition seams on its `CROSS_GROUP_COMPOSITION_OWNERS` list.

## GNC aerobraking boundary

1. Typed energy-depletion guidance (EDG) is owned by
   `src/models/gnc/guidance/aerobraking/` (`typed_edg/`,
   `energy_depletion_lifecycle.jl` and `edg_guidance_family.jl`).
2. The older E-EDG and T-EDG strategy routes are gone: their guidance folders
   (`aerobraking/e_edg/`, `t_edg/`, `common/`), their control files
   (`control/aerobraking/constraint_tracking.jl`, `control_commands.jl` and
   `tracking_executor.jl`; that folder now holds the typed EDG control files),
   `src/models/gnc/internal/`, the EDG compatibility bindings and the
   mission-owned `src/orchestration/mission/operations/aerobraking_policy/`
   selector. `ci_path_policy_gate.jl` keeps all of these retired.
3. Guidance files do not include control source files.
   `ci_gnc_aerobraking_boundary_gate.jl` checks this under
   `guidance/aerobraking/`, and `ci_no_guidance_control_cross_include_gate.jl`
   checks it across the guidance, planning, companion and contract trees.
4. EDG control code in `src/models/gnc/control/` defines no numerical EDG
   solver; it calls the companion's control decision through `EDGAlgorithms`.
5. Propulsive maneuver guidance and control are coupled only through typed
   command buffers; guidance does not change thruster runtime state.
   `ci_gnc_typed_command_boundary_gate.jl` checks that guidance writes
   `PropulsiveManeuverCommand` values into `shared_buffers.maneuver_commands`
   and control reads them from there.
6. The typed EDG numerical kernels live in the optional `SpaceAGORAEDG.jl`
   package. Core `typed_edg/algorithms.jl` is not an aggregator: it holds
   exactly the two package forwarders, `guidance_decision!` and
   `control_decision!`, with no include, no name list and no numerical code.
   The facade's `IMPLEMENTED_FUNCTIONS` in
   `src/models/gnc/companions/edg/services.jl` lists exactly those two names,
   and `typed_edg/services.jl` includes only `adapter_queries.jl`, the
   capability wrapper through which the companion receives copied values rather
   than raw integrator state or mutable vehicle geometry. These checks live in
   `test/contracts/edg_ownership_checks.jl`. `ci_gnc_aerobraking_boundary_gate.jl`
   runs all of them, including a check that no retired kernel file
   (`heat_rate.jl`, `targeting.jl` and the like) reappears in `typed_edg/`;
   `ci_src_completeness_contract_gate.jl` and
   `ci_no_legacy_include_chains_gate.jl` also check the two forwarders in
   `algorithms.jl`.

## Experimental interface set

1. `experimental/vehicle/resources/resources.jl`
2. `experimental/vehicle/resources/battery/battery_model.jl`
3. `experimental/vehicle/resources/solar_array/solar_array_model.jl`
4. `experimental/vehicle/resources/power_bus/power_bus_model.jl`
5. `experimental/vehicle/resources/loads/load_model.jl`
6. `experimental/mission/constellation/network/network.jl`
7. `experimental/gnc/estimation/estimation.jl`

## Gates

1. `test/gates/ci_src_completeness_contract_gate.jl`
2. `test/gates/ci_no_legacy_include_chains_gate.jl`
3. `test/gates/ci_scaffold_interface_gate.jl`
4. `test/gates/ci_translational_ownership_gate.jl`
5. `test/gates/ci_path_policy_gate.jl`
6. `test/gates/ci_vehicle_structure_boundary_gate.jl`
7. `test/gates/ci_gnc_aerobraking_boundary_gate.jl`
8. `test/gates/ci_no_guidance_control_cross_include_gate.jl`
9. `test/gates/ci_gnc_typed_command_boundary_gate.jl`
10. `test/gates/ci_architecture_contract_gate.jl`
