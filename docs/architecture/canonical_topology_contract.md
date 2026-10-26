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

## Repository organization

The root README maps the repository's entry points. Organize files by their
responsibility and intended use, rather than their author, date or language.

- Keep a capability's runner, configuration, analysis and instructions together.
  Study-specific tools belong with their study under `benchmarks/`; reusable
  maintenance and operational tools belong in `scripts/`. Share helpers when
  independent consumers need them, rather than creating a general utility area
  in anticipation of future use.
- Reuse an existing home first. Add a directory only when it gives a distinct
  responsibility, dependency environment or lifecycle a clear boundary. A file
  count alone is not a reason to split or combine directories.
- Runnable demonstrations belong in `examples/`. Copy-and-complete extension
  skeletons belong in `templates/`. Automated pass/fail verification and its
  required fixtures belong in `test/`; long research campaigns belong with
  their benchmarks, even when they call tested verification code.
- Use descriptive names. Preserve Julia package conventions and supported
  command names; avoid names based only on a person, date or temporary revision.
- Remove redundant pointers, empty placeholders and superseded implementations
  when their useful content has a maintained home. Before retiring an executable,
  inspect its behavior, callers and published usage: lack of callers alone does
  not establish that it is obsolete. Update references and checks in the same
  change. Keep compatibility launchers only for a supported usage requirement.
- Keep generated run products in ignored output locations. Retain necessary
  fixtures, reproducibility inputs, licensed assets and documentation figures
  deliberately, with their purpose and provenance. Historical research records
  need an explicit retention or archive decision before removal.
- Use the root map and existing catalogs for discovery. Add a local README when
  it explains a real workflow or boundary; do not create one just to fill an
  otherwise empty folder. Extend existing ownership checks rather than creating
  a separate gate for each placement rule.

Optional companions share supported interface and compatibility requirements,
but their loading mechanisms can differ. `ext/` contains Julia-discovered glue;
`packages/` locates external source and installation metadata. Changes to source
pins must preserve both Git-checkout and documented standalone installation
routes. Source grouping does not move companion implementations into core.

## Locked Ownership Decisions
1. Benchmarks and studies are owned only by top-level `benchmarks/` (never under `src/`).
2. Runnable examples are owned only by top-level `examples/` (never under `src/`).
3. Reusable operational plotting scripts are owned by top-level `scripts/plotting/`; study-specific plotting belongs beside its study under `benchmarks/`. Installed analysis APIs remain under `src/interfaces/analysis/`.
4. Domain dynamics ownership is:
   - `src/models/dynamics/translational/*`
   - `src/models/dynamics/rotational/*`
   - `src/models/dynamics/effectors/*` (every force and torque the simulator applies)
   - `src/models/dynamics/multibody/*` (articulated and compliant multibody dynamics)
   - `src/models/dynamics/external_propagation.jl` (the external-propagator hook)
   and `src/models/dynamics/models/` is retired.
5. Mission configuration ownership: initial-state definitions live with the
   spacecraft model (`src/models/vehicle/spacecraft/model.jl`); mission plans and
   schedules live under `src/orchestration/mission/operations/*`.
   One-run scenario configuration lives under `src/orchestration/simulation/config/`; the existing
   spacecraft collection and selected effectors are defined in
   `constellation_configuration.jl`, retaining the `SpacecraftModels.DynamicsModel`
   identity. Physical spacecraft assembly remains under vehicle.
6. Simulation execution ownership belongs to `src/orchestration/simulation/engine/*`; legacy `src/orchestration/simulation/events/*`, `src/orchestration/simulation/execution/*`, and `src/orchestration/simulation/solver_orchestration/*` are retired.
7. Shared runtime serialization services are owned by `src/orchestration/simulation/runtime_services.jl`, not by `SimulationModel`.
8. Actuator ownership: the actuator contract lives in
   `src/models/vehicle/actuators/actuators.jl`, mounting and the evaluation of a
   mounted set in `src/models/vehicle/actuators/mounting.jl`, and each actuator
   kind in one file beside them (`thruster.jl`); thruster firing dynamics live in
   `src/models/vehicle/actuators/thruster_hooks.jl`.
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
    - `src/models/dynamics/external_propagation.jl`: the `SimulationModel.ExternalPropagation`
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
4. `src/models/dynamics/effectors/force_torque_models.jl`
5. `src/models/vehicle/actuators/thruster_hooks.jl`
6. `src/models/vehicle/structure/structure_models.jl`
7. `src/models/vehicle/structure/mass_properties.jl`
7a. `src/models/vehicle/structure/mesh_geometry.jl`
8. `src/models/vehicle/structure/geometry_properties.jl`

## Path Policy
1. This contract documents canonical owners only.
2. Removed or forbidden legacy paths are intentionally not listed here as a
   tombstone inventory.
3. Canonical-path enforcement belongs to CI gates, not to user-facing topology
   documentation.

## Canonical GNC and Runtime IO Owners
1. GNC hook owners: `src/models/gnc/control/hooks.jl`, `src/models/gnc/guidance/hooks.jl`, `src/models/gnc/navigation/hooks.jl`.
2. Runtime IO owners: `src/interfaces/io/paths/io_paths.jl`, `src/interfaces/io/checkpoints/checkpoint_storage.jl`, `src/interfaces/io/outputs/io_outputs.jl`; `src/orchestration/simulation/engine/*` orchestrates IO but does not own serialization or output implementations.
3. No `legacy_` identifiers in `src/`.

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
  Settings, mission, environment, spacecraft and campaign-sample copies use
  `SimConfig._with_fields` the same way: they name only the fields that change,
  and every other field passes through the type's constructor.
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
