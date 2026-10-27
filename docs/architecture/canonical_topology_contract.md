# Canonical topology contract

This contract says where each implementation responsibility lives in the
repository. The gates listed under [Enforcement](#enforcement) check the parts
that can be checked mechanically; the rest is applied in review.

## Source layout

`src/` has four physical groups, each holding a fixed set of domain
directories:

| Group | Domains |
|---|---|
| `src/interfaces/` | `cli`, `io`, `analysis` |
| `src/orchestration/` | `simulation`, `mission`, `parallel` |
| `src/models/` | `vehicle`, `environment`, `dynamics`, `gnc` |
| `src/foundation/` | `core`, `assets` |

The groups organize files only. They add no Julia modules, keep the existing
modules and their include order, and impose no dependency hierarchy:
`src/orchestration/simulation/composition/simulation_model.jl` still composes
definitions from every group. `src/SpaceAGORA.jl` and
`src/precompile_workload.jl` stay at the source root, and package extensions
stay in the repository-root `ext/`, where Julia discovers them. Optional
companions keep their own repositories and source trees outside core `src/`.

## Repository organization

The root README maps the repository's entry points. Organize files by
responsibility and intended use, not by author, date or language.

- Keep a capability's runner, configuration, analysis and instructions
  together. Study-specific tools belong with their study under `benchmarks/`;
  reusable maintenance and operational tools belong in `scripts/`. Share a
  helper once independent consumers need it; don't create a general utility
  area in anticipation.
- Reuse an existing home first. Add a directory only when it gives a distinct
  responsibility, dependency environment or lifecycle a clear boundary. File
  count alone is not a reason to split or merge directories.
- Runnable demonstrations belong in `examples/`, copy-and-complete extension
  skeletons in `templates/`, and automated pass/fail verification with its
  fixtures in `test/`. Long research campaigns stay with their benchmarks even
  when they call tested verification code.
- Use descriptive names. Keep Julia package conventions and supported command
  names; avoid names based only on a person, date or temporary revision.
- Remove redundant pointers, empty placeholders and superseded implementations
  once their useful content has a maintained home. Before retiring an
  executable, check its behavior, callers and published usage; having no
  callers does not by itself make it obsolete. Update references and checks in
  the same change, and keep a compatibility launcher only for a supported use.
- Keep generated run products in ignored output locations. Keep necessary
  fixtures, reproducibility inputs, licensed assets and documentation figures
  deliberately, with their purpose and provenance. Historical research records
  need an explicit retention or archive decision before removal.
- Use the root map and existing catalogs for discovery. Add a local README when
  it explains a real workflow or boundary, not to fill an otherwise empty
  folder. Extend an existing ownership check rather than adding a gate per
  placement rule.

Optional companions share the supported interface and compatibility rules, but
they load in different ways: `ext/` holds the glue Julia loads as package
extensions, and `packages/` records where each companion's source lives and how
it is installed. A change to a source pin must keep both the Git-checkout route
and the documented standalone installation working. Grouping source never moves
a companion's implementation into core.

## Ownership decisions

1. Benchmarks and studies live only in top-level `benchmarks/`, never under
   `src/`.
2. Runnable examples live only in top-level `examples/`, never under `src/`.
3. Reusable operational plotting scripts live in `scripts/plotting/`;
   study-specific plotting stays beside its study under `benchmarks/`.
   Installed analysis APIs stay under `src/interfaces/analysis/`.
4. Dynamics is owned by `src/models/dynamics/translational/`,
   `src/models/dynamics/rotational/`, `src/models/dynamics/effectors/` (every
   force and torque the simulator applies), `src/models/dynamics/multibody/`
   (articulated and compliant multibody dynamics) and
   `src/models/dynamics/external_propagation.jl` (the external-propagator
   hooks). `src/models/dynamics/models/` is retired.
5. Initial-state definitions live with the spacecraft model
   (`src/models/vehicle/spacecraft/model.jl`); mission plans and schedules live
   under `src/orchestration/mission/operations/`. One-run scenario configuration
   lives under `src/orchestration/simulation/config/`, where
   `constellation_configuration.jl` defines the spacecraft collection and the
   selected effectors while keeping the `SpacecraftModels.DynamicsModel` type
   identity. Physical spacecraft assembly stays under `src/models/vehicle/`.
6. Simulation execution belongs to `src/orchestration/simulation/engine/`. The
   former `src/orchestration/simulation/events/`,
   `src/orchestration/simulation/execution/` and
   `src/orchestration/simulation/solver_orchestration/` are retired.
7. Shared runtime serialization services (`SPICE_LOCK`, `GRAM_LOCK`) are owned
   by `src/orchestration/simulation/runtime_services.jl`, not by
   `SimulationModel`.
8. The actuator contract lives in `src/models/vehicle/actuators/actuators.jl`,
   mounting and evaluation of a mounted set in `mounting.jl`, and each actuator
   kind in its own file beside them (for example `thruster.jl`). Thruster firing
   dynamics live in `thruster_hooks.jl`.
9. Within the vehicle, `src/models/vehicle/spacecraft/` composes and integrates
   a spacecraft, and `src/models/vehicle/structure/` computes mass, inertia,
   geometry and structure, including the CAD and mesh readers
   (`mesh_geometry.jl`) that the visualization scene also uses.
10. Visualization: `src/interfaces/analysis/visualization/scene/` owns the
    renderer-independent scene layer (scene types, spacecraft geometry
    conversion, planet rotation sampling, the JSON sidecar) as a
    `SimulationModel` submodule, and
    `src/interfaces/analysis/visualization/rpo/` owns the PlotlyJS RPO plots.
    The browser viewer and HTML assembly belong to the independent
    `SpaceAGORAVisualizer.jl` package. Core keeps portable payload export, and
    `ext/SpaceAGORAVisualizerExt.jl` connects the optional renderer. No renderer
    source is copied under `src/viewer/` or a top-level `viewer/`.
11. Externally propagated spacecraft (shadow entries kept in the state vector
    while an optional package integrates them) are owned in four places, none
    of which knows a concrete integrator:
    - `src/models/dynamics/external_propagation.jl`: the
      `SimulationModel.ExternalPropagation` hook identities. Default hooks throw
      `Not implemented`, except `external_preflight`, whose default accepts.
    - `src/orchestration/simulation/engine/external_propagation.jl`: setup-time
      refusals, per-run runtimes and the shadow right-hand side.
    - `src/orchestration/simulation/callbacks/external_propagation_callbacks.jl`:
      the sync callback that overwrites shadow entries after every accepted
      step, installed before every other discrete callback.
    - `src/orchestration/simulation/external_propagation_services.jl`: the
      versioned `SpaceAGORA.ExternalPropagationServices` facade, which aliases
      the hook identities and supplies narrow configuration, state and frame
      services. Optional adapters use only this facade.

    `SimulationConfiguration.external_propagators` carries the owners. Concrete
    owners live in their own repositories, such as the
    [SpaceAGORAMuJoCo.jl companion](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl),
    and core never depends on them. Native bindings, adapter tests and the
    Basilisk/MuJoCo comparison harness belong to that companion; the generic
    lifecycle, conformance tests and shared physics stay in core.

## GNC and runtime IO owners

- GNC hooks: `src/models/gnc/control/hooks.jl`,
  `src/models/gnc/guidance/hooks.jl` and `src/models/gnc/navigation/hooks.jl`.
- Runtime IO: `src/interfaces/io/paths/io_paths.jl`,
  `src/interfaces/io/checkpoints/checkpoint_storage.jl` and
  `src/interfaces/io/outputs/io_outputs.jl`. The engine under
  `src/orchestration/simulation/engine/` orchestrates IO but owns no
  serialization or output implementation.
- No identifier in `src/` contains `legacy_`.

## Shared configuration and mathematics

- `SimulationModel.SimConfig._with_configuration` copies a simulation
  configuration with named overrides. Unchanged fields keep their references,
  and the constructor infers model types again. Callers stay responsible for
  nested settings and deliberate resets; for example, the telemetry runner
  clears the typed solver (`solver_config=nothing`) because it owns solver
  selection and its per-attempt iteration retry. Mutable-state isolation is the engine's job, not this
  copy's. Settings, mission, environment, spacecraft and campaign-sample copies
  use `SimConfig._with_fields` the same way: name only the fields that change,
  and every other field passes through the type's constructor.
- `SimulationModel.QuaternionMath` owns the quaternion functions in
  `src/foundation/core/numerics/quaternion_utils.jl`. Consumers import this
  module rather than including another copy. Quaternions are scalar-last, and
  `rot(q)` maps inertial components to body components.
- `SimulationModel.FrameTransforms` owns the coordinate functions in
  `src/foundation/core/interfaces/reference_system.jl`. The engine, callbacks,
  GNC, aerodynamics and telemetry import these same bindings.
  `rotate_vector_by_quaternion` is an active vector rotation, so it is not
  interchangeable with `rot(q)`.
- These are internal owners, not root-level public APIs. Existing qualified
  bindings stay available, and `ReferenceSystems` still owns the
  reference-system types.

## Enforcement

The retired paths are listed in `RETIRED_PATHS` in
`test/gates/ci_path_policy_gate.jl`, not on this page; record a new retirement
there instead of adding a gate.

- `test/gates/ci_path_policy_gate.jl` checks the source layout above (exactly
  these groups and domains, real directories, each with Julia files; no former
  `src/<domain>/` root; the two root files; an entry file in `ext/` for every
  extension registered in `Project.toml`). It fails if a retired path exists or
  a required owner file in `REQUIRED_PATHS` is missing, and it scans source,
  tests, docs, examples and workflows for forbidden tokens: retired paths, the
  `legacy_` prefix and the retired `ControlEffectors`, `GuidanceEffectors` and
  `NavigationEffectors` names in `src/`. This page,
  `src_completeness_contract.md` and `src_canonical_owner_audit.md` are exempt
  from the token scan because they name retired paths on purpose.
- `test/gates/ci_architecture_contract_gate.jl` checks module composition and
  include order, lock ownership in `RuntimeServices`, the typed
  `run_simulation` boundary and the effector split. It also fails if any of
  these owner files is missing:
  `src/orchestration/simulation/engine/public_api.jl`,
  `src/orchestration/simulation/engine/execution.jl`,
  `src/orchestration/simulation/runtime_services.jl`,
  `src/models/dynamics/effectors/force_torque_models.jl`,
  `src/models/vehicle/actuators/thruster_hooks.jl`,
  `src/models/vehicle/structure/structure_models.jl`,
  `src/models/vehicle/structure/mass_properties.jl` and
  `src/models/vehicle/structure/geometry_properties.jl`.
- `test/gates/ci_vehicle_structure_boundary_gate.jl` requires the structure
  owner files and fails if a structure function (`get_COM`,
  `get_inertia_tensor` and the area, normal-vector and tangent-vector helpers)
  is defined under `src/models/vehicle/spacecraft/`.
- `test/gates/ci_thin_entry_files_gate.jl` keeps the listed module entry files
  and benchmark launchers thin: at most 80 non-comment lines and no function or
  struct definitions, and the two benchmark launchers forward to their split
  implementations.
