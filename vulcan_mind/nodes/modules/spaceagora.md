---
id: module.spaceagora
label: SpaceAGORA
kind: module
source:
  file: src/SpaceAGORA.jl
  symbol: SpaceAGORA
inputs:
- id: analysis
  type: Module
  units: n/a
  required: true
  description: TelemetryVerification and RPOVisualization, included at lines 13 and
    15, supplying VerificationRequest, run_verification, run_study, rpo_path_plot
    and rpo_tracking_plot.
- id: assets
  type: Module
  units: n/a
  required: true
  description: RPOStationAssets, included at line 14, supplying station_geometry_path,
    station_cad_path and the load_rpo_station_pointcloud family of loaders.
- id: cli
  type: Module
  units: n/a
  required: true
  description: SpaceAGORACLI, included at line 16, supplying AssetCheckItem, AssetCheckReport,
    check_assets, render_asset_report and run_cli.
- id: core
  type: Module
  units: n/a
  required: true
  description: SimulationModel, included at line 10, the aggregator that owns spacecraft,
    environment, dynamics, guidance and control types re-exported by this package.
- id: ext
  type: Module
  units: n/a
  required: true
  description: SpaceAGORAGRAMSuiteExt, a Julia package extension loaded by the runtime
    when GRAMSuite is present, weakly extending the density-model surface.
- id: parallel
  type: Module
  units: n/a
  required: true
  description: ParallelProfiles and ParallelProcess, included at lines 7 and 8, supplying
    ParallelProfile, with_parallel_profile, ProcessPool and campaign_process_pool.
- id: simulation
  type: Module
  units: n/a
  required: true
  description: RuntimeServices, SimulationEngine and SimulationCampaigns, included
    at lines 9, 11 and 12, supplying SimulationEngineConfig, run_simulation and run_monte_carlo.
outputs:
- id: api
  type: Module
  units: n/a
  description: The flat public namespace declared by the export list at lines 435
    to 489, plus the forwarding entrypoints run_simulation, run_cli, check_assets,
    render_asset_report, prewarm_nbody_ephemeris_cache and load_nbody_ephemeris_cache!.
tags:
- module
charts:
- master
origin: agent
---

# SpaceAGORA

## Purpose
`SpaceAGORA` is the top-level Julia package module. It owns no physics: it includes the ten
implementation sub-modules in dependency order, re-binds their exported symbols into a single
flat namespace, attaches documentation to those bindings, and defines six thin forwarding
entrypoints that constitute the package's supported call surface.

## Model & Assumptions
- Include order at lines 7 to 16 is load-bearing. `parallel/routing/parallel_profiles.jl` must
  precede `simulation/runtime_services.jl`, which must precede `core/simulation_model.jl`,
  because each later file resolves names through `parentmodule(@__MODULE__)`.
- Every symbol re-exported here is assumed to be owned by exactly one sub-module. The `using
  .Submodule: name` form is explicit rather than blanket, so a name collision between two
  sub-modules surfaces as a load-time error instead of a silent shadow.
- `@doc (@doc Submodule.X) X` at lines 76 to 184 assumes the docs build resolves `@docs
  SpaceAGORA.X` against this module's own doc metadata rather than following the import alias.
- `run_simulation` assumes `isolate_state=true` by default, deep-copying the configuration so
  repeated or concurrent runs cannot mutate shared campaign state.

## Design & Implementation
The module body is three phases. Phase one, lines 7 to 16, evaluates ten `include` calls whose
paths are built with `joinpath(@__DIR__, ...)`. Phase two, lines 18 to 86, imports named symbols
from `.ParallelProfiles`, `.ParallelProcess`, `.SimulationEngine`, `.SimulationCampaigns`,
`.SimulationModel`, `.TelemetryVerification`, `.RPOStationAssets`, `.RPOVisualization` and
`.SpaceAGORACLI`; `prewarm_nbody_ephemeris_cache` and `load_nbody_ephemeris_cache!` use `import`
rather than `using` because lines 533 and 544 add package-level methods to them. Phase three
declares the exports and the forwarding methods. `run_simulation` has two methods, one generic
and one specialised on `SimulationEngineConfig`, both delegating to
`SimulationEngine.run_simulation`. The final `include` at line 575 pulls in
`precompile_workload.jl`, which runs a five-second Mars aerobraking trajectory under
`@compile_workload` so the package ships with a warm method cache.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `analysis` | Module | n/a | yes | TelemetryVerification and RPOVisualization, included at lines 13 and 15, supplying VerificationRequest, run_verification, run_study, rpo_path_plot and rpo_tracking_plot. |
| in | `assets` | Module | n/a | yes | RPOStationAssets, included at line 14, supplying station_geometry_path, station_cad_path and the load_rpo_station_pointcloud family of loaders. |
| in | `cli` | Module | n/a | yes | SpaceAGORACLI, included at line 16, supplying AssetCheckItem, AssetCheckReport, check_assets, render_asset_report and run_cli. |
| in | `core` | Module | n/a | yes | SimulationModel, included at line 10, the aggregator that owns spacecraft, environment, dynamics, guidance and control types re-exported by this package. |
| in | `ext` | Module | n/a | yes | SpaceAGORAGRAMSuiteExt, a Julia package extension loaded by the runtime when GRAMSuite is present, weakly extending the density-model surface. |
| in | `parallel` | Module | n/a | yes | ParallelProfiles and ParallelProcess, included at lines 7 and 8, supplying ParallelProfile, with_parallel_profile, ProcessPool and campaign_process_pool. |
| in | `simulation` | Module | n/a | yes | RuntimeServices, SimulationEngine and SimulationCampaigns, included at lines 9, 11 and 12, supplying SimulationEngineConfig, run_simulation and run_monte_carlo. |
| out | `api` | Module | n/a | — | The flat public namespace declared by the export list at lines 435 to 489, plus the forwarding entrypoints run_simulation, run_cli, check_assets, render_asset_report, prewarm_nbody_ephemeris_cache and load_nbody_ephemeris_cache!. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.analysis|TelemetryVerification]] · `api` → `analysis` · call · `src/SpaceAGORA.jl:12-14`
- [[module.assets|RPOStationAssets]] · `api` → `assets` · call · `src/SpaceAGORA.jl:13-13`
- [[module.cli|SpaceAGORACLI]] · `api` → `cli` · call · `src/SpaceAGORA.jl:15-15`
- [[module.core|SimulationModel]] · `api` → `core` · call · `src/SpaceAGORA.jl:9-9`
- [[module.ext|SpaceAGORAGRAMSuiteExt]] · `api` → `ext` · call · `Project.toml:54-54`
- [[module.parallel|ParallelProfiles]] · `api` → `parallel` · call · `src/SpaceAGORA.jl:6-7`
- [[module.simulation|RuntimeServices]] · `api` → `simulation` · call · `src/SpaceAGORA.jl:8-11`

**Downstream**

- `api` → [[spaceagora.spaceagora_load_nbody_ephemeris_cache_bang|load_nbody_ephemeris_cache!]] · `module_api` · call · `src/SpaceAGORA.jl`
- `api` → [[spaceagora.spaceagora_prewarm_nbody_ephemeris_cache|prewarm_nbody_ephemeris_cache]] · `module_api` · call · `src/SpaceAGORA.jl`
- `api` → [[spaceagora.spaceagora_run_cli|run_cli]] · `module_api` · call · `src/SpaceAGORA.jl`
- `api` → [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `module_api` · call · `src/SpaceAGORA.jl`
<!-- vulcan:connections:end -->

## Limitations
Because the sub-modules are `include`d rather than loaded as separate packages, a syntax or
type error anywhere under `src/` fails the whole package load; there is no partial-import path.
The flat namespace means renaming a sub-module symbol silently drops it from the package API
unless the corresponding `using` and `export` lines are updated together. The precompile
workload writes into `tempdir()` and shells `cd` into a `mktempdir`, so package precompilation
requires a writable temporary directory. The GRAM extension is weak: on a machine without
`GRAMSuite` the GRAM-backed density models are absent, and only the analytic
`NoAtmosphereModel`, `ExponentialAtmosphereModel`, `PiecewiseExponentialAtmosphereModel` and
`NRLMSISE00AtmosphereModel` paths remain callable.

## Provenance
Mapped from `src/SpaceAGORA.jl` and `src/precompile_workload.jl`.
