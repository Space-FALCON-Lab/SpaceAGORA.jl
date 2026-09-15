---
id: simx.engine_config_artifact_config_artifactconfig
label: ArtifactConfig
kind: struct
source:
  file: src/simulation/engine/config/artifact_config.jl
  symbol: ArtifactConfig
  lines:
  - 7
  - 10
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: engine_api
  type: Module
  units: n/a
  required: true
  description: RuntimeServices namespace and the SimulationEngine/SimulationCampaigns
    scope it aggregates, supplying the SPICE lock and the SimulationModel types used
    here.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: artifact_policy
  type: ArtifactConfig
  units: n/a
  description: Two-field policy record controlling result-bundle emission and deprecated-configuration
    warnings.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simx
origin: agent
---
# ArtifactConfig

## Purpose
`ArtifactConfig` is the typed record for what a run leaves on disk. It carries `save_bundle`, which decides whether the results bundle is written alongside the CSV, and `warn_deprecated_config`, which controls whether legacy configuration shapes emit a warning when they are still accepted.

## Model & Assumptions
Both fields are `Bool` and both default to true through `Base.@kwdef`, so a default-constructed value reproduces the historical behaviour of always writing the bundle and always warning. The struct is immutable, which means a run cannot flip artifact policy part-way through and produce a half-written bundle.

## Design & Implementation
The declaration is a four-line `Base.@kwdef struct` preceded by its docstring, and it is included by the engine module before `simulation_engine_config.jl` so the aggregate configuration can name it as a field type. Population happens in `simulation_engine_config_from_env`, which maps `SPACEAGORA_SAVE_BUNDLE` to `save_bundle` and `SPACEAGORA_WARN_DEPRECATED_CONFIG` to `warn_deprecated_config`. Keeping artifact policy in its own struct rather than as loose fields on the engine configuration means the persistence layer can be handed exactly the policy it needs without seeing solver or parallelism settings.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `artifact_policy` | ArtifactConfig | n/a | — | Two-field policy record controlling result-bundle emission and deprecated-configuration warnings. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.simulation_engine_config|SimulationEngineConfig]] · `callees` → `callers` · call · `src/simulation/engine/config/simulation_engine_config.jl:12-12`
- [[simx.engine_adapters_from_env_simulation_engine_config_from_env|simulation_engine_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:186-186`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The record says whether to write a bundle but not where; output directories and file naming are resolved separately through the `IOConfig` path helpers. There is no field for bundle compression, retention or schema selection, so the bundle schema version is a module constant rather than a configurable value.

## Provenance
Mapped from `src/simulation/engine/config/artifact_config.jl:7-10`; consumed by `simulation_engine_config_from_env` at `src/simulation/engine/adapters/from_env.jl:187`.
