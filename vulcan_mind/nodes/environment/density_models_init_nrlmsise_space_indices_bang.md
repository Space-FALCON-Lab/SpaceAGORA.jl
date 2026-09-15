---
id: environment.density_models_init_nrlmsise_space_indices_bang
label: init_nrlmsise_space_indices!
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: init_nrlmsise_space_indices!
  lines:
  - 279
  - 279
inputs:
- id: force_download
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `force_download` (default `false`).
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: Nothing
  units: n/a
  description: Return value of `init_nrlmsise_space_indices!`; mutates `force_download`
    in place. Returns `nothing`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# init_nrlmsise_space_indices!

## Purpose
Initialises the CelesTrak space-weather dataset once per process so NRLMSISE-00 evaluations do not pay download or parse cost mid-solve.

## Design & Implementation
Under `_NRLMSISE00_SPACE_INDICES_LOCK`, calls `SpaceIndices.init(SpaceIndices.Celestrak; force_download)` if forced or not yet ready, then sets the ready flag. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `force_download` | Bool | n/a | no | Keyword argument `force_download` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `init_nrlmsise_space_indices!`; mutates `force_download` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_provider_indices|_nrlmsise_provider_indices]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:610-610`
- [[environment.density_models_nrlmsise00spaceindicesprovider|NRLMSISE00SpaceIndicesProvider]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:269-269`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:232-232`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The ready flag is process-global, so a forced refresh in one run affects every other run in the process; there is no way to load a pinned historical dataset instead of the live one.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 279.
