---
id: environment.density_models_nrlmsise00spaceindicesprovider
label: NRLMSISE00SpaceIndicesProvider
kind: struct
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: NRLMSISE00SpaceIndicesProvider
  lines:
  - 264
  - 264
inputs:
- id: force_download
  type: Bool
  units: n/a
  required: true
  description: Field `force_download`.
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
  type: NRLMSISE00SpaceIndicesProvider
  units: n/a
  description: Constructed `NRLMSISE00SpaceIndicesProvider`.
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

# NRLMSISE00SpaceIndicesProvider

## Purpose
The built-in index provider that feeds NRLMSISE-00 with real F10.7 and Ap history from CelesTrak through the `SpaceIndices` package.

## Design & Implementation
A mutable struct with a single `force_download` flag. It is callable with `(instant, h, lat, lon)`: below 80 km it returns the standard defaults of 150, 150 and 4 without touching the dataset; otherwise it initialises the dataset once through `init_nrlmsise_space_indices!`, clears its own force flag, and returns the adjusted centred 81-day average, the previous day's flux and the seven-slot Ap vector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `force_download` | Bool | n/a | yes | Field `force_download`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NRLMSISE00SpaceIndicesProvider | n/a | — | Constructed `NRLMSISE00SpaceIndicesProvider`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_flux_value|_nrlmsise_flux_value]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:372-372`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.density_models_init_nrlmsise_space_indices_bang|init_nrlmsise_space_indices!]] · `callers` · call · `src/environment/atmosphere/density_models.jl:269-269`
- `callees` → [[environment.density_models_nrlmsise00atmospheremodel|NRLMSISE00AtmosphereModel]] · `callers` · call · `src/environment/atmosphere/density_models.jl:272-272`
<!-- vulcan:connections:end -->

## Limitations
The first evaluation above 80 km can trigger a network download inside the integrator unless the dataset was prewarmed; the 80 km cutoff is a literal constant.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 264.
