---
id: environment.density_models__nrlmsise_space_indices_indices
label: _nrlmsise_space_indices_indices
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_space_indices_indices
  lines:
  - 557
  - 557
inputs:
- id: lookup
  type: Any
  units: n/a
  required: true
  description: Positional argument `lookup`.
- id: instant
  type: DateTime
  units: n/a
  required: true
  description: Positional argument `instant`.
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
  type: Any
  units: n/a
  description: Return value of `_nrlmsise_space_indices_indices`. Returns `(`.
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

# _nrlmsise_space_indices_indices

## Purpose
Bundles the three NRLMSISE-00 geophysical inputs read from the space-indices dataset into the named tuple the model's index resolver consumes.

## Design & Implementation
Calls `_nrlmsise_space_indices_f107a`, `_nrlmsise_space_indices_f107` and `_nrlmsise_space_indices_ap_vector` with the same lookup function and instant, and returns `(f107a=..., f107=..., ap=...)`. Taking the lookup as an argument rather than calling `SpaceIndices` directly is what lets tests inject a deterministic fake dataset.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `lookup` | Any | n/a | yes | Positional argument `lookup`. |
| in | `instant` | DateTime | n/a | yes | Positional argument `instant`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_nrlmsise_space_indices_indices`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_provider_indices|_nrlmsise_provider_indices]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:612-612`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.density_models__nrlmsise_space_indices_ap_vector|_nrlmsise_space_indices_ap_vector]] · `callers` · call · `src/environment/atmosphere/density_models.jl:561-561`
- `callees` → [[environment.density_models__nrlmsise_space_indices_f107|_nrlmsise_space_indices_f107]] · `callers` · call · `src/environment/atmosphere/density_models.jl:560-560`
- `callees` → [[environment.density_models__nrlmsise_space_indices_f107a|_nrlmsise_space_indices_f107a]] · `callers` · call · `src/environment/atmosphere/density_models.jl:559-559`
<!-- vulcan:connections:end -->

## Limitations
It performs no caching, so an evaluation at every RK stage repeats all eighteen underlying dataset queries even when successive instants fall in the same three-hour bin.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 557.
