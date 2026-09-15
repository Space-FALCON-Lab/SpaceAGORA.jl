---
id: environment.density_models__nrlmsise_space_indices_f107a
label: _nrlmsise_space_indices_f107a
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_space_indices_f107a
  lines:
  - 522
  - 522
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
  type: Float64
  units: n/a
  description: Return value of `_nrlmsise_space_indices_f107a`.
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

# _nrlmsise_space_indices_f107a

## Purpose
Fetches the adjusted 81-day centred average F10.7 flux for the evaluation instant.

## Design & Implementation
Calls the lookup with `Val(:F10adj_avg_center81)` at `instant` and validates through `_nrlmsise_flux_value`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `lookup` | Any | n/a | yes | Positional argument `lookup`. |
| in | `instant` | DateTime | n/a | yes | Positional argument `instant`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_nrlmsise_space_indices_f107a`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_space_indices_indices|_nrlmsise_space_indices_indices]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:559-559`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.density_models__nrlmsise_flux_value|_nrlmsise_flux_value]] · `callers` · call · `src/environment/atmosphere/density_models.jl:523-523`
<!-- vulcan:connections:end -->

## Limitations
A centred average needs data forty days beyond the instant, so evaluations near the present day rely on the dataset's provisional values.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 522.
