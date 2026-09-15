---
id: environment.density_models__nrlmsise_space_indices_f107
label: _nrlmsise_space_indices_f107
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_space_indices_f107
  lines:
  - 518
  - 518
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
  description: Return value of `_nrlmsise_space_indices_f107`.
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

# _nrlmsise_space_indices_f107

## Purpose
Fetches the previous day's adjusted F10.7 solar flux for an NRLMSISE-00 evaluation, following the model's convention that the daily flux input lags the evaluation date by one day.

## Design & Implementation
Calls the injected lookup with `Val(:F10adj)` at `instant - Day(1)` and passes the result through `_nrlmsise_flux_value`, which enforces finiteness and non-negativity and names the field in any error. Declared `@inline` with a `::Float64` return so the provider's index assembly has no call overhead.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `lookup` | Any | n/a | yes | Positional argument `lookup`. |
| in | `instant` | DateTime | n/a | yes | Positional argument `instant`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_nrlmsise_space_indices_f107`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_space_indices_indices|_nrlmsise_space_indices_indices]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:560-560`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.density_models__nrlmsise_flux_value|_nrlmsise_flux_value]] · `callers` · call · `src/environment/atmosphere/density_models.jl:519-519`
<!-- vulcan:connections:end -->

## Limitations
A date outside the CelesTrak dataset's coverage raises from inside `SpaceIndices` rather than falling back to the 150 sfu default, so a simulation epoch in the far future fails at the first above-80-km evaluation.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 518.
