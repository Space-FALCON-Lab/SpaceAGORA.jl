---
id: environment.density_models__nrlmsise_space_indices_ap_bin
label: _nrlmsise_space_indices_ap_bin
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_space_indices_ap_bin
  lines:
  - 541
  - 541
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
  description: Return value of `_nrlmsise_space_indices_ap_bin`.
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

# _nrlmsise_space_indices_ap_bin

## Purpose
Fetches the three-hour Ap value in effect at a given instant.

## Design & Implementation
Looks up the day's eight bins with `Val(:Ap)`, validates them, and indexes by `_nrlmsise_ap_slot_index`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `lookup` | Any | n/a | yes | Positional argument `lookup`. |
| in | `instant` | DateTime | n/a | yes | Positional argument `instant`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_nrlmsise_space_indices_ap_bin`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_space_indices_ap_vector|_nrlmsise_space_indices_ap_vector]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:548-548`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.density_models__nrlmsise_ap_bins|_nrlmsise_ap_bins]] · `callers` · call · `src/environment/atmosphere/density_models.jl:542-542`
- `callees` → [[environment.density_models__nrlmsise_ap_slot_index|_nrlmsise_ap_slot_index]] · `callers` · call · `src/environment/atmosphere/density_models.jl:543-543`
<!-- vulcan:connections:end -->

## Limitations
Each call re-fetches and re-validates the whole day's record, so the seven-slot builder that calls this sixteen times per evaluation performs sixteen dataset lookups.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 541.
