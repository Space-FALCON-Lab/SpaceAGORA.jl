---
id: environment.density_models__nrlmsise_ap_slot_index
label: _nrlmsise_ap_slot_index
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_ap_slot_index
  lines:
  - 526
  - 526
inputs:
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
  type: Int
  units: n/a
  description: Return value of `_nrlmsise_ap_slot_index`.
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

# _nrlmsise_ap_slot_index

## Purpose
Maps a time of day onto the one-based index of its three-hour Ap bin.

## Design & Implementation
Returns `fld(hour, 3) + 1`, giving one through eight. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `instant` | DateTime | n/a | yes | Positional argument `instant`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_nrlmsise_ap_slot_index`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_space_indices_ap_bin|_nrlmsise_space_indices_ap_bin]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:543-543`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Uses the `DateTime`'s own hour without timezone adjustment; the bins are defined in UT, which is correct only if the instant is UT.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 526.
