---
id: environment.density_models__nrlmsise_space_indices_ap_vector
label: _nrlmsise_space_indices_ap_vector
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_space_indices_ap_vector
  lines:
  - 546
  - 546
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
  type: SVector{7,
  units: n/a
  description: Return value of `_nrlmsise_space_indices_ap_vector`.
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

# _nrlmsise_space_indices_ap_vector

## Purpose
Assembles the seven-element Ap history NRLMSISE-00 uses for geomagnetic heating: daily, current, three prior bins and two eight-bin averages.

## Design & Implementation
Reads the daily Ap, the bins at the instant and three, six and nine hours earlier, and averages the bins at twelve through thirty-three and thirty-six through fifty-seven hours earlier in three-hour steps. Returns an `SVector{7,Float64}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `lookup` | Any | n/a | yes | Positional argument `lookup`. |
| in | `instant` | DateTime | n/a | yes | Positional argument `instant`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{7, | n/a | — | Return value of `_nrlmsise_space_indices_ap_vector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_space_indices_indices|_nrlmsise_space_indices_indices]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:561-561`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.density_models__nrlmsise_ap_value|_nrlmsise_ap_value]] · `callers` · call · `src/environment/atmosphere/density_models.jl:547-547`
- `callees` → [[environment.density_models__nrlmsise_space_indices_ap_bin|_nrlmsise_space_indices_ap_bin]] · `callers` · call · `src/environment/atmosphere/density_models.jl:548-548`
<!-- vulcan:connections:end -->

## Limitations
Sixteen bin lookups per evaluation with no caching; a density model evaluated at every RK stage pays this repeatedly for instants that share a day.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 546.
