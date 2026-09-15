---
id: environment.density_models__nrlmsise_ap_bins
label: _nrlmsise_ap_bins
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_ap_bins
  lines:
  - 530
  - 530
inputs:
- id: values
  type: Any
  units: n/a
  required: true
  description: Positional argument `values`.
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
  type: NTuple{8,
  units: n/a
  description: Return value of `_nrlmsise_ap_bins`.
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

# _nrlmsise_ap_bins

## Purpose
Validates a day's Ap record as exactly eight finite non-negative three-hour values.

## Design & Implementation
Requires length eight, converts to an `NTuple{8,Float64}`, and requires every entry finite and non-negative, raising `ArgumentError` otherwise.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `values` | Any | n/a | yes | Positional argument `values`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NTuple{8, | n/a | — | Return value of `_nrlmsise_ap_bins`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_space_indices_ap_bin|_nrlmsise_space_indices_ap_bin]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:542-542`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/atmosphere/density_models.jl:534-534`
<!-- vulcan:connections:end -->

## Limitations
Missing bins in the source data, represented as negative sentinels by some providers, are rejected as an error rather than interpolated.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 530.
