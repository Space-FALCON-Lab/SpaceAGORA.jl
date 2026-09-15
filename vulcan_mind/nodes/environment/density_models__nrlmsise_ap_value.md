---
id: environment.density_models__nrlmsise_ap_value
label: _nrlmsise_ap_value
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_ap_value
  lines:
  - 319
  - 319
inputs:
- id: ap
  type: Real
  units: n/a
  required: true
  description: Positional argument `ap`.
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
  description: Return value of `_nrlmsise_ap_value`.
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

# _nrlmsise_ap_value

## Purpose
Validates and normalises the Ap input, accepting either a single daily value or the full seven-slot history NRLMSISE-00 can use.

## Design & Implementation
Two methods. A `Real` is converted to `Float64` and required finite and non-negative. A vector must have length seven and every element finite and non-negative, and is returned as an `SVector{7,Float64}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ap` | Real | n/a | yes | Positional argument `ap`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_nrlmsise_ap_value`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_flux_value|_nrlmsise_flux_value]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:371-371`
- [[environment.density_models__nrlmsise_provider_indices|_nrlmsise_provider_indices]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:572-572`
- [[environment.density_models__nrlmsise_space_indices_ap_vector|_nrlmsise_space_indices_ap_vector]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:547-547`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/atmosphere/density_models.jl:320-320`
<!-- vulcan:connections:end -->

## Limitations
The seven-slot ordering — daily, current three-hour, three, six and nine hours prior, and two eight-bin averages — is assumed from the caller and not documented in the error messages.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 319.
