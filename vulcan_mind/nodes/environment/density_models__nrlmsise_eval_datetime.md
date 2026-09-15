---
id: environment.density_models__nrlmsise_eval_datetime
label: _nrlmsise_eval_datetime
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_eval_datetime
  lines:
  - 502
  - 502
inputs:
- id: initial_time
  type: Any
  units: n/a
  required: true
  description: Positional argument `initial_time`.
- id: el_time
  type: Float64
  units: n/a
  required: true
  description: Positional argument `el_time`.
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
  type: DateTime
  units: n/a
  description: Return value of `_nrlmsise_eval_datetime`.
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

# _nrlmsise_eval_datetime

## Purpose
Converts the run's initial time plus elapsed seconds into the `DateTime` NRLMSISE-00 and the space-indices lookup require.

## Design & Implementation
Builds a `DateTime` from the integer year, month, day, hour and minute with zero seconds, adds the fractional seconds as milliseconds, then adds elapsed time as milliseconds. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `initial_time` | Any | n/a | yes | Positional argument `initial_time`. |
| in | `el_time` | Float64 | n/a | yes | Positional argument `el_time`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | DateTime | n/a | — | Return value of `_nrlmsise_eval_datetime`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__gram_point_density|_gram_point_density]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:1077-1077`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/atmosphere/density_models.jl:510-510`
<!-- vulcan:connections:end -->

## Limitations
Millisecond rounding of both the epoch seconds and elapsed time limits temporal resolution to one millisecond, which is irrelevant for atmosphere but worth knowing for any reuse.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 502.
