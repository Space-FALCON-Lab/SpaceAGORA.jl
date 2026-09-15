---
id: environment.density_models__batch_elapsed_time
label: _batch_elapsed_time
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _batch_elapsed_time
  lines:
  - 835
  - 835
inputs:
- id: el_time
  type: Float64
  units: n/a
  required: true
  description: Positional argument `el_time`.
- id: idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `idx`.
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
  description: Return value of `_batch_elapsed_time`.
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

# _batch_elapsed_time

## Purpose
Lets the batch density interface accept either one elapsed time for all satellites or a per-satellite vector.

## Design & Implementation
Two `@inline` methods: a `Float64` returns itself regardless of index, and a vector returns its `idx`-th element as `Float64`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `el_time` | Float64 | n/a | yes | Positional argument `el_time`. |
| in | `idx` | Int | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_batch_elapsed_time`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models_getdensitybatch_bang|getDensityBatch!]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:1010-1010`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/atmosphere/density_models.jl:840-840`
<!-- vulcan:connections:end -->

## Limitations
No bounds check on the vector form beyond what `_validate_density_batch_lengths` performed earlier.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 835.
