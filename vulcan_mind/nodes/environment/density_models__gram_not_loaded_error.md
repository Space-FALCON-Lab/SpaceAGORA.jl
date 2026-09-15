---
id: environment.density_models__gram_not_loaded_error
label: _gram_not_loaded_error
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _gram_not_loaded_error
  lines:
  - 438
  - 438
inputs:
- id: fn_name
  type: String
  units: n/a
  required: true
  description: Positional argument `fn_name`.
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
  description: Return value of `_gram_not_loaded_error`.
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

# _gram_not_loaded_error

## Purpose
The uniform failure raised by every GRAM entry point when GRAMSuite has not been loaded.

## Design & Implementation
Calls `error` with a message naming the attempted function and instructing the user to add and load GRAMSuite. Used as the body of the core-package stubs for precompute, cache clearing and density evaluation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `fn_name` | String | n/a | yes | Positional argument `fn_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gram_not_loaded_error`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__gram_core_density_state|_gram_core_density_state]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:471-471`
- [[environment.density_models_precompute_gram_static_grids_bang|precompute_gram_static_grids!]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:446-446`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It raises `ErrorException` rather than a dedicated exception type, so callers cannot catch the not-loaded case specifically.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 438.
