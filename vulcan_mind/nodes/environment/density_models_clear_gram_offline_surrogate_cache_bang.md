---
id: environment.density_models_clear_gram_offline_surrogate_cache_bang
label: clear_gram_offline_surrogate_cache!
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: clear_gram_offline_surrogate_cache!
  lines:
  - 459
  - 459
inputs:
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
  description: Return value of `clear_gram_offline_surrogate_cache!`. Returns `_CLEAR_GRAM_OFFLINE_SURROGATE_CACHE_FN[]()`.
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

# clear_gram_offline_surrogate_cache!

## Purpose
Clears GRAM's loaded offline surrogate grids, forcing the next surrogate query to re-read its file.

## Design & Implementation
Calls through `_CLEAR_GRAM_OFFLINE_SURROGATE_CACHE_FN[]`, filled by the extension with `GRAMSuite.clear_gram_offline_surrogate_cache!` and defaulting to the not-loaded error.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `clear_gram_offline_surrogate_cache!`. Returns `_CLEAR_GRAM_OFFLINE_SURROGATE_CACHE_FN[]()`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[ext.spaceagoragramsuiteext___init__|__init__]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:25-25`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Clearing while a solve is using a surrogate forces a re-read on the next call, which takes the GRAM lock and stalls other threads.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 459.
