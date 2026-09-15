---
id: environment.density_models_clear_gram_static_grid_cache_bang
label: clear_gram_static_grid_cache!
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: clear_gram_static_grid_cache!
  lines:
  - 458
  - 458
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
  description: Return value of `clear_gram_static_grid_cache!`. Returns `_CLEAR_GRAM_STATIC_GRID_CACHE_FN[]()`.
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

# clear_gram_static_grid_cache!

## Purpose
Clears GRAM's process-wide static-grid cache so memory can be released between campaigns or a rebuilt grid can take effect.

## Design & Implementation
A one-line call through `_CLEAR_GRAM_STATIC_GRID_CACHE_FN[]`, a `Ref{Function}` slot that defaults to raising the not-loaded error and is overwritten by the extension's `__init__` with `GRAMSuite.clear_gram_static_grid_cache!`. The slot pattern lets the core package export the function name without depending on GRAMSuite.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `clear_gram_static_grid_cache!`. Returns `_CLEAR_GRAM_STATIC_GRID_CACHE_FN[]()`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[ext.spaceagoragramsuiteext___init__|__init__]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:24-24`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the target is resolved through a mutable reference at call time, a call site cannot know statically whether the operation is available, and clearing while another thread is mid-evaluation is only safe if the extension's implementation takes the GRAM lock.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 458.
