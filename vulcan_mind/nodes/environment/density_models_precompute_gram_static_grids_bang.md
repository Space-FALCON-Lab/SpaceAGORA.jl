---
id: environment.density_models_precompute_gram_static_grids_bang
label: precompute_gram_static_grids!
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: precompute_gram_static_grids!
  lines:
  - 445
  - 445
inputs:
- id: kwargs
  type: Vararg{Any}
  units: n/a
  required: false
  description: Keyword argument `kwargs` (variadic).
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
  description: Return value of `precompute_gram_static_grids!`; mutates `kwargs` in
    place.
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

# precompute_gram_static_grids!

## Purpose
The core-package stub for static-grid precomputation, replaced by the extension's real implementation when GRAMSuite loads.

## Design & Implementation
A method on `AbstractDensityModel` accepting any keywords that immediately raises the not-loaded error. The extension defines a more specific method on `GRAMAtmosphereModel` that forwards to `GRAMSuite.precompute_gram_static_grids!` with the process-wide lock.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `precompute_gram_static_grids!`; mutates `kwargs` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[ext.spaceagoragramsuiteext__gram_spice_ephemeris_state|_gram_spice_ephemeris_state]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:149-149`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.density_models__gram_not_loaded_error|_gram_not_loaded_error]] · `callers` · call · `src/environment/atmosphere/density_models.jl:446-446`
<!-- vulcan:connections:end -->

## Limitations
Calling it on a non-GRAM model with the extension loaded still hits this stub and reports GRAM as not loaded, which is misleading.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 445.
