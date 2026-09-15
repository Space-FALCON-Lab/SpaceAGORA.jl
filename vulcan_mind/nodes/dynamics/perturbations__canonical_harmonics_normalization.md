---
id: dynamics.perturbations__canonical_harmonics_normalization
label: _canonical_harmonics_normalization
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _canonical_harmonics_normalization
  lines:
  - 154
  - 154
inputs:
- id: normalization
  type: Any
  units: n/a
  required: true
  description: Positional argument `normalization`.
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
  type: Symbol
  units: n/a
  description: Return value of `_canonical_harmonics_normalization`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _canonical_harmonics_normalization

## Purpose
Maps the many spellings of a harmonics normalisation convention onto the three symbols the coefficient converter understands.

## Design & Implementation
Converts to `Symbol` and returns `:full` for `full` or `fully_normalized`, `:schmidt` for `schmidt` and its two quasi-normalised spellings, and `:unnormalized` for `unnormalized` or `none`, raising `ArgumentError` listing the accepted values otherwise. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `normalization` | Any | n/a | yes | Positional argument `normalization`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_canonical_harmonics_normalization`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__convert_harmonics_coefficients_to_full_bang|_convert_harmonics_coefficients_to_full!]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:184-184`
- [[dynamics.perturbations__harmonics_model_cache_key|_harmonics_model_cache_key]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:725-725`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Case-sensitive; `:Full` is rejected.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 154.
