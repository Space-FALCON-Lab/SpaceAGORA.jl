---
id: dynamics.perturbations__convert_harmonics_coefficients_to_full_bang
label: _convert_harmonics_coefficients_to_full!
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _convert_harmonics_coefficients_to_full!
  lines:
  - 177
  - 177
inputs:
- id: C
  type: AbstractMatrix{Float64}
  units: n/a
  required: true
  description: Positional argument `C`.
- id: S
  type: AbstractMatrix{Float64}
  units: n/a
  required: true
  description: Positional argument `S`.
- id: L
  type: Int
  units: n/a
  required: true
  description: Positional argument `L`.
- id: M
  type: Int
  units: n/a
  required: true
  description: Positional argument `M`.
- id: normalization
  type: Symbol
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
  type: Nothing
  units: n/a
  description: Return value of `_convert_harmonics_coefficients_to_full!`; mutates
    `C` in place.
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

# _convert_harmonics_coefficients_to_full!

## Purpose
Rescales a loaded coefficient set from Schmidt or unnormalised convention into fully normalised form in place, so the evaluator only ever sees one convention.

## Design & Implementation
Canonicalises the normalisation and returns immediately for `:full`. Otherwise for every degree `l` and order `m` up to `min(M, l)` it divides `C` and `S` by `sqrt(2l + 1)` for Schmidt or by `_fully_normalized_legendre_scale(l, m)` for unnormalised. `@inline`, mutates both matrices.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `C` | AbstractMatrix{Float64} | n/a | yes | Positional argument `C`. |
| in | `S` | AbstractMatrix{Float64} | n/a | yes | Positional argument `S`. |
| in | `L` | Int | n/a | yes | Positional argument `L`. |
| in | `M` | Int | n/a | yes | Positional argument `M`. |
| in | `normalization` | Symbol | n/a | yes | Positional argument `normalization`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_convert_harmonics_coefficients_to_full!`; mutates `C` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__harmonics_model_cache_key|_harmonics_model_cache_key]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:835-835`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.perturbations__canonical_harmonics_normalization|_canonical_harmonics_normalization]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:184-184`
- `callees` → [[dynamics.perturbations__fully_normalized_legendre_scale|_fully_normalized_legendre_scale]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:194-194`
<!-- vulcan:connections:end -->

## Limitations
Applies the scale once, so calling it twice on the same matrices double-converts; the model constructor is the only caller.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 177.
