---
id: dynamics.perturbations__fully_normalized_legendre_scale
label: _fully_normalized_legendre_scale
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _fully_normalized_legendre_scale
  lines:
  - 169
  - 169
inputs:
- id: l
  type: Int
  units: n/a
  required: true
  description: Positional argument `l`.
- id: m
  type: Int
  units: n/a
  required: true
  description: Positional argument `m`.
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
  description: Return value of `_fully_normalized_legendre_scale`.
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

# _fully_normalized_legendre_scale

## Purpose
Computes the factor relating an unnormalised associated Legendre coefficient to its fully normalised counterpart, for one degree and order.

## Theory & Math
$$
N_{lm} = \sqrt{(2 - \delta_{0m})(2l + 1)\,\frac{(l - m)!}{(l + m)!}}
$$

## Design & Implementation
Validates `l >= 0` and `0 <= m <= l`, then returns `sqrt((2 - δ_0m)(2l + 1) (l - m)! / (l + m)!)` with the factorial ratio computed via `loggamma` to avoid overflow at high degree. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `l` | Int | n/a | yes | Positional argument `l`. |
| in | `m` | Int | n/a | yes | Positional argument `m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_fully_normalized_legendre_scale`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__convert_harmonics_coefficients_to_full_bang|_convert_harmonics_coefficients_to_full!]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:194-194`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Accurate to floating-point roundoff of the log-gamma difference, which for degree 165 and order near 165 is a few ulps but not exact.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 169.
