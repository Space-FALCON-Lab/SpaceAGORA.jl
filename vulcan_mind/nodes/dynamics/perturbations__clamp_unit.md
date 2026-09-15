---
id: dynamics.perturbations__clamp_unit
label: _clamp_unit
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _clamp_unit
  lines:
  - 2301
  - 2301
inputs:
- id: x
  type: Float64
  units: n/a
  required: true
  description: Positional argument `x`.
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
  description: Return value of `_clamp_unit`. Returns `clamp(x, -1.0, 1.0)`.
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

# _clamp_unit

## Purpose
Clamps a value into the closed unit interval so `asin` and `acos` in the eclipse geometry never receive an argument a hair outside their domain.

## Design & Implementation
A local `@inline` closure inside `eclipse_area_calc` returning `clamp(x, -1, 1)`. Applied to every apparent-radius and separation ratio.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Float64 | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_clamp_unit`. Returns `clamp(x, -1.0, 1.0)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Clamping masks a genuinely invalid geometry, such as a satellite inside the planet, by returning an edge angle instead of failing.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 2301.
