---
id: core.reference_system__safe_acos
label: _safe_acos
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: _safe_acos
  lines:
  - 155
  - 155
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
  type: Float64
  units: n/a
  description: Return value of `_safe_acos`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# _safe_acos

## Purpose
Evaluates the arccosine after clamping its argument to the valid domain, so rounding in a dot product cannot produce `NaN`.

## Design & Implementation
Clamps `x` into `[-1, 1]` and calls `acos`. `@inline` with a `::Float64` return. Used for inclination, argument of periapsis and true anomaly where normalised dot products can exceed one by a few ulps.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Float64 | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_safe_acos`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system__rvtoorbitalelement_core|_rvtoorbitalelement_core]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:179-179`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Clamping hides a genuinely out-of-range input, such as a non-normalised vector, by silently returning zero or pi instead of failing.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 155.
