---
id: core.reference_system__wrap_2pi
label: _wrap_2pi
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: _wrap_2pi
  lines:
  - 150
  - 150
inputs:
- id: theta
  type: Float64
  units: n/a
  required: true
  description: Positional argument `θ`.
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
  description: Return value of `_wrap_2pi`.
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

# _wrap_2pi

## Purpose
Reduces an angle into the half-open interval from zero to two pi, used to canonicalise orbital element angles.

## Design & Implementation
Applies `mod(θ, 2π)` and, defensively, adds `2π` if the result is negative. In Julia `mod` with a positive divisor already returns a non-negative result, so the second step is a guard against a differently behaving implementation. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `theta` | Float64 | n/a | yes | Positional argument `θ`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_wrap_2pi`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system__rvtoorbitalelement_core|_rvtoorbitalelement_core]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:189-189`
- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:516-516`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:516-516`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Not applied to inclination, which is correctly left in zero to pi; a caller wrapping inclination through this would produce a wrong value for retrograde orbits.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 150.
