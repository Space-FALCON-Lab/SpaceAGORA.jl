---
id: dynamics.torque_models_body_torque
label: body_torque
kind: function
source:
  file: src/dynamics/rotational/torque_models.jl
  symbol: body_torque
  lines:
  - 1
  - 1
inputs:
- id: torque
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `torque`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `body_torque`.
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

# body_torque

## Purpose
Normalises an arbitrary real-valued torque vector into the fixed-size, double-precision form the rotational dynamics expect.

## Design & Implementation
Identical in shape to the angular-velocity conversion: three `Float64` conversions into an `SVector{3,Float64}`, `@inline` with a declared return type so it allocates nothing on the hot path. Having a distinct function for torque rather than reusing one generic converter keeps the call sites self-describing about which physical quantity is being marshalled.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `torque` | AbstractVector{<:Real} | n/a | yes | Positional argument `torque`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `body_torque`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/rotational/torque_models.jl`
- [[simulation.dynamics_rhs__assign_orientation_rhs_bang|_assign_orientation_rhs!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1595-1595`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/rotational/torque_models.jl:2-2`
<!-- vulcan:connections:end -->

## Limitations
No length check and no unit checking; a caller passing a body-frame torque where an inertial one is expected gets no diagnostic from this layer.

## Provenance
Mapped from `src/dynamics/rotational/torque_models.jl` line 1.
