---
id: dynamics.torque_models_body_angular_velocity
label: body_angular_velocity
kind: function
source:
  file: src/dynamics/rotational/torque_models.jl
  symbol: body_angular_velocity
  lines:
  - 5
  - 5
inputs:
- id: omega
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `omega`.
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
  description: Return value of `body_angular_velocity`.
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

# body_angular_velocity

## Purpose
Normalises an arbitrary real-valued angular velocity vector into the fixed-size, double-precision form the rotational dynamics expect.

## Design & Implementation
Reads the first three elements of `omega`, converts each through `Float64`, and returns an `SVector{3,Float64}`. Marked `@inline` with a declared return type, so on the integration hot path it lowers to three loads and no heap allocation, whatever concrete `AbstractVector` the caller passed.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `omega` | AbstractVector{<:Real} | n/a | yes | Positional argument `omega`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `body_angular_velocity`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/rotational/torque_models.jl`
- [[simulation.dynamics_rhs__assign_orientation_rhs_bang|_assign_orientation_rhs!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1594-1594`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/rotational/torque_models.jl:6-6`
<!-- vulcan:connections:end -->

## Limitations
It indexes elements one through three without checking length, so a shorter vector raises a bounds error and a longer one is silently truncated; unit consistency is the caller's responsibility.

## Provenance
Mapped from `src/dynamics/rotational/torque_models.jl` line 5.
