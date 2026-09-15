---
id: simulation.runtime__extract_pos_vel
label: _extract_pos_vel
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/runtime.jl
  symbol: _extract_pos_vel
  lines:
  - 1
  - 1
inputs:
- id: x
  type: Any
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
  type: SVector
  units: n/a
  description: Return value of `_extract_pos_vel`. Returns `SVector{3, Float64}(x[1],
    x[2], x[3]), SVector{3, Float64}(x[4], x[5], x[6])`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _extract_pos_vel

## Purpose
Splits a per-spacecraft state vector into its inertial position and velocity as fixed-size static vectors for the density callback.

## Design & Implementation
Reads elements one to three into an `SVector{3,Float64}` for position and four to six for velocity, returning them as a tuple. Declared `@inline` so the callback's per-satellite loop lowers to six loads with no allocation regardless of the concrete state container type.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `_extract_pos_vel`. Returns `SVector{3, Float64}(x[1], x[2], x[3]), SVector{3, Float64}(x[4], x[5], x[6])`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/runtime.jl`
- [[simulation.runtime__stage_environment_kinematics|_stage_environment_kinematics]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:58-58`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No length check and no unit check; a state shorter than six elements raises a bounds error, and any layout other than position-then-velocity in metres and metres per second is silently misread.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/runtime.jl` line 1.
