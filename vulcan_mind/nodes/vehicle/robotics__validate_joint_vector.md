---
id: vehicle.robotics__validate_joint_vector
label: _validate_joint_vector
kind: function
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: _validate_joint_vector
  lines:
  - 161
  - 161
inputs:
- id: model
  type: ClothArmModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: q
  type: Any
  units: n/a
  required: true
  description: Positional argument `q`.
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
  type: Float64.
  units: n/a
  description: Return value of `_validate_joint_vector`. Returns `Float64.(collect(q))`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# _validate_joint_vector

## Purpose
Checks a caller's joint vector against the model's joint count and converts it to a fresh `Vector{Float64}`.

## Design & Implementation
Compares `length(q)` with the number of joints, throwing `ArgumentError` that reports both numbers, then returns `Float64.(collect(q))`. The `collect` accepts tuples and ranges; the broadcast conversion returns a new vector so callers such as `cloth_ik` can mutate it without aliasing the input.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64. | n/a | — | Return value of `_validate_joint_vector`. Returns `Float64.(collect(q))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[vehicle.cloth_fk|cloth_fk]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:169-169`
- [[vehicle.robotics_cloth_fk_state|cloth_fk_state]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:218-218`
- [[vehicle.robotics_cloth_ik|cloth_ik]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:271-271`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It validates length only, not finiteness or joint limits, so NaN angles pass through to FK.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 161.
