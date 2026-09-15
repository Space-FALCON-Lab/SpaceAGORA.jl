---
id: dynamics.aerodynamic_wrench_models_collect_and_reset_link_wrenches_bang
label: collect_and_reset_link_wrenches!
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: collect_and_reset_link_wrenches!
  lines:
  - 211
  - 211
inputs:
- id: bodies
  type: Any
  units: n/a
  required: true
  description: Positional argument `bodies`.
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
  description: Return value of `collect_and_reset_link_wrenches!`; mutates `bodies`
    in place. Returns `SVector{3, Float64}(force_acc), SVector{3, Float64}(torque_acc)`.
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

# collect_and_reset_link_wrenches!

## Purpose
Sums each link's accumulated `net_force` and `net_torque`, zeroes them in place, and returns the totals as static vectors.

## Design & Implementation
Uses `MVector{3,Float64}` accumulators, iterates `bodies` with `@inbounds`, adds `SVector(b.net_force)` and `SVector(b.net_torque)`, then broadcasts zero into both link fields. Returns `(SVector(force_acc), SVector(torque_acc))`. Fresh accumulators avoid aliasing when there is only one link.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `bodies` | Any | n/a | yes | Positional argument `bodies`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `collect_and_reset_link_wrenches!`; mutates `bodies` in place. Returns `SVector{3, Float64}(force_acc), SVector{3, Float64}(torque_acc)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:694-694`
- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:999-999`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Mixes frames: `net_force` is inertial while `net_torque` is body-frame per the legacy convention, and the sum of body-frame torques across differently oriented links is not physically meaningful. Only used by the broken legacy `calcForceTorque` methods.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 211.
