---
id: dynamics.aerodynamic_wrench_models__attitude_link_alpha
label: _attitude_link_alpha
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _attitude_link_alpha
  lines:
  - 260
  - 260
inputs:
- id: body
  type: Any
  units: n/a
  required: true
  description: Positional argument `body`.
- id: root
  type: Any
  units: n/a
  required: true
  description: Positional argument `root`.
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
  description: Return value of `_attitude_link_alpha`.
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

# _attitude_link_alpha

## Purpose
`:attitude`-mode incidence: composes the child quaternion with the root attitude so a rigidly mounted child shares the root's incidence relative to a flow-aligned frame.

## Design & Implementation
With `e1 = (1,0,0)`, computes `flow_body = body.root ? rot(body.q) * e1 : rot(body.q) * (rot(root.q) * e1)` and returns `atan(flow_body[1], flow_body[3])`, matching the `orientation_sim=true` composition `rot(child.q) * rot(root.q) * v`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `body` | Any | n/a | yes | Positional argument `body`. |
| in | `root` | Any | n/a | yes | Positional argument `root`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_attitude_link_alpha`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__aero_link_angles|_aero_link_angles]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:312-312`
- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:795-795`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:262-262`
<!-- vulcan:connections:end -->

## Limitations
Sideslip is always zero in this mode, so yaw-only attitudes (flow into body ±y) collapse to the `:max_drag` result. The attitude is held relative to the velocity direction, not inertially fixed, which is a modelling choice the docstring documents.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 260.
