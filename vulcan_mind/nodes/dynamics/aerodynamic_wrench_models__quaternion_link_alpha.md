---
id: dynamics.aerodynamic_wrench_models__quaternion_link_alpha
label: _quaternion_link_alpha
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _quaternion_link_alpha
  lines:
  - 251
  - 251
inputs:
- id: body
  type: Any
  units: n/a
  required: true
  description: Positional argument `body`.
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
  description: Return value of `_quaternion_link_alpha`.
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

# _quaternion_link_alpha

## Purpose
Historical `:max_drag` incidence for a non-root link: the angle of the assumed flow direction (reference +x) expressed in the link frame.

## Design & Implementation
Computes `flow_body = rot(body.q) * (1, 0, 0)` and returns `atan(flow_body[1], flow_body[3])`. Deliberately does not compose with the root attitude, preserving bit-identical legacy behaviour as the comment insists.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `body` | Any | n/a | yes | Positional argument `body`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_quaternion_link_alpha`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__aero_link_angles|_aero_link_angles]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:316-316`
- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:802-802`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:252-252`
<!-- vulcan:connections:end -->

## Limitations
Ignores that child quaternions are root-relative, so a child mounted on a rotated root reports the wrong incidence. Uses the reference-to-body reading of `rot(q)`; the transpose interpretation would flip lift sign.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 251.
