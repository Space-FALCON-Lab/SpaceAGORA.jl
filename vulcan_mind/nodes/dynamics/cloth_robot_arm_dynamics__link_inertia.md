---
id: dynamics.cloth_robot_arm_dynamics__link_inertia
label: _link_inertia
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: _link_inertia
  lines:
  - 86
  - 86
inputs:
- id: link
  type: Any
  units: n/a
  required: true
  description: Positional argument `link`.
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
  type: SMatrix
  units: n/a
  description: Return value of `_link_inertia`. Returns `SMatrix{3, 3, Float64}(`.
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

# _link_inertia

## Purpose
Approximates each arm link as a uniform solid cylinder and returns its body-frame inertia tensor for the compliant multibody model.

## Theory & Math
$I = \operatorname{diag}\left(\tfrac{1}{2} m r^2,\; \tfrac{1}{12} m (3 r^2 + \ell^2),\; \tfrac{1}{12} m (3 r^2 + \ell^2)\right)$ with mass $m$ (kg), radius $r$ (m), and length $\ell = \lVert \mathbf{v}_{parent} \rVert$ (m).

## Design & Implementation
Reads `len = norm(link.vector_parent)`, `r = link.radius_m`, and `m = link.mass_kg`. Returns a diagonal `SMatrix{3,3,Float64}` with `0.5*m*r²` about the x axis and `(1/12)*m*(3r² + len²)` about y and z. The link's own axis is therefore assumed to be body x.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `link` | Any | n/a | yes | Positional argument `link`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix | n/a | — | Return value of `_link_inertia`. Returns `SMatrix{3, 3, Float64}(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_multibody|cloth_robot_arm_multibody]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:235-235`
- [[dynx.multibody_cloth_assign_coupled_cloth_robot_arm_rhs_bang|assign_coupled_cloth_robot_arm_rhs!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:416-416`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Despite the docstring, no explicit inertia field on `link` is ever consulted. If a link's `vector_parent` is not along body x the tensor is misaligned with the actual geometry. Zero mass or radius produces a singular tensor, and `J \ (...)` in the RHS then fails or returns Inf.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 86.
