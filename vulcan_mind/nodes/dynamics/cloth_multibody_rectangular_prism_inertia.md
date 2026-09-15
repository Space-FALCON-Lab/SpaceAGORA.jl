---
id: dynamics.cloth_multibody_rectangular_prism_inertia
label: rectangular_prism_inertia
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: rectangular_prism_inertia
  lines:
  - 187
  - 187
inputs:
- id: mass_kg
  type: Real
  units: n/a
  required: true
  description: Positional argument `mass_kg`.
- id: dimensions_m
  type: Any
  units: n/a
  required: true
  description: Positional argument `dimensions_m`.
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
  type: SMatrix{3,
  units: n/a
  description: Return value of `rectangular_prism_inertia`.
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

# rectangular_prism_inertia

## Purpose
Gives the body-frame inertia tensor of a uniform rectangular block, the building block for panel inertias.

## Theory & Math
$$
I = \frac{m}{12}\,\operatorname{diag}\left(d_y^2 + d_z^2,\; d_x^2 + d_z^2,\; d_x^2 + d_y^2\right)
$$

## Design & Implementation
Validates all three dimensions positive, then returns the diagonal tensor with entries `m/12 (d_j² + d_k²)`. `thin_panel_inertia` wraps it with a default thickness of one millimetre.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `mass_kg` | Real | n/a | yes | Positional argument `mass_kg`. |
| in | `dimensions_m` | Any | n/a | yes | Positional argument `dimensions_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `rectangular_prism_inertia`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:189-189`
- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:217-217`
- `callees` → [[dynamics.cloth_multibody_compliantjointactuator|CompliantJointActuator]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:249-249`
- `callees` → [[dynamics.cloth_multibody_complianttopologyedge|CompliantTopologyEdge]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:223-223`
- `callees` → [[dynamics.cloth_multibody_complianttopologynode|CompliantTopologyNode]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:201-201`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:217-217`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:217-217`
<!-- vulcan:connections:end -->

## Limitations
Principal axes are assumed aligned with the body frame; a panel modelled at an angle in its own frame needs the tensor rotated.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 187.
