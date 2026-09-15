---
id: dynamics.cloth_multibody_compliantbody
label: CompliantBody
kind: struct
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: CompliantBody
  lines:
  - 20
  - 20
inputs:
- id: name
  type: Symbol
  units: n/a
  required: true
  description: Field `name`.
- id: mass_kg
  type: Float64
  units: n/a
  required: true
  description: Field `mass_kg`.
- id: inertia_body_kg_m2
  type: SMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Field `inertia_body_kg_m2`.
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
  type: CompliantBody
  units: n/a
  description: Constructed `CompliantBody`.
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

# CompliantBody

## Purpose
One rigid panel in the compliant model: its name, mass and body-frame inertia tensor.

## Design & Implementation
An immutable struct with `name`, `mass_kg` and a three-by-three static `inertia_body_kg_m2`. Constructed by `build_compliant_topology` from topology nodes; typically the inertia comes from `thin_panel_inertia` or `rectangular_prism_inertia`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | Symbol | n/a | yes | Field `name`. |
| in | `mass_kg` | Float64 | n/a | yes | Field `mass_kg`. |
| in | `inertia_body_kg_m2` | SMatrix{3, 3, Float64} | n/a | yes | Field `inertia_body_kg_m2`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CompliantBody | n/a | — | Constructed `CompliantBody`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_build_compliant_topology|build_compliant_topology]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:291-291`
- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_multibody|cloth_robot_arm_multibody]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:235-235`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No validation that the inertia is symmetric positive definite; a mis-entered tensor produces a singular `J \` solve in the dynamics rather than a construction error.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 20.
