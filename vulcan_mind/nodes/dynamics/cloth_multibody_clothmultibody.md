---
id: dynamics.cloth_multibody_clothmultibody
label: ClothMultibody
kind: module
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: ClothMultibody
  lines:
  - 2
  - 2
inputs:
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
  type: Any
  units: n/a
  description: Value produced by this symbol.
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

# ClothMultibody

## Purpose
The compliant multibody module modelling cloth-like panels as rigid bodies connected by six-axis spring-damper joints, with topology builders, joint-load diagnostics and two fixed-step integrators.

## Design & Implementation
Depends on LinearAlgebra and StaticArrays only. It defines value types for bodies, joints, models, trajectories, topology nodes and edges, actuators and load diagnostics, and exports builders (`build_compliant_topology`, `build_rectangular_compliant_grid`), state packing (`compliant_state_vector`, `compliant_state_parts`), the derivative (`compliant_multibody_dynamics`), diagnostics (`compliant_joint_loads`), RK4 and implicit-midpoint steppers, and `simulate_compliant_multibody`. State is a flat vector of thirteen entries per body: position, scalar-last quaternion, velocity, body angular rate. Quaternion helpers are private and normalise defensively.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The model is a maximal-coordinate penalty formulation — joints are stiff springs, not constraints — so stiffness sets the stable step size and the RK4 stepper is only usable at small `dt`; the implicit stepper compensates with a dense finite-difference Jacobian whose cost is quadratic in body count.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 2.
