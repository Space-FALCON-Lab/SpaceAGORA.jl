---
id: flow.forces
label: Force & torque models
kind: group
inputs:
- id: force_requests
  type: StateSample + EnvironmentSample
  units: n/a
  description: From the RHS.
- id: spacecraft_geometry
  type: SpacecraftModel
  units: n/a
  description: Links, facets, areas and attitudes.
outputs:
- id: wrench
  type: (force_ii, torque_body)
  units: n/a
  description: Inertial force and body torque per effector.
tags:
- master-flow
charts:
- master
origin: agent
opens: dynamics
---

# Force & torque models

## Purpose
The physics: point-mass, J2 and spherical-harmonics gravity, third-body perturbations, solar radiation pressure with eclipse, planetary albedo and infrared, free-molecular aerodynamics over the vehicle's links, magnetic torque and eddy damping, and the compliant multibody dynamics for cloth and robot arms.

## Design & Implementation
Each model is an `AbstractForceTorqueModel` in `src/dynamics/` declaring its environment requirements and its solver partition (implicit or explicit) and gravity-backbone role; the harmonics model carries precomputed recurrence tables and per-satellite scratch workspaces with a batched Pines kernel for constellations. Aerodynamic models integrate Hart free-molecular coefficients over each link at its incidence.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `force_requests` | StateSample + EnvironmentSample | n/a | — | From the RHS. |
| in | `spacecraft_geometry` | SpacecraftModel | n/a | — | Links, facets, areas and attitudes. |
| out | `wrench` | (force_ii, torque_body) | n/a | — | Inertial force and body torque per effector. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.environment|Environment sampling]] · `frames_and_bodies` → `force_requests` · dataflow · `src/simulation/engine/effector_sampling.jl`
- [[flow.rhs|Dynamics right-hand side]] · `force_requests` → `force_requests` · dataflow · `src/simulation/engine/dynamics_rhs.jl`
- [[flow.vehicle|Spacecraft model]] · `spacecraft_geometry` → `spacecraft_geometry` · dataflow · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Cannonball radiation models ignore attitude; the harmonics field is treated as static within an integrator step; and the flat-path vectorised kernels exist only for four effector types, so anything else runs per satellite.
