---
id: flow.vehicle
label: Spacecraft model
kind: group
inputs:
- id: run_config
  type: SimulationConfiguration
  units: n/a
  description: Vehicle definitions in the dynamics model.
- id: actuation
  type: panel angles / wheel torques
  units: n/a
  description: Commands from GNC that change geometry or momentum.
  required: false
outputs:
- id: spacecraft_geometry
  type: SpacecraftModel
  units: n/a
  description: Links, joints, facets, thrusters, wheels, magnets and their attitudes.
tags:
- master-flow
charts:
- master
origin: agent
opens: vehicle
---

# Spacecraft model

## Purpose
The vehicle as the physics sees it: a tree of rigid links with dimensions, masses, inertias, reference areas and attitude quaternions, carrying components — thrusters, reaction wheels, magnets, facets — and the kinematics that place each link in the inertial frame.

## Design & Implementation
`src/vehicle/` defines `Link`, `Joint` and `SpacecraftModel`, the component structs, forward kinematics and geometry accessors (reference areas, lengths, normals and tangents), thermal properties, and the cloth-arm kinematics used by the robotics models. The standard example vehicle is a bus with two panel wings whose reference-area convention is chosen per scenario.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `run_config` | SimulationConfiguration | n/a | — | Vehicle definitions in the dynamics model. |
| in | `actuation` | panel angles / wheel torques | n/a | no | Commands from GNC that change geometry or momentum. |
| out | `spacecraft_geometry` | SpacecraftModel | n/a | — | Links, joints, facets, thrusters, wheels, magnets and their attitudes. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.configure|Configure a run]] · `run_config` → `run_config` · dataflow · `src/vehicle/spacecraft/model.jl`
- [[flow.gnc|Guidance, navigation & control]] · `actuation` → `actuation` · dataflow · `src/gnc/control/targeting_control.jl`

**Downstream**

- `spacecraft_geometry` → [[flow.forces|Force & torque models]] · `spacecraft_geometry` · dataflow · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
<!-- vulcan:connections:end -->

## Limitations
Links are mutable and hold live state such as thrust level and panel angle, so a `SpacecraftModel` cannot be shared between concurrently integrating runs without deep-copying; the panel-area convention is a documented trap that doubles drag when misused.
