---
id: vehicle.model_spacecraftmodels
label: SpacecraftModels
kind: module
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: SpacecraftModels
  lines:
  - 1
  - 1
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
- vehicle
charts:
- vehicle
origin: agent
---

# SpacecraftModels

## Purpose
Module defining the vehicle data model: `Link` and `Joint` for multibody assemblies, `SpacecraftModel` aggregating them with mass and initial-condition metadata, the Keplerian and Cartesian initial-condition types, and the `DynamicsModel`, `GuidanceModel`, `NavigationModel` and `ControlModel` containers that bind effector tuples to spacecraft.

## Design & Implementation
`SpacecraftModels` uses `StaticArrays`, `LinearAlgebra`, the sibling `Components` module (for `ReactionWheelAssembly`, `Facet`, `Thruster`, `Magnet`) and `EphemeridesModels` (for `SpiceEphemeridesModel`, `ephemerides_time_seconds`, `planet_frame_lpi`). It defines constants `I3`, `DEFAULT_INITIAL_CONDITION_Q = (0,0,0,1)` and zero `DEFAULT_INITIAL_CONDITION_ANG_VEL`. `SpacecraftModel(; joints, links, root, ...)` appends the root to `links` if absent and sums `link.m` into `dry_mass`; `prop_mass`, `inertia_tensor`, `n_reaction_wheels`, `n_thrusters`, `initial_condition` and `id` are stored directly. The oblate `InitialCondition(planet; ra, hp)` path is supported by four private helpers for frame resolution, apsis direction, ellipsoid altitude and radius bisection.

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

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`SpacecraftModel.inertia_tensor` defaults to zeros and is never derived from the links, so callers must supply it. `dry_mass` is computed once at construction and not updated if links are modified later. The private oblate-altitude helper deviates from the standard geodetic formula, and the module duplicates that geometry code from the thruster guidance module instead of sharing it.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 1.
