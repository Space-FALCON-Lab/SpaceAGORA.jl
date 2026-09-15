---
id: vehicle.assembly_add_thruster_bang
label: add_thruster!
kind: function
source:
  file: src/vehicle/spacecraft/assembly.jl
  symbol: add_thruster!
  lines:
  - 74
  - 74
inputs:
- id: model
  type: SpacecraftModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: link
  type: Link
  units: n/a
  required: true
  description: Positional argument `link`.
- id: thruster
  type: Thruster
  units: n/a
  required: true
  description: Positional argument `thruster`.
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
  description: Return value of `add_thruster!`; mutates `model` in place.
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

# add_thruster!

## Purpose

`add_thruster!(model::SpacecraftModel, link::Link, thruster::Thruster)` mounts one thruster on a link and incrementally extends the link's thrust-torque Jacobian, so that the control allocation layer can map commanded thruster levels to body torque.

## Design & Implementation

The function first calls `normalize!(thruster.direction)`, mutating the supplied thruster so its direction is a unit vector. If `link.thrusters` is empty it sets `link.J_thruster = cross(thruster.location, thruster.direction)`, the single r x F column; otherwise it appends that column with `hcat(link.J_thruster, cross(...))`, growing the Jacobian to 3 x n. The thruster is then pushed onto `link.thrusters` and `model.n_thrusters` is incremented, so three objects are mutated per call: the thruster, the link and the model.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | SpacecraftModel | n/a | yes | Positional argument `model`. |
| in | `link` | Link | n/a | yes | Positional argument `link`. |
| in | `thruster` | Thruster | n/a | yes | Positional argument `thruster`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `add_thruster!`; mutates `model` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/assembly.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/vehicle/spacecraft/assembly.jl:90-90`
<!-- vulcan:connections:end -->

## Theory & Math
Each column of the Jacobian is the unit-thrust torque of one nozzle about the link origin,

$$\mathbf{J}_i = \mathbf{r}_i \times \hat{\mathbf{u}}_i,$$

with $\mathbf{r}_i$ the mount location and $\hat{\mathbf{u}}_i$ the normalised thrust direction. The resulting body torque for a thrust magnitude vector $\mathbf{f}$ is $\boldsymbol{\tau} = \mathbf{J}\mathbf{f}$.

## Limitations

The Jacobian is grown by `hcat`, which reallocates the whole matrix on every call, making assembly of a large thruster cluster quadratic in the number of thrusters. Normalising in place silently rewrites the caller's `Thruster` and fails on a zero direction vector. The r x F column ignores the thruster's own plume geometry and any offset between link origin and centre of mass, and `model.n_thrusters` is only ever incremented, so removal is unsupported.

## Provenance
Mapped from `src/vehicle/spacecraft/assembly.jl` line 74.
