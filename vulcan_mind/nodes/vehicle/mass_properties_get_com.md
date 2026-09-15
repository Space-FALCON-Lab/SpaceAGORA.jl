---
id: vehicle.mass_properties_get_com
label: get_COM
kind: function
source:
  file: src/vehicle/structure/mass_properties.jl
  symbol: get_COM
  lines:
  - 6
  - 6
inputs:
- id: model
  type: SpacecraftModel
  units: n/a
  required: true
  description: Positional argument `model`.
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
  description: Return value of `get_COM`. Returns `get_COM(bodies)`.
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

# get_COM

## Purpose

`get_COM` computes the mass-weighted centre of mass of a set of links. One method takes a `SpacecraftModel` and forwards `model.links`; the other takes a `Vector{Link}` directly and does the arithmetic.

## Design & Implementation

The vector method accumulates into an `MVector{3,Float64}` initialised to zero and a scalar `total_mass`, looping over each `body` to add `body.r * body.m` and `body.m`, then returns `COM / total_mass`. Positions `body.r` are taken as-is in whatever frame the links store them, so no rotation is applied. The model method exists purely so callers can pass the whole model rather than extracting the link vector themselves.

## Theory & Math

The centre of mass returned is the standard mass-weighted mean of the link positions:

$$\mathbf{r}_{\text{COM}} = \frac{\sum_i m_i \mathbf{r}_i}{\sum_i m_i}$$

where $m_i$ is `body.m` in kilograms and $\mathbf{r}_i$ is `body.r` in metres for link $i$.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | SpacecraftModel | n/a | yes | Positional argument `model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `get_COM`. Returns `get_COM(bodies)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/structure/mass_properties.jl`

**Downstream**

- `callees` → [[vehx.structure_mass_properties_update_inertia_tensor_bang|update_inertia_tensor!]] · `callers` · call · `src/vehicle/structure/mass_properties.jl:27-27`
<!-- vulcan:connections:end -->

## Limitations

No guard exists for an empty `bodies` vector or for a total mass of zero, either of which produces `NaN` components rather than an error. Propellant mass tracked separately in `model.prop_mass` is ignored, so the result is a dry-structure centre of mass. The model method uses all of `model.links` regardless of which assembly root they belong to, so a multi-root model yields a blended centroid rather than a per-assembly one.

## Provenance
Mapped from `src/vehicle/structure/mass_properties.jl` line 6.
