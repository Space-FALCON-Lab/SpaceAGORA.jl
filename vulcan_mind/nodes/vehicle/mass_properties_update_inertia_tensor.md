---
id: vehicle.mass_properties_update_inertia_tensor
label: update_inertia_tensor
kind: function
source:
  file: src/vehicle/structure/mass_properties.jl
  symbol: update_inertia_tensor
  lines:
  - 53
  - 53
inputs:
- id: bodies
  type: Vector{Link}
  units: n/a
  required: true
  description: Positional argument `bodies`.
- id: prop_mass
  type: Float64
  units: n/a
  required: false
  description: Positional argument `prop_mass` (default `0.0`).
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
  description: Return value of `update_inertia_tensor`. Returns `inertia_tensor`.
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

# update_inertia_tensor

## Purpose

`update_inertia_tensor(bodies::Vector{Link}, prop_mass::Float64 = 0.0)` computes the combined inertia tensor of a set of links about the assembly origin, applying the parallel axis theorem to each link and optionally adding shared propellant mass to the root link's offset term.

## Design & Implementation

Starting from a zero `SMatrix{3,3,Float64}`, the loop builds `R` as the identity for a root link and `rot(b.q)` otherwise, rotates the link's own tensor with the similarity transform `R * b.inertia * R'`, and adds the offset term `(b.m + fuel_mass) * hat(r) * hat(r)'`, where `hat` is the skew-symmetric cross-product matrix of the link position `b.r` and `fuel_mass` is `prop_mass` only for the root link. The related `update_inertia_tensor!` wraps this: it calls `traverse_bodies`, passes `model.prop_mass[root_index]`, then either `push!`es the result onto `model.inertia_tensors` or overwrites the entry at `root_index`.

## Theory & Math

For each link $i$ with body-frame inertia $I_i$, attitude rotation $R_i$, position $\mathbf{r}_i$ and effective mass $m_i^{\text{eff}} = m_i + m_{\text{prop}}\delta_{i,\text{root}}$, the accumulated tensor is

$$I = \sum_i \left( R_i I_i R_i^{\mathsf{T}} + m_i^{\text{eff}} \, [\mathbf{r}_i]_\times [\mathbf{r}_i]_\times^{\mathsf{T}} \right)$$

where $[\mathbf{r}]_\times$ is the skew-symmetric matrix such that $[\mathbf{r}]_\times \mathbf{v} = \mathbf{r} \times \mathbf{v}$, and $[\mathbf{r}]_\times [\mathbf{r}]_\times^{\mathsf{T}} = \lVert\mathbf{r}\rVert^2 I_{3\times3} - \mathbf{r}\mathbf{r}^{\mathsf{T}}$ is the usual parallel-axis displacement term. Units are $\text{kg}\,\text{m}^2$.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `bodies` | Vector{Link} | n/a | yes | Positional argument `bodies`. |
| in | `prop_mass` | Float64 | n/a | no | Positional argument `prop_mass` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `update_inertia_tensor`. Returns `inertia_tensor`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[vehx.structure_mass_properties_update_inertia_tensor_bang|update_inertia_tensor!]] · `callees` → `callers` · call · `src/vehicle/structure/mass_properties.jl:36-36`

**Downstream**

- `callees` → [[core.quaternion_utils_hat|hat]] · `callers` · call · `src/vehicle/structure/mass_properties.jl:61-61`
- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/vehicle/structure/mass_properties.jl:56-56`
- `callees` → [[vehicle.mass_properties_get_inertia_tensor|get_inertia_tensor]] · `callers` · call · `src/vehicle/structure/mass_properties.jl:67-67`
<!-- vulcan:connections:end -->

## Limitations

The tensor is taken about the frame origin the link positions are expressed in, not about the assembly centre of mass, so callers wanting a COM-referenced tensor must shift it themselves using `get_COM`. Propellant is lumped as a point mass at the root link's position, an approximation that ignores tank geometry and slosh. The code comment marks this as a legacy structure model, and the mutating wrapper's `push!` branch assumes root indices arrive in increasing order, otherwise the appended tensor lands at the wrong index.

## Provenance
Mapped from `src/vehicle/structure/mass_properties.jl` line 53.
