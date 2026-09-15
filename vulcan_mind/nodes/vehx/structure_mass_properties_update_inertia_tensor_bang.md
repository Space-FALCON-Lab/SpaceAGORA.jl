---
id: vehx.structure_mass_properties_update_inertia_tensor_bang
label: update_inertia_tensor!
kind: function
source:
  file: src/vehicle/structure/mass_properties.jl
  symbol: update_inertia_tensor!
  lines:
  - 32
  - 46
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Vehicle package namespace that loads this file and exports the symbol.
- id: model
  type: SpacecraftModel
  units: n/a
  required: true
  description: Model whose stored inertia tensor for the relevant assembly is refreshed.
- id: body
  type: Link
  units: n/a
  required: true
  description: Any link of the assembly whose composite inertia is to be recomputed.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: inertia_tensor
  type: SMatrix{3,3,Float64}
  units: kg*m^2
  description: Composite inertia tensor of the assembly about its centre of mass,
    also cached on the model.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
- structure
- mass
charts:
- vehx
origin: agent
---

# update_inertia_tensor!

## Purpose
`update_inertia_tensor!` refreshes the composite inertia of an assembly and stores it back on the model. Because articulation, propellant depletion and deployment all change the mass distribution, the attitude dynamics cannot rely on a value fixed at build time; this routine is the sanctioned way to bring the cached tensor back in step with the current configuration.

## Theory & Math
For links $k$ with mass $m_k$, body inertia $I_k$ and offset $\mathbf{d}_k$ from the assembly centre of mass, the parallel-axis theorem gives $I = \sum_k \left[ C_k I_k C_k^{\mathsf T} + m_k\left((\mathbf{d}_k^{\mathsf T}\mathbf{d}_k) \mathbf{1}_3 - \mathbf{d}_k \mathbf{d}_k^{\mathsf T}\right)\right]$, where $C_k$ rotates link axes into assembly axes.

## Model & Assumptions
The composite tensor is the sum over every link of its own body-frame inertia rotated into the assembly frame plus the parallel-axis contribution of its mass at its offset from the assembly centre of mass. Propellant is treated as an additional mass associated with the root, passed through to the pure computation as a scalar. Links are assumed rigid and their individual tensors are assumed expressed about their own centres of mass.

## Design & Implementation
The mutating method is deliberately thin. It calls `traverse_bodies` to obtain the assembly and its root index, forwards the body list and the propellant mass to the non-mutating `update_inertia_tensor`, and then writes the result into `model.inertia_tensors`. The write is guarded: if the root index exceeds the current length of the cache vector the tensor is pushed, otherwise the existing slot is overwritten. That grow-or-replace pattern lets assemblies be registered in any order. The computed tensor is also returned, so a caller can use it directly without a second lookup.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Vehicle package namespace that loads this file and exports the symbol. |
| in | `model` | SpacecraftModel | n/a | yes | Model whose stored inertia tensor for the relevant assembly is refreshed. |
| in | `body` | Link | n/a | yes | Any link of the assembly whose composite inertia is to be recomputed. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `inertia_tensor` | SMatrix{3,3,Float64} | kg*m^2 | — | Composite inertia tensor of the assembly about its centre of mass, also cached on the model. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[vehicle.mass_properties_get_com|get_COM]] · `callees` → `callers` · call · `src/vehicle/structure/mass_properties.jl:27-27`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/vehicle/structure/mass_properties.jl:41-41`
- `callees` → [[vehicle.mass_properties_update_inertia_tensor|update_inertia_tensor]] · `callers` · call · `src/vehicle/structure/mass_properties.jl:36-36`
- `callees` → [[vehx.structure_assembly_graph_traverse_bodies|traverse_bodies]] · `callers` · call · `src/vehicle/structure/mass_properties.jl:33-33`
<!-- vulcan:connections:end -->

## Limitations
The routine reads `model.prop_mass[root_index]` and writes `model.inertia_tensors`, both of which are vector-valued fields belonging to a multi-assembly model variant; against the current scalar `prop_mass` and single `inertia_tensor` fields of `SpacecraftModel` this path does not apply unchanged. Nothing recomputes the tensor automatically, so a caller that mutates a link mass and forgets this call will propagate with a stale tensor. Fuel slosh, flexible modes and reaction wheel rotor inertia are outside the model.

## Provenance
Mapped from `src/vehicle/structure/mass_properties.jl:32-46`; the pure summation it delegates to begins at line 53 of the same file, and the centre of mass helpers precede it.
