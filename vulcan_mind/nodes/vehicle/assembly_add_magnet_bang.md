---
id: vehicle.assembly_add_magnet_bang
label: add_magnet!
kind: function
source:
  file: src/vehicle/spacecraft/assembly.jl
  symbol: add_magnet!
  lines:
  - 66
  - 66
inputs:
- id: link
  type: Link
  units: n/a
  required: true
  description: Positional argument `link`.
- id: magnet
  type: Magnet
  units: n/a
  required: true
  description: Positional argument `magnet`.
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
  description: Return value of `add_magnet!`; mutates `link` in place.
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

# add_magnet!

## Purpose

`add_magnet!` attaches magnetic dipole sources to a `Link`, populating the `magnets` collection that the magnetic torque model integrates against the ambient field. It covers both permanent magnets and magnetorquer-style dipoles represented by the `Magnet` component type.

## Design & Implementation

The symbol has a scalar and a vector method, mirroring `add_facet!`. `add_magnet!(link::Link, magnet::Magnet)` performs `push!(link.magnets, magnet)`; `add_magnet!(link::Link, magnets::Vector{Magnet})` splats with `push!(link.magnets, magnets...)` so an array of dipoles is appended in order. Both mutate the link's `magnets` vector in place and return that vector. No aggregate dipole is precomputed and no counter on the parent `SpacecraftModel` is incremented.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `link` | Link | n/a | yes | Positional argument `link`. |
| in | `magnet` | Magnet | n/a | yes | Positional argument `magnet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `add_magnet!`; mutates `link` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/assembly.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/vehicle/spacecraft/assembly.jl:67-67`
<!-- vulcan:connections:end -->

## Limitations

Dipole moment direction and magnitude are accepted verbatim, with no unit-vector normalisation and no check that the mounting location lies on the link. Repeated calls with the same `Magnet` silently sum its contribution twice. Since the total dipole is recomputed from the vector at every torque evaluation, very large magnet counts cost time in the inner dynamics loop.

## Provenance
Mapped from `src/vehicle/spacecraft/assembly.jl` line 66.
