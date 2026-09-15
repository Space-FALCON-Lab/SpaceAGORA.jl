---
id: vehicle.assembly_add_joint_bang
label: add_joint!
kind: function
source:
  file: src/vehicle/spacecraft/assembly.jl
  symbol: add_joint!
  lines:
  - 48
  - 48
inputs:
- id: model
  type: SpacecraftModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: joint
  type: Joint
  units: n/a
  required: true
  description: Positional argument `joint`.
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
  description: Return value of `add_joint!`; mutates `model` in place.
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

# add_joint!

## Purpose

`add_joint!(model::SpacecraftModel, joint::Joint)` registers one articulation joint with the spacecraft model, making it part of the multibody topology that the dynamics later walks when propagating link states.

## Design & Implementation

The implementation is a single `push!(model.joints, joint)`, so the joint is appended to the model's `joints` vector in call order and that ordering becomes the joint indexing used downstream. The function mutates `model` and returns the result of `push!`, namely the `joints` vector itself, rather than a meaningful value. A commented-out call to `update_inertia_tensor!(model, joint.link1)` records that the mass-property update was intentionally moved out of assembly to break a dependency cycle with the analysis module.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | SpacecraftModel | n/a | yes | Positional argument `model`. |
| in | `joint` | Joint | n/a | yes | Positional argument `joint`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `add_joint!`; mutates `model` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/assembly.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/vehicle/spacecraft/assembly.jl:54-54`
<!-- vulcan:connections:end -->

## Limitations

No check is made that the joint's parent and child links have already been added with `add_body!`, that the joint does not duplicate an existing one, or that the resulting graph stays acyclic; a malformed topology surfaces only during propagation. Inertia properties are not refreshed, so callers must recompute them after assembly is complete.

## Provenance
Mapped from `src/vehicle/spacecraft/assembly.jl` line 48.
