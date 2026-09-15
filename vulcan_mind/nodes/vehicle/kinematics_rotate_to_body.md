---
id: vehicle.kinematics_rotate_to_body
label: rotate_to_body
kind: function
source:
  file: src/vehicle/kinematics/kinematics.jl
  symbol: rotate_to_body
  lines:
  - 22
  - 22
inputs:
- id: body
  type: Link
  units: n/a
  required: true
  description: Positional argument `body`.
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
  type: I
  units: n/a
  description: Return value of `rotate_to_body`. Returns `I(3)` or `rot(body.q)'`.
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

# rotate_to_body

## Purpose

`rotate_to_body(body::Link)` returns the 3x3 rotation matrix that maps a vector expressed in the given link's own frame into the root body frame of the assembly it belongs to.

## Design & Implementation

The function branches on the `body.root` flag. For a root link the link frame and the body frame are by definition the same, so it returns the identity `I(3)`. For a child link, whose quaternion `body.q` is stored relative to the root body frame, it returns the transpose `rot(body.q)'` of the direction cosine matrix built by `rot`, which is the inverse rotation because a DCM is orthonormal. The result is used alongside `rotate_to_inertial` when composing inertia tensors and offsets across an articulated assembly.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `body` | Link | n/a | yes | Positional argument `body`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | I | n/a | — | Return value of `rotate_to_body`. Returns `I(3)` or `rot(body.q)'`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[vehx.actuators_thruster_hooks_update_thrusters_bang|update_thrusters!]] · `callees` → `callers` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl:24-24`

**Downstream**

- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/vehicle/kinematics/kinematics.jl:29-29`
<!-- vulcan:connections:end -->

## Limitations

The identity branch returns `I(3)`, a `UniformScaling`-backed `Diagonal`-like object rather than the `SMatrix{3,3,Float64}` returned by the child branch, so the method is not type-stable and callers that store the result in a concrete matrix field must convert. Correctness depends entirely on `body.q` being unit norm and on the child-relative-to-root storage convention being honoured by whoever wrote the link.

## Provenance
Mapped from `src/vehicle/kinematics/kinematics.jl` line 22.
