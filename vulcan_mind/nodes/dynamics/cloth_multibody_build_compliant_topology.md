---
id: dynamics.cloth_multibody_build_compliant_topology
label: build_compliant_topology
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: build_compliant_topology
  lines:
  - 281
  - 281
inputs:
- id: nodes
  type: AbstractVector{CompliantTopologyNode}
  units: n/a
  required: true
  description: Positional argument `nodes`.
- id: edges
  type: AbstractVector{CompliantTopologyEdge}
  units: n/a
  required: true
  description: Positional argument `edges`.
- id: base_position
  type: Any
  units: n/a
  required: false
  description: Keyword argument `base_position` (default `(0.0, 0.0, 0.0)`).
- id: base_quaternion
  type: Any
  units: n/a
  required: false
  description: Keyword argument `base_quaternion` (default `_Q_IDENTITY`).
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
  type: CompliantTopologyBuild
  units: n/a
  description: Return value of `build_compliant_topology`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# build_compliant_topology

## Purpose
Turns topology nodes and edges into a runnable model plus initial state, validating connectivity indices and resolving joint rest orientations.

## Design & Implementation
Requires at least one node. Bodies are copied from nodes. For each edge it checks the child index is in range, the parent is in zero to n, and parent differs from child, raising `ArgumentError` naming the edge. The rest quaternion is taken from the edge or computed from the parent's (or base's) and child's initial quaternions. The initial state is packed from node positions, quaternions and rates. Returns a `CompliantTopologyBuild`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `nodes` | AbstractVector{CompliantTopologyNode} | n/a | yes | Positional argument `nodes`. |
| in | `edges` | AbstractVector{CompliantTopologyEdge} | n/a | yes | Positional argument `edges`. |
| in | `base_position` | Any | n/a | no | Keyword argument `base_position` (default `(0.0, 0.0, 0.0)`). |
| in | `base_quaternion` | Any | n/a | no | Keyword argument `base_quaternion` (default `_Q_IDENTITY`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CompliantTopologyBuild | n/a | — | Return value of `build_compliant_topology`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_build_rectangular_compliant_grid|build_rectangular_compliant_grid]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:432-432`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__rest_child_parent_quat|_rest_child_parent_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:301-301`
- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:289-289`
- `callees` → [[dynamics.cloth_multibody_compliant_state_vector|compliant_state_vector]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:316-316`
- `callees` → [[dynamics.cloth_multibody_compliantbody|CompliantBody]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:291-291`
- `callees` → [[dynamics.cloth_multibody_compliantjoint|CompliantJoint]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:303-303`
- `callees` → [[dynamics.cloth_multibody_compliantmultibodymodel|CompliantMultibodyModel]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:323-323`
- `callees` → [[dynamics.cloth_multibody_complianttopologybuild|CompliantTopologyBuild]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:322-322`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__rest_child_parent_quat|_rest_child_parent_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:301-301`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:289-289`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:303-303`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:289-289`
<!-- vulcan:connections:end -->

## Limitations
Connectivity is validated per edge only; nothing checks that the graph is connected or acyclic, so an isolated body or a loop is accepted and simply behaves as its springs dictate.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 281.
