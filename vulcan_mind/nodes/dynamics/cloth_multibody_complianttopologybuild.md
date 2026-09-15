---
id: dynamics.cloth_multibody_complianttopologybuild
label: CompliantTopologyBuild
kind: struct
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: CompliantTopologyBuild
  lines:
  - 80
  - 80
inputs:
- id: model
  type: CompliantMultibodyModel
  units: n/a
  required: true
  description: Field `model`.
- id: initial_state
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `initial_state`.
- id: node_names
  type: Vector{Symbol}
  units: n/a
  required: true
  description: Field `node_names`.
- id: joint_names
  type: Vector{Symbol}
  units: n/a
  required: true
  description: Field `joint_names`.
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
  description: Constructed `CompliantTopologyBuild`.
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

# CompliantTopologyBuild

## Purpose
The product of a topology build: the model, its initial state vector, and name-to-index tables for bodies and joints.

## Design & Implementation
Immutable with `model`, `initial_state`, `node_names` and `joint_names`. The name vectors let callers find indices by symbol without threading them through separately.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | CompliantMultibodyModel | n/a | yes | Field `model`. |
| in | `initial_state` | Vector{Float64} | n/a | yes | Field `initial_state`. |
| in | `node_names` | Vector{Symbol} | n/a | yes | Field `node_names`. |
| in | `joint_names` | Vector{Symbol} | n/a | yes | Field `joint_names`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CompliantTopologyBuild | n/a | — | Constructed `CompliantTopologyBuild`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_build_compliant_topology|build_compliant_topology]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:322-322`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Lookup by name is a linear `findfirst` over the vectors, not a dictionary.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 80.
