---
id: dynamics.cloth_multibody_compliantmultibodytrajectory
label: CompliantMultibodyTrajectory
kind: struct
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: CompliantMultibodyTrajectory
  lines:
  - 49
  - 49
inputs:
- id: t_s
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `t_s`.
- id: states
  type: Vector{Vector{Float64}}
  units: n/a
  required: true
  description: Field `states`.
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
  type: CompliantMultibodyTrajectory
  units: n/a
  description: Constructed `CompliantMultibodyTrajectory`.
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

# CompliantMultibodyTrajectory

## Purpose
The time history produced by `simulate_compliant_multibody`: sample times and the full flat state at each.

## Design & Implementation
Immutable with `t_s::Vector{Float64}` and `states::Vector{Vector{Float64}}`, one thirteen-per-body vector per time. Keeping raw state vectors rather than unpacked poses lets callers use `compliant_state_parts` on any sample.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t_s` | Vector{Float64} | n/a | yes | Field `t_s`. |
| in | `states` | Vector{Vector{Float64}} | n/a | yes | Field `states`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CompliantMultibodyTrajectory | n/a | — | Constructed `CompliantMultibodyTrajectory`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_simulate_compliant_multibody|simulate_compliant_multibody]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:729-729`
- [[dynamics.cloth_robot_arm_dynamics_simulate_cloth_robot_arm_plan|simulate_cloth_robot_arm_plan]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:518-518`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Memory is `13 n_bodies` floats per sample with no thinning option, so a long run of a large grid at small `dt` is heavy.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 49.
