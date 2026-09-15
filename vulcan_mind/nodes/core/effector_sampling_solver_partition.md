---
id: core.effector_sampling_solver_partition
label: solver_partition
kind: function
source:
  file: src/core/types/effector_sampling.jl
  symbol: solver_partition
  lines:
  - 184
  - 184
inputs:
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
  type: Symbol
  units: n/a
  description: Return value of `solver_partition`. Returns `:explicit`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# solver_partition

## Purpose
Optional declaration hook telling the `split_imex` solver mode whether an effector's contribution should be integrated on the stiff implicit side or the non-stiff explicit side, so stiff effectors like high-rate attitude control or dense-atmosphere drag can be handled implicitly without changing the other effectors.

## Design & Implementation
`@inline solver_partition(::Any) = :explicit` is the only method defined here. Effector packages return `:implicit` from an overload to move their force onto the implicit operator. The engine evaluates the symbol once when partitioning the effector list at problem construction time.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `solver_partition`. Returns `:explicit`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.effector_sampling_wrench_caching_bang|wrench_caching!]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:175-175`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/effector_sampling.jl`
- [[simulation.dynamics_rhs__solver_partition_validated|_solver_partition_validated]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:23-23`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:295-295`

**Downstream**

- `callees` → [[core.effector_sampling_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/core/types/effector_sampling.jl:187-187`
- `callees` → [[dynamics.perturbations_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/core/types/effector_sampling.jl:187-187`
- `callees` → [[environment.gravity_models_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/core/types/effector_sampling.jl:187-187`
<!-- vulcan:connections:end -->

## Limitations
The return value is an unchecked `Symbol`; anything other than `:implicit` or `:explicit` is not rejected here and its treatment depends on the engine's partitioning code. Partition is static per model type, so an effector cannot switch sides as the trajectory enters or leaves the atmosphere. Only translational and rotational wrench contributions are partitioned; mass-flow or thermal states are outside this hook.

## Provenance
Mapped from `src/core/types/effector_sampling.jl` line 184.
