---
id: dynamics.aerodynamic_wrench_models_solver_partition
label: solver_partition
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: solver_partition
  lines:
  - 227
  - 227
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
  description: Return value of `solver_partition`. Returns `:implicit`.
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

# solver_partition

## Purpose
Assigns all aerodynamic effectors to the `:implicit` partition so split solvers treat drag as the potentially stiff component.

## Design & Implementation
Three `@inline` methods returning the symbol `:implicit` for `AerodynamicCoefficientConstant`, `AerodynamicCoefficientfM`, and `AerodynamicCoefficientNoBallisticFlight`. The split-IMEX and multirate solver modes use this to route the effector into `f1` versus `f2`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `solver_partition`. Returns `:implicit`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.effector_sampling_wrench_caching_bang|wrench_caching!]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:175-175`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- [[simulation.dynamics_rhs__solver_partition_validated|_solver_partition_validated]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:23-23`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:295-295`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The partition is fixed regardless of altitude; in near-vacuum the aerodynamic term is negligible yet still forces the implicit stage. Effectors that ignore the partition (plain `Tsit5`) are unaffected.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 227.
