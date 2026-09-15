---
id: simulation.dynamics_rhs__solver_partition_validated
label: _solver_partition_validated
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _solver_partition_validated
  lines:
  - 22
  - 22
inputs:
- id: effector
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector`.
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
  description: Return value of `_solver_partition_validated`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _solver_partition_validated

## Purpose
Reads an effector's declared solver partition and rejects anything other than `:implicit` or `:explicit`, so a mis-declared effector fails loudly instead of being silently dropped from both partitions.

## Design & Implementation
Calls `solver_partition(effector)` and returns the symbol if valid, else raises `ArgumentError` naming the effector type and the bad value. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_solver_partition_validated`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`

**Downstream**

- `callees` → [[core.effector_sampling_solver_partition|solver_partition]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:23-23`
- `callees` → [[dynamics.aerodynamic_wrench_models_solver_partition|solver_partition]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:23-23`
<!-- vulcan:connections:end -->

## Limitations
Called on every partitioned evaluation rather than once at setup.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 22.
