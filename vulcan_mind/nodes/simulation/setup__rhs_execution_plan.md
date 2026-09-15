---
id: simulation.setup__rhs_execution_plan
label: _rhs_execution_plan
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_execution_plan
  lines:
  - 1023
  - 1023
inputs:
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
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
  type: SimulationModel.RhsExecutionPlan
  units: n/a
  description: Return value of `_rhs_execution_plan`.
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

# _rhs_execution_plan

## Purpose
Returns the execution plan for the current RHS call, serving it from the per-accepted-step cache when that cache is enabled.

## Design & Implementation
If shared buffers exist and `_rhs_plan_step_cache_enabled()`, returns the cached plan if present, otherwise computes it through `_rhs_execution_plan_uncached` and stores it. The planet-frame callback clears the cache once per accepted step. Without the cache it always recomputes. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.RhsExecutionPlan | n/a | — | Return value of `_rhs_execution_plan`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs_spacecraft_dynamics_explicit_remainder_bang|spacecraft_dynamics_explicit_remainder!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2085-2085`
- [[simulation.dynamics_rhs_spacecraft_dynamics_implicit_atmosphere_bang|spacecraft_dynamics_implicit_atmosphere!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1986-1986`
- [[simulation.dynamics_rhs_spacecraft_dynamics_slow_bang|spacecraft_dynamics_slow!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1838-1838`
- [[simx.engine_dynamics_rhs_spacecraft_dynamics_bang|spacecraft_dynamics!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1724-1724`

**Downstream**

- `callees` → [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callers` · call · `src/simulation/engine/setup.jl:1032-1032`
- `callees` → [[simulation.setup__rhs_plan_step_cache_enabled|_rhs_plan_step_cache_enabled]] · `callers` · call · `src/simulation/engine/setup.jl:1029-1029`
<!-- vulcan:connections:end -->

## Limitations
`_rhs_plan_step_cache_enabled()` is itself a live environment read on every RHS call, partly defeating the purpose; the snapshot mechanism was not extended to this flag.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1023.
