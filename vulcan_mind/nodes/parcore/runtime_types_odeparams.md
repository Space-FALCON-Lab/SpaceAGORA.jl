---
id: parcore.runtime_types_odeparams
label: ODEParams
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: ODEParams
  lines:
  - 848
  - 869
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Configuration, cache and workspace records declared earlier in the
    ConfigTypes module.
- id: config
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Typed simulation configuration stored in the args field and used as
    the struct's type parameter.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: params
  type: ODEParams
  units: n/a
  description: Integrator parameter object carrying the configuration, shared buffers,
    per-satellite activity flags, orbit counters and the save cache through every
    derivative and callback evaluation.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- parcore
origin: agent
---

# ODEParams

## Purpose
`ODEParams` is the parameter object passed to the differential-equation solver. It is the single handle through which the right-hand side, the callbacks and the saving logic reach the configuration, the preallocated buffers and the mutable per-satellite bookkeeping of a run.

## Model & Assumptions
The struct is parameterised on `A <: SimulationConfiguration`, so the configuration's concrete model types propagate into the parameter object and the derivative function specialises on them. `n_sats` fixes the width of the per-satellite vectors: `is_active` marks which satellites are still being propagated, and `orbit_counter` accumulates completed orbits. `shared_buffers` and `save_cache` are mutable and are reused across steps, which is what keeps the inner loop allocation-free.

## Design & Implementation
The record sits at the end of a long `ConfigTypes` module that also declares the initial condition, aerodynamics, engine, model, configuration, solution, ephemeris cache, scratch workspace and right-hand-side planning types. A comment above the declaration explains that the struct's default positional constructor and the keyword constructor defined immediately below it coexist deliberately, because the keyword form takes zero positional arguments while the default takes six, so the signatures cannot collide. The keyword constructor supplies `SimulationConfiguration()`, a `SharedBuffers` sized by `n_sats`, an all-true activity vector and a fresh save cache as defaults.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Configuration, cache and workspace records declared earlier in the ConfigTypes module. |
| in | `config` | SimulationConfiguration | n/a | yes | Typed simulation configuration stored in the args field and used as the struct's type parameter. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `params` | ODEParams | n/a | — | Integrator parameter object carrying the configuration, shared buffers, per-satellite activity flags, orbit counters and the save cache through every derivative and callback evaluation. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:184-184`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The mutable fields make the parameter object unsafe to share across concurrently propagated satellites unless the outer route gives each its own instance; the process route sidesteps this by construction, the thread route relies on per-task copies. Deactivating a satellite through `is_active` does not shrink the state vector, so the solver keeps integrating a frozen block. Buffer sizes are fixed at construction, so a configuration change that alters the state width requires a rebuilt parameter object rather than an in-place edit.

## Provenance
Mapped from `src/core/types/runtime_types.jl:848-869`.
