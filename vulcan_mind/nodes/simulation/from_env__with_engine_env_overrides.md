---
id: simulation.from_env__with_engine_env_overrides
label: _with_engine_env_overrides
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _with_engine_env_overrides
  lines:
  - 273
  - 273
inputs:
- id: config
  type: SimulationEngineConfig
  units: n/a
  required: true
  description: Positional argument `config`.
- id: f
  type: Function
  units: n/a
  required: true
  description: Positional argument `f`.
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
  description: Return value of `_with_engine_env_overrides`. Returns `f()`.
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

# _with_engine_env_overrides

## Purpose
Runs a function inside a scope where a `SimulationEngineConfig` is the authoritative source of `SPACEAGORA_*` settings, by installing its override dictionary both in the engine's global refs and in the process `ENV`, and restoring everything afterwards.

## Design & Implementation
`_with_engine_env_overrides(config::SimulationEngineConfig, f::Function)` (plus a `do`-block friendly method with arguments swapped) computes `overrides = _engine_env_overrides(config)`, saves the previous `_engine_active_config_ref[]` and `_engine_active_overrides_ref[]`, and sets them to `config` and `overrides`. If `overrides` is empty it runs `f()` in a `try`/`finally` that restores the refs. Otherwise it records each key's prior `ENV` value (or `nothing`), assigns `ENV[k] = v` for every override, runs `f()`, and in `finally` deletes keys that were previously absent, restores the others, and resets both refs. The return value of `f()` is propagated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | SimulationEngineConfig | n/a | yes | Positional argument `config`. |
| in | `f` | Function | n/a | yes | Positional argument `f`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_with_engine_env_overrides`. Returns `f()`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.public_api_prewarm_nbody_ephemeris_cache|prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/public_api.jl:36-36`
- [[simulation.run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/public_api.jl:20-20`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:280-280`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:280-280`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:280-280`
- `callees` → [[simulation.from_env__engine_env_overrides|_engine_env_overrides]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:274-274`
<!-- vulcan:connections:end -->

## Limitations
`ENV` and the two `Ref`s are process-global, so nested or concurrent scopes on different tasks corrupt each other; the restore logic is correct only for strictly nested single-threaded use. The empty-override branch is effectively unreachable since `_engine_env_overrides` always emits at least 23 keys. Restoring `ENV` key by key is not atomic, and an exception during restore would leave a partially modified environment.

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 273.
