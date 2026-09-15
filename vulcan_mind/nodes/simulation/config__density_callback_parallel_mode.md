---
id: simulation.config__density_callback_parallel_mode
label: _density_callback_parallel_mode
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _density_callback_parallel_mode
  lines:
  - 46
  - 46
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
  description: Return value of `_density_callback_parallel_mode`.
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

# _density_callback_parallel_mode

## Purpose
Resolves the requested threading mode for the per-satellite density callback into a `Symbol` consumed by the shared parallel policy layer.

## Design & Implementation
Forwards to `ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_DENSITY_CALLBACK_PARALLEL")`, which yields the tri-state `:off`, `:on` or the automatic mode. The returned symbol is stored as `CallbackEnvConfig.density_parallel_mode` and later passed as the `mode` keyword of `ParallelPolicy.thread_policy_decision` inside `_density_callback_thread_decision`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_density_callback_parallel_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__thermal_callback_parallel_mode|_thermal_callback_parallel_mode]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:114-114`
- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:188-188`

**Downstream**

- `callees` → [[parallel.env_config_parse_parallel_mode_env|parse_parallel_mode_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:47-47`
<!-- vulcan:connections:end -->

## Limitations
The mode expresses intent only; it can be overridden downstream by the thread-safety gate on the density model and by the outer-parallel interlock, so requesting `:on` does not guarantee threads are used. Unsupported spellings throw inside the parallel policy parser.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 46.
