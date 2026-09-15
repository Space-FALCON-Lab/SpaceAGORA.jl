---
id: simulation.config__thermal_callback_thread_threshold
label: _thermal_callback_thread_threshold
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _thermal_callback_thread_threshold
  lines:
  - 117
  - 117
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
  type: Int
  units: n/a
  description: Return value of `_thermal_callback_thread_threshold`.
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

# _thermal_callback_thread_threshold

## Purpose
Supplies the satellite count at which the thermal callback may be threaded, inheriting the density threshold when unset.

## Design & Implementation
If `SPACEAGORA_THERMAL_CALLBACK_THREAD_THRESHOLD` is present in `ENV` it is parsed by `ParallelPolicy.parse_thread_threshold_env` with a default of 8; otherwise the function returns `_density_callback_thread_threshold()`. The resolved integer becomes `CallbackEnvConfig.thermal_thread_threshold`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_thermal_callback_thread_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:202-202`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:119-119`
- `callees` → [[simulation.config__density_callback_thread_threshold|_density_callback_thread_threshold]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:121-121`
<!-- vulcan:connections:end -->

## Limitations
The inherited value tracks the density threshold, which is tuned for atmosphere-model cost rather than thermal-node integration cost, so the crossover is unlikely to be optimal for either when they differ. Presence-based inheritance means an empty export disables the fallback and throws.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 117.
