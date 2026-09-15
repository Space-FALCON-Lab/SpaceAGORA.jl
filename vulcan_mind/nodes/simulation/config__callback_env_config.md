---
id: simulation.config__callback_env_config
label: _callback_env_config
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _callback_env_config
  lines:
  - 209
  - 209
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  type: CallbackEnvConfig
  units: n/a
  description: Return value of `_callback_env_config`.
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

# _callback_env_config

## Purpose
Run-scoped accessor that returns the `CallbackEnvConfig` snapshot for the current simulation, so hot callback paths read plain struct fields instead of re-parsing environment variables.

## Design & Implementation
Takes the ODE parameter object `p`. When `p` is not `nothing` and `hasproperty(p, :shared_buffers)`, it fetches `shared_buffers`, checks `hasproperty(sb, :callback_env_config)`, and dereferences the `Ref` with `sb.callback_env_config[]`. A non-`nothing` snapshot is returned directly. In every other case — `p === nothing`, a hand-constructed parameter object lacking the field, or an unset `Ref` — it falls back to `_snapshot_callback_env_config()`, which re-parses the whole environment.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CallbackEnvConfig | n/a | — | Return value of `_callback_env_config`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__control_callback_thread_decision|_control_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:302-302`
- [[simulation.config__density_callback_thread_decision|_density_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:257-257`
- [[simulation.config__gram_track_trajectory_supported|_gram_track_trajectory_supported]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:175-175`
- [[simulation.config__thermal_callback_thread_decision|_thermal_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:333-333`
- [[simulation.effector_sampling__buffered_atmosphere_valid|_buffered_atmosphere_valid]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:135-135`
- [[simulation.model_selection__gram_isolated_pool_batch_model_for_callback|_gram_isolated_pool_batch_model_for_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:71-71`
- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:94-94`
- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:254-254`
- [[simulation.runtime_update_density_sat_bang|update_density_sat!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:223-223`
- [[simulation_a.model_selection_gram_isolated_pool_batch_eval__gram_isolated_pool_batch_eval_bang|_gram_isolated_pool_batch_eval!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:185-185`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:223-223`

**Downstream**

- `callees` → [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:217-217`
<!-- vulcan:connections:end -->

## Limitations
The fallback is silent and expensive: a parameter object missing the field re-reads and re-validates roughly two dozen environment variables on every callback invocation, which can dominate step time in unit tests and `withenv` probes. The `hasproperty` probes are duck-typed, so a field of the wrong type is only detected when the returned value is used.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 209.
