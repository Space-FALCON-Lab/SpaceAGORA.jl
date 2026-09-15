---
id: simulation.config__control_callback_parallel_mode
label: _control_callback_parallel_mode
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _control_callback_parallel_mode
  lines:
  - 98
  - 98
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
  description: Return value of `_control_callback_parallel_mode`.
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

# _control_callback_parallel_mode

## Purpose
Resolves the requested threading mode for the per-satellite control callback.

## Design & Implementation
Forwards to `ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_CONTROL_CALLBACK_PARALLEL")` and returns the resulting `Symbol`. It is snapshotted as `CallbackEnvConfig.control_parallel_mode` and supplied as the `mode` keyword of `ParallelPolicy.thread_policy_decision` from `_control_callback_thread_decision`, whose policy source tag is `:control_callback`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_control_callback_parallel_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:197-197`

**Downstream**

- `callees` → [[parallel.env_config_parse_parallel_mode_env|parse_parallel_mode_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:99-99`
<!-- vulcan:connections:end -->

## Limitations
The mode is advisory: `control_model_threadsafe` defaults to `false` for every type except `BaseThrusterModel`, so requesting `:on` with a custom controller produces a serial run unless the assume-threadsafe override is also set. Unsupported spellings throw from the parser.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 98.
