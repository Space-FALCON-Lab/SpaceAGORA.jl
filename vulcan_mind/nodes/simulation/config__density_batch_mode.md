---
id: simulation.config__density_batch_mode
label: _density_batch_mode
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _density_batch_mode
  lines:
  - 58
  - 58
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
  description: Return value of `_density_batch_mode`.
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

# _density_batch_mode

## Purpose
Resolves whether density evaluations for multiple satellites should be batched into a single vectorised model call.

## Design & Implementation
Delegates to `ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_DENSITY_BATCH_PARALLEL")`, returning `:off`, `:on`, or the automatic mode that defers to `_density_batch_threshold`. The result is stored as `CallbackEnvConfig.density_batch_mode` and read by both methods of `_density_batch_enabled`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_density_batch_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__density_batch_enabled|_density_batch_enabled]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:67-67`
- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:192-192`

**Downstream**

- `callees` → [[parallel.env_config_parse_parallel_mode_env|parse_parallel_mode_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:59-59`
<!-- vulcan:connections:end -->

## Limitations
Batching is orthogonal to threading but shares the same parser and vocabulary, which makes the two families of environment variables easy to confuse. The mode says nothing about whether the selected density model actually implements a batched entry point.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 58.
