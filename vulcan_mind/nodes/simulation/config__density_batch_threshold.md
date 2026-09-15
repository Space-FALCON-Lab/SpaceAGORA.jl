---
id: simulation.config__density_batch_threshold
label: _density_batch_threshold
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _density_batch_threshold
  lines:
  - 62
  - 62
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
  description: Return value of `_density_batch_threshold`.
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

# _density_batch_threshold

## Purpose
Supplies the satellite count at which automatic mode switches density evaluation over to batched calls.

## Design & Implementation
Returns `ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_DENSITY_BATCH_THRESHOLD", 2)`, so batching engages from two satellites upward by default — much lower than the threading threshold of eight, reflecting that batching has no thread-spawn overhead to amortise.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_density_batch_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__density_batch_enabled|_density_batch_enabled]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:73-73`
- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:193-193`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:63-63`
<!-- vulcan:connections:end -->

## Limitations
The default of two means almost every constellation run batches, so a model whose batched path diverges numerically from its scalar path will change results as soon as a second satellite is added. The threshold is ignored under explicit `:on` or `:off`.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 62.
