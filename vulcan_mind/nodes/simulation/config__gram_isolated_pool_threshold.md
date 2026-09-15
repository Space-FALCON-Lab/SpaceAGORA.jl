---
id: simulation.config__gram_isolated_pool_threshold
label: _gram_isolated_pool_threshold
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _gram_isolated_pool_threshold
  lines:
  - 80
  - 80
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
  description: Return value of `_gram_isolated_pool_threshold`.
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

# _gram_isolated_pool_threshold

## Purpose
Supplies the item count at which automatic mode routes GRAM evaluations to the isolated worker pool.

## Design & Implementation
Returns `ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD", 4)`, defaulting to four items. `_gram_isolated_pool_enabled` combines it with a `Threads.nthreads() > 1` check so a single-threaded session never pays pool overhead.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_gram_isolated_pool_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__gram_isolated_pool_enabled|_gram_isolated_pool_enabled]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:95-95`
- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:195-195`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:81-81`
<!-- vulcan:connections:end -->

## Limitations
The threshold counts items without weighting their cost, so four cheap surrogate queries and four native GRAM calls are treated identically. Because the enclosing mode defaults to off, this threshold is inert unless the pool mode is explicitly switched to automatic.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 80.
