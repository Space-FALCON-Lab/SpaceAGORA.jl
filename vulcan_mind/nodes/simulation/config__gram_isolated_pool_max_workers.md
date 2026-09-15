---
id: simulation.config__gram_isolated_pool_max_workers
label: _gram_isolated_pool_max_workers
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _gram_isolated_pool_max_workers
  lines:
  - 84
  - 84
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
  description: Return value of `_gram_isolated_pool_max_workers`.
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

# _gram_isolated_pool_max_workers

## Purpose
Caps the number of workers the isolated GRAM pool may occupy.

## Design & Implementation
Returns `ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS", max(1, Threads.nthreads()))`, so the default tracks the thread count of the running Julia session and is never below one. Reusing the threshold parser means the same integer parsing and validation rules apply.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_gram_isolated_pool_max_workers`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:196-196`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:85-85`
<!-- vulcan:connections:end -->

## Limitations
The default is captured from `Threads.nthreads()` at call time, which reflects the session's thread count rather than any budget already claimed by an outer parallel construct, so combining the pool with a threaded ensemble can request more workers than exist. No upper bound is enforced against the machine's core count.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 84.
