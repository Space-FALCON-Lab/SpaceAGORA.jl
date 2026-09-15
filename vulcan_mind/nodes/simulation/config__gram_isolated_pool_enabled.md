---
id: simulation.config__gram_isolated_pool_enabled
label: _gram_isolated_pool_enabled
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _gram_isolated_pool_enabled
  lines:
  - 88
  - 88
inputs:
- id: num_items
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_items`.
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
  type: Bool
  units: n/a
  description: Return value of `_gram_isolated_pool_enabled`.
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

# _gram_isolated_pool_enabled

## Purpose
Decides whether a batch of GRAM evaluations of a given size should be routed through the isolated worker pool.

## Design & Implementation
Two methods mirror the batch pair. The single-argument form reads the mode live; the `(env::CallbackEnvConfig, num_items)` form reads `env.gram_isolated_pool_mode` and `env.gram_isolated_pool_threshold`. Mode `:off` returns `false`; `:on` returns `num_items > 0`; automatic mode additionally requires `Threads.nthreads() > 1` before comparing `num_items` against the threshold.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_gram_isolated_pool_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__policy_env_config|_policy_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:242-242`
- [[simulation.model_selection__gram_isolated_pool_batch_model_for_callback|_gram_isolated_pool_batch_model_for_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:51-51`
- [[simulation_a.model_selection_gram_isolated_pool_batch_eval__gram_isolated_pool_batch_eval_bang|_gram_isolated_pool_batch_eval!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:186-186`

**Downstream**

- `callees` → [[simulation.config__gram_isolated_pool_mode|_gram_isolated_pool_mode]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:89-89`
- `callees` → [[simulation.config__gram_isolated_pool_threshold|_gram_isolated_pool_threshold]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:95-95`
<!-- vulcan:connections:end -->

## Limitations
The `Threads.nthreads() > 1` guard applies only in automatic mode, so an explicit `:on` enables the pool even in a single-threaded session, where it adds dispatch overhead for no concurrency. The function does not check that the pool has been constructed, leaving that error to surface later at dispatch time.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 88.
