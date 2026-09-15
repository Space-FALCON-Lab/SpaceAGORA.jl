---
id: simulation.registry__gram_runtime_stats_reset_bang
label: _gram_runtime_stats_reset!
kind: function
source:
  file: src/simulation/callbacks/registry.jl
  symbol: _gram_runtime_stats_reset!
  lines:
  - 69
  - 69
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
  type: Nothing
  units: n/a
  description: Return value of `_gram_runtime_stats_reset!`. Returns `nothing`.
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

# _gram_runtime_stats_reset!

## Purpose
Clears the process-wide GRAM runtime profiling counters by replacing the global `GramRuntimeStats` record with a freshly constructed default instance.

## Design & Implementation
Acquires `_gram_runtime_stats_lock` (a `ReentrantLock`) and assigns `_gram_runtime_stats[] = GramRuntimeStats()`, which resets all `Int64` counters (`density_calls`, `cache_hits`, `refresh_calls`, ...) to 0 and the `Float64` accumulators (`refresh_elapsed_s`, the `*_err_abs_max_*` and `*_err_abs_sum_*` fields) to 0.0. Returns `nothing`. It is intended to be called at the start of a profiled run so that `_gram_runtime_stats_snapshot` reports statistics for that run only.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_gram_runtime_stats_reset!`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/registry.jl`

**Downstream**

- `callees` → [[simulation_a.registry_gramruntimestats|GramRuntimeStats]] · `callers` · call · `src/simulation/callbacks/registry.jl:71-71`
<!-- vulcan:connections:end -->

## Limitations
The stats object is a single module-level `Ref`, so a reset from one task discards counters accumulated by any other concurrently running simulation in the same process. There is no history: the previous record is dropped rather than archived, and only the lock prevents a reset while another thread is inside an update closure.

## Provenance
Mapped from `src/simulation/callbacks/registry.jl` line 69.
