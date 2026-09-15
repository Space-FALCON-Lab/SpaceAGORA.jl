---
id: simulation.registry__gram_runtime_stats_snapshot
label: _gram_runtime_stats_snapshot
kind: function
source:
  file: src/simulation/callbacks/registry.jl
  symbol: _gram_runtime_stats_snapshot
  lines:
  - 83
  - 83
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
  type: Any
  units: n/a
  description: Return value of `_gram_runtime_stats_snapshot`. Returns `(`.
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

# _gram_runtime_stats_snapshot

## Purpose
Produces an immutable copy of the current GRAM profiling counters as a `NamedTuple`, so reporting code can read consistent values without holding the lock or aliasing the mutable global.

## Design & Implementation
Takes no arguments. Inside `lock(_gram_runtime_stats_lock) do ... end` it reads `_gram_runtime_stats[]` and copies all nineteen fields of `GramRuntimeStats` into a `NamedTuple` with identical names: call counters (`density_calls`, `cache_enabled_calls`, `cache_hits`, `cache_misses`, `miss_time_window`, `miss_state_tolerance`, `direct_calls`), refresh bookkeeping (`refresh_calls`, `refresh_points_total`, `refresh_points_max`, `refresh_failures`, `refresh_elapsed_s` in seconds) and state-error tracking (`state_error_samples`, max and summed absolute altitude error in m and latitude/longitude error in degrees). Because the copy is taken under the lock, the returned fields are mutually consistent.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gram_runtime_stats_snapshot`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/registry.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The field list is duplicated by hand; adding a field to `GramRuntimeStats` without updating this function silently omits it from the snapshot. Mean errors must be derived by the caller as `*_sum / state_error_samples`, and a snapshot taken when `state_error_samples == 0` will give division by zero if the caller does not guard it.

## Provenance
Mapped from `src/simulation/callbacks/registry.jl` line 83.
