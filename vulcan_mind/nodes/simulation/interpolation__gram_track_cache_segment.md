---
id: simulation.interpolation__gram_track_cache_segment
label: _gram_track_cache_segment
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/interpolation.jl
  symbol: _gram_track_cache_segment
  lines:
  - 49
  - 49
inputs:
- id: cache
  type: GramTrackCache
  units: n/a
  required: true
  description: Positional argument `cache`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: ignore_time_window
  type: Bool
  units: n/a
  required: false
  description: Positional argument `ignore_time_window` (default `_gram_track_cache_ignore_time_window()`).
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_gram_track_cache_segment`.
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

# _gram_track_cache_segment

## Purpose
Locates the cache segment bracketing a query time and returns its index together with the normalised position within it, the lookup step every cached atmosphere evaluation starts from.

## Design & Implementation
Returns `Union{Nothing, Tuple{Int, Float64}}`. It bails out with `nothing` if `!cache.valid` or fewer than two samples exist, and, unless `ignore_time_window` is set, if `t` falls outside `[cache.t0, cache.t1]`; when that flag is set the query is clamped into range instead. The search starts from `cache.index_hint` clamped to `1:n-1`. If `tq < times[idx]`, time moved backwards — a rejected integrator step — and it falls back to `searchsortedlast`. If `tq > times[idx+1]`, it first tries advancing the hint by one, the O(1) case for a monotonically advancing integrator, and only then binary searches. The fraction is `x = (tq - t_lo) / (t_hi - t_lo)`, guarded so that a zero-width segment yields `0.0`, and `cache.index_hint` is mutated to the found index before returning.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | GramTrackCache | n/a | yes | Positional argument `cache`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `ignore_time_window` | Bool | n/a | no | Positional argument `ignore_time_window` (default `_gram_track_cache_ignore_time_window()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_gram_track_cache_segment`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:127-127`
- [[simulation_a.interpolation_gram_track_cache_ready|_gram_track_cache_ready]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:102-102`

**Downstream**

- `callees` → [[simulation.config__gram_track_cache_ignore_time_window|_gram_track_cache_ignore_time_window]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:52-52`
<!-- vulcan:connections:end -->

## Limitations
Mutating `cache.index_hint` makes this function stateful and not thread-safe: two tasks evaluating the same cache concurrently will fight over the hint, degrading it to a binary search each time or, worse, racing on the write. The `ignore_time_window` default is read from a global, `_gram_track_cache_ignore_time_window()`, so behaviour depends on process-wide state not visible at the call site. Clamping under that flag silently extrapolates the endpoints rather than reporting that the query was out of range, and a non-monotonic `cache.times` array breaks the `searchsortedlast` fallback without any check.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/interpolation.jl` line 49.
