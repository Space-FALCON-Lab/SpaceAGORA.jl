---
id: simulation.config_gramtrackcache
label: GramTrackCache
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/config.jl
  symbol: GramTrackCache
  lines:
  - 1
  - 1
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
  type: GramTrackCache
  units: n/a
  description: Return value of `GramTrackCache`. Returns `GramTrackCache(     false,     0.0,     0.0,     1,     Float64[],     Float64[]`.
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

# GramTrackCache

## Purpose
This zero-argument `GramTrackCache()` constructor produces an empty, explicitly-invalid GRAM ground-track cache ready to be filled in by the track-cache refresh callback. It is the allocation point for the per-run cache object that lets the atmosphere model interpolate density, temperature and wind along a precomputed ground track instead of calling the GRAM backend at every RHS evaluation.

## Design & Implementation
Marked `@inline`, it calls the eleven-field inner constructor of the mutable `GramTrackCache` defined in `src/core/types/runtime_types.jl` with `valid = false`, `t0 = t1 = 0.0`, `index_hint = 1`, six empty `Float64[]` vectors for `times`, `alts`, `lats`, `lons`, `rhos` and `Ts`, and an empty `SVector{3, Float64}[]` for `winds`. Because the type is mutable, the refresh path resizes these vectors in place rather than reallocating a new cache per segment; `index_hint` seeds the monotonic search used by the interpolator so successive queries at increasing time cost O(1).

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | GramTrackCache | n/a | — | Return value of `GramTrackCache`. Returns `GramTrackCache(     false,     0.0,     0.0,     1,     Float64[],     Float64[]`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/gram_track_cache/config.jl`
- [[simulation.model_selection__gram_density_cache_for_sat_bang|_gram_density_cache_for_sat!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:122-122`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The `valid = false` flag is the only guard against reading the empty vectors — any consumer that indexes `times` or `rhos` without checking it hits a `BoundsError`. The constructor imposes no invariant tying the six sample vectors to a common length, so a partially-completed refresh that is interrupted can leave the cache internally ragged while still marked valid. The vectors start at zero capacity, so the first refresh incurs growth reallocations for every field.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/config.jl` line 1.
