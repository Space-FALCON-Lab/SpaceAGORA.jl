---
id: parallel.outer_route_state__route_stats_payload
label: _route_stats_payload
kind: function
source:
  file: src/parallel/routing/outer_route_state.jl
  symbol: _route_stats_payload
  lines:
  - 91
  - 91
inputs:
- id: stats
  type: OuterRouteStats
  units: n/a
  required: true
  description: Positional argument `stats`.
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
  type: Dict{String,
  units: n/a
  description: Return value of `_route_stats_payload`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# _route_stats_payload

## Purpose
Converts one `OuterRouteStats` into a plain `Dict{String, Any}` suitable for TOML serialisation by `save_outer_route_state`. It is the write-side counterpart of `_route_payload_stats`.

## Design & Implementation
Marked `@inline`. Produces five string-keyed entries: `"samples"`, `"successes"`, `"failures"` as `Int` clamped with `max(0, ...)`, and `"elapsed_sum_s"`, `"elapsed_sq_sum_s"` as `Float64` clamped with `max(0.0, ...)`. The clamping guarantees the persisted file never contains negative counts even if in-memory stats were corrupted. The function is pure and never throws for a well-typed `OuterRouteStats`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `stats` | OuterRouteStats | n/a | yes | Positional argument `stats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Dict{String, | n/a | — | Return value of `_route_stats_payload`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_state.jl`
- [[parallel.outer_route_state_save_outer_route_state|save_outer_route_state]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_state.jl:166-166`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/parallel/routing/outer_route_state.jl:96-96`
<!-- vulcan:connections:end -->

## Limitations
Clamping hides corruption rather than reporting it; a negative `samples` value becomes `0` and the row is then dropped by the loader because `samples <= 0`. `NaN` values in `elapsed_sq_sum_s` pass through `max(0.0, NaN)` as `NaN` under Julia's semantics and are written to TOML, relying on the loader's `isfinite` fallback. The schema version is not embedded here but in the caller.

## Provenance
Mapped from `src/parallel/routing/outer_route_state.jl` line 91.
