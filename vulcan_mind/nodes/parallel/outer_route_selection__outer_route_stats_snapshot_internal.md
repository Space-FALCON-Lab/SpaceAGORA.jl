---
id: parallel.outer_route_selection__outer_route_stats_snapshot_internal
label: _outer_route_stats_snapshot_internal
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _outer_route_stats_snapshot_internal
  lines:
  - 213
  - 213
inputs:
- id: state
  type: OuterRouteState
  units: n/a
  required: true
  description: Positional argument `state`.
- id: signature
  type: String
  units: n/a
  required: true
  description: Positional argument `signature`.
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
  type: Dict{Symbol,
  units: n/a
  description: Return value of `_outer_route_stats_snapshot_internal`.
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

# _outer_route_stats_snapshot_internal

## Purpose
Takes a thread-safe snapshot of per-route statistics for one signature from `OuterRouteState.history`, converting raw sums into `(samples, mean_s, success_rate, std_s)` tuples that the selection heuristics can read without holding the lock.

## Design & Implementation
Signature `_outer_route_stats_snapshot_internal(state::OuterRouteState, signature::String)`. Inside `lock(state.lock) do ... end` it fetches `entry = get(state.history, signature, nothing)`, returning an empty typed `Dict{Symbol, NamedTuple}` when absent. For each `(route, stats)` with `stats.samples > 0` it computes `success_rate = successes / max(1, samples)` and calls `_route_elapsed_stats(stats)`, storing the four-field tuple. Routes with zero samples are omitted, so an empty result means no usable history.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `state` | OuterRouteState | n/a | yes | Positional argument `state`. |
| in | `signature` | String | n/a | yes | Positional argument `signature`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Dict{Symbol, | n/a | — | Return value of `_outer_route_stats_snapshot_internal`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.outer_route_selection_outer_route_stats_snapshot|outer_route_stats_snapshot]] · `callees` → `callers` · feedback · `src/parallel/routing/outer_route_selection.jl:250-250`
- [[parallel.select_outer_route_select_outer_route_bang|select_outer_route!]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:571-571`

**Downstream**

- `callees` → [[parallel.outer_route_selection__route_elapsed_stats|_route_elapsed_stats]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:228-228`
- `callees` → [[parallel.outer_route_selection_outer_route_stats_snapshot|outer_route_stats_snapshot]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:241-241`
<!-- vulcan:connections:end -->

## Limitations
Holds `state.lock` for the whole iteration, so a large history entry briefly blocks concurrent feedback writers. The returned Dict is a copy; later updates to `state.history` are not reflected. `success_rate` uses `max(1, samples)` even though the `samples <= 0` branch already filtered, a redundant guard.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 213.
