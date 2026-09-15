---
id: parallel.outer_route_selection_outer_route_stats_snapshot
label: outer_route_stats_snapshot
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: outer_route_stats_snapshot
  lines:
  - 246
  - 246
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
  description: Return value of `outer_route_stats_snapshot`.
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

# outer_route_stats_snapshot

## Purpose
Public, reduced-field view of routing statistics for a signature, exposing `samples`, `mean_s`, and `success_rate` per route while hiding the `std_s` field that is internal to the confidence-bound selector.

## Design & Implementation
`outer_route_stats_snapshot(state::OuterRouteState, signature::String)` calls `_outer_route_stats_snapshot_internal`, then rebuilds a new `Dict{Symbol, NamedTuple{(:samples, :mean_s, :success_rate), ...}}` by copying the three public fields from each entry. The return type annotation fixes the tuple layout for downstream consumers such as diagnostics and tests.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `state` | OuterRouteState | n/a | yes | Positional argument `state`. |
| in | `signature` | String | n/a | yes | Positional argument `signature`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Dict{Symbol, | n/a | — | Return value of `outer_route_stats_snapshot`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__outer_route_stats_snapshot_internal|_outer_route_stats_snapshot_internal]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:241-241`

**Downstream**

- `callees` → [[parallel.outer_route_selection__outer_route_stats_snapshot_internal|_outer_route_stats_snapshot_internal]] · `callers` · feedback · `src/parallel/routing/outer_route_selection.jl:250-250`
<!-- vulcan:connections:end -->

## Limitations
Allocates a second Dict purely to drop one field; callers needing `std_s` must use the internal function. Locking is delegated to the internal call, so the returned snapshot can be stale immediately after return. Routes with zero recorded samples never appear, so callers cannot distinguish "never tried" from "unknown route".

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 246.
