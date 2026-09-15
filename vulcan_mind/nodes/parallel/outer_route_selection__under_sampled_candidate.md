---
id: parallel.outer_route_selection__under_sampled_candidate
label: _under_sampled_candidate
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _under_sampled_candidate
  lines:
  - 424
  - 424
inputs:
- id: candidates
  type: Vector{Symbol}
  units: n/a
  required: true
  description: Positional argument `candidates`.
- id: snapshot
  type: Any
  units: n/a
  required: true
  description: Positional argument `snapshot`.
- id: default_route
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `default_route`.
- id: min_samples
  type: Int
  units: n/a
  required: true
  description: Positional argument `min_samples`.
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
  description: Return value of `_under_sampled_candidate`.
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

# _under_sampled_candidate

## Purpose
Implements the exploration phase of the adaptive router: returns the first ranked candidate whose recorded sample count is below `min_samples`, forcing every feasible route to be tried before exploitation begins.

## Design & Implementation
`@inline _under_sampled_candidate(candidates, snapshot, default_route, min_samples::Int)::Union{Nothing, Symbol}` ranks candidates via `_route_ranked_candidates`, then for each route reads `get(snapshot, route, (samples=0, mean_s=Inf, success_rate=0.0))` and returns the route when `info.samples < max(1, min_samples)`. Returns `nothing` when all candidates are sufficiently sampled. `select_outer_route!` labels this outcome `"explore_hier"`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `candidates` | Vector{Symbol} | n/a | yes | Positional argument `candidates`. |
| in | `snapshot` | Any | n/a | yes | Positional argument `snapshot`. |
| in | `default_route` | Symbol | n/a | yes | Positional argument `default_route`. |
| in | `min_samples` | Int | n/a | yes | Positional argument `min_samples`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_under_sampled_candidate`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.select_outer_route_select_outer_route_bang|select_outer_route!]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:583-583`

**Downstream**

- `callees` → [[parallel.outer_route_selection__route_ranked_candidates|_route_ranked_candidates]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:430-430`
<!-- vulcan:connections:end -->

## Limitations
Exploration ignores observed failures, so a route that crashed on every prior attempt is still retried until `min_samples` is reached. Since the snapshot may come from a coarser fallback signature, the sample count can belong to a different workload class. `min_samples <= 0` is lifted to 1 silently.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 424.
