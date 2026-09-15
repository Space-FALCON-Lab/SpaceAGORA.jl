---
id: parallel.outer_route_selection_outer_route_candidates
label: outer_route_candidates
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: outer_route_candidates
  lines:
  - 371
  - 371
inputs:
- id: f
  type: OuterRouteFeatures
  units: n/a
  required: true
  description: Positional argument `f`.
- id: tuning
  type: OuterRouteTuning
  units: n/a
  required: false
  description: Keyword argument `tuning` (default `OuterRouteTuning()`).
- id: machine_class
  type: Symbol
  units: n/a
  required: false
  description: Keyword argument `machine_class` (default `:small`).
- id: threads_available
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `threads_available` (default `true`).
- id: parallel_enabled
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `parallel_enabled` (default `true`).
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
  type: Vector{Symbol}
  units: n/a
  description: Return value of `outer_route_candidates`.
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

# outer_route_candidates

## Purpose
Public function enumerating the feasible outer parallel routes for a workload, forming the action set that the adaptive selector explores and ranks; it also enforces hard constraints such as GRAM lock limitations and the absence of threads.

## Design & Implementation
Signature `outer_route_candidates(f::OuterRouteFeatures; tuning=OuterRouteTuning(), machine_class=:small, threads_available=true, parallel_enabled=true)::Vector{Symbol}`. Returns `[:none]` when parallelism is disabled, `[:none, :process]` for native GRAM point density, and otherwise starts from `[:none]`, appends `:threads` if `threads_available`, and appends `:process` when `allow_process` holds, where `allow_process` is `_priority_outer_route_montecarlo(...) == :process` for the "montecarlo" category and `_feature_heavy_for_process(f, tuning)` otherwise. The final `unique` guards against duplicates.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | OuterRouteFeatures | n/a | yes | Positional argument `f`. |
| in | `tuning` | OuterRouteTuning | n/a | no | Keyword argument `tuning` (default `OuterRouteTuning()`). |
| in | `machine_class` | Symbol | n/a | no | Keyword argument `machine_class` (default `:small`). |
| in | `threads_available` | Bool | n/a | no | Keyword argument `threads_available` (default `true`). |
| in | `parallel_enabled` | Bool | n/a | no | Keyword argument `parallel_enabled` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{Symbol} | n/a | — | Return value of `outer_route_candidates`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.select_outer_route_select_outer_route_bang|select_outer_route!]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:554-554`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:387-387`
- `callees` → [[parallel.outer_route_selection__feature_heavy_for_process|_feature_heavy_for_process]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:398-398`
- `callees` → [[parallel.outer_route_selection__is_native_gram_point_density|_is_native_gram_point_density]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:381-381`
- `callees` → [[parallel.outer_route_selection__priority_outer_route_montecarlo|_priority_outer_route_montecarlo]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:390-390`
- `callees` → [[parallel.outer_route_state_outerroutetuning|OuterRouteTuning]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:373-373`
<!-- vulcan:connections:end -->

## Limitations
`:none` is always a candidate, so the adaptive selector will spend `adaptive_min_samples` runs on serial execution even for obviously heavy workloads. The `unique` call cannot trigger given the construction and merely allocates. Ordering of the returned vector is `[:none, :threads, :process]`, which `_route_ranked_candidates` later reorders.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 371.
