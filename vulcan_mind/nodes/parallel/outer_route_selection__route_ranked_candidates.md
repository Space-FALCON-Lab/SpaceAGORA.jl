---
id: parallel.outer_route_selection__route_ranked_candidates
label: _route_ranked_candidates
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _route_ranked_candidates
  lines:
  - 406
  - 406
inputs:
- id: candidates
  type: Vector{Symbol}
  units: n/a
  required: true
  description: Positional argument `candidates`.
- id: default_route
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `default_route`.
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
  description: Return value of `_route_ranked_candidates`.
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

# _route_ranked_candidates

## Purpose
Imposes a deterministic exploration order on the candidate routes: the heuristic default first, then `:threads`, `:none`, `:process`, then any remaining candidates, so under-sampling checks and tie-breaking favour cheaper-to-try routes.

## Design & Implementation
`@inline _route_ranked_candidates(candidates::Vector{Symbol}, default_route::Symbol)::Vector{Symbol}` builds `ranked` by pushing `default_route` if present, then each of `(:threads, :none, :process)` not already ranked, then any leftover candidate in original order. Membership checks use linear `in` on vectors of at most three or four symbols. Consumed by `_under_sampled_candidate` and `_best_candidate_confidence`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `candidates` | Vector{Symbol} | n/a | yes | Positional argument `candidates`. |
| in | `default_route` | Symbol | n/a | yes | Positional argument `default_route`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{Symbol} | n/a | — | Return value of `_route_ranked_candidates`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__best_candidate_confidence|_best_candidate_confidence]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:498-498`
- [[parallel.outer_route_selection__under_sampled_candidate|_under_sampled_candidate]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:430-430`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:409-409`
<!-- vulcan:connections:end -->

## Limitations
The fixed preference `:threads` before `:none` before `:process` is hard-coded, not tunable. Because `_best_candidate_confidence` iterates in this order and only replaces on a strict improvement, rank order acts as a hidden tie-breaker. Allocates a new vector per call.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 406.
