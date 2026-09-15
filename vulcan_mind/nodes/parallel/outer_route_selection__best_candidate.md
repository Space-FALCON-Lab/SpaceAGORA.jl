---
id: parallel.outer_route_selection__best_candidate
label: _best_candidate
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _best_candidate
  lines:
  - 440
  - 440
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
  description: Return value of `_best_candidate`.
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

# _best_candidate

## Purpose
Greedy exploitation rule that picks the candidate route with the lowest recorded mean elapsed time, breaking near-ties by higher success rate; it is the non-confidence-adjusted counterpart of `_best_candidate_confidence`.

## Design & Implementation
`@inline _best_candidate(candidates::Vector{Symbol}, snapshot)::Union{Nothing, Symbol}` iterates candidates in the given order, skipping entries with `samples <= 0` or non-finite `mean_s`. It replaces the incumbent when `info.mean_s < best_mean - 1e-12`, or when `isapprox(info.mean_s, best_mean; atol=1e-12, rtol=0.0)` and `info.success_rate > best_success_rate`. Returns `nothing` when no candidate has usable history.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `candidates` | Vector{Symbol} | n/a | yes | Positional argument `candidates`. |
| in | `snapshot` | Any | n/a | yes | Positional argument `snapshot`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_best_candidate`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Not called by `select_outer_route!` in the current file, which uses the UCB variant instead; it remains available for callers wanting pure greedy selection. Success rate only matters on an exact tie within 1e-12 s, so a route that is 1 ms faster but fails half the time still wins. Candidate order acts as the final tie-breaker.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 440.
