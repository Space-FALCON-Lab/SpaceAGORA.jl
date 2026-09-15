---
id: parallel.outer_route_selection__best_candidate_confidence
label: _best_candidate_confidence
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _best_candidate_confidence
  lines:
  - 477
  - 477
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
- id: exploration_c
  type: Float64
  units: n/a
  required: true
  description: Positional argument `exploration_c`.
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
  type: NamedTuple{(:route,
  units: n/a
  description: Return value of `_best_candidate_confidence`.
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

# _best_candidate_confidence

## Purpose
Exploitation rule of the adaptive router: selects the candidate with the lowest optimistic score `mean_s - width` (a lower-confidence-bound on elapsed time), reporting the chosen route, its confidence width, and the regret relative to the best observed mean.

## Theory & Math
$$\text{route}^* = \arg\min_r \left(\mu_r - w_r\right),\qquad \text{regret} = \mu_{\text{route}^*} - \min_r \mu_r$$ with $\mu_r$ the mean elapsed seconds of route $r$ and $w_r$ its confidence width.

## Design & Implementation
Signature `_best_candidate_confidence(candidates, snapshot, default_route, exploration_c::Float64)` returning `(route, confidence_s, regret_s)`. It first sums `total_samples` over candidates present in the snapshot (floored at 1), then iterates `_route_ranked_candidates(candidates, default_route)`, skipping routes with no samples or non-finite mean. For each it tracks `best_observed_mean`, computes `width` via `_candidate_confidence_width(info.std_s, info.samples, total_samples, exploration_c)`, and `score = mean_s - width`; a candidate wins on `score < best_score - 1e-12` or on an approximate tie with strictly higher success rate. Returns `(nothing, 0.0, 0.0)` if nothing qualifies, else `regret_s = max(0, best_mean - best_observed_mean)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `candidates` | Vector{Symbol} | n/a | yes | Positional argument `candidates`. |
| in | `snapshot` | Any | n/a | yes | Positional argument `snapshot`. |
| in | `default_route` | Symbol | n/a | yes | Positional argument `default_route`. |
| in | `exploration_c` | Float64 | n/a | yes | Positional argument `exploration_c`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NamedTuple{(:route, | n/a | — | Return value of `_best_candidate_confidence`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.select_outer_route_select_outer_route_bang|select_outer_route!]] · `callees` → `callers` · feedback · `src/parallel/routing/outer_route_selection.jl:588-588`

**Downstream**

- `callees` → [[parallel.outer_route_selection__candidate_confidence_width|_candidate_confidence_width]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:506-506`
- `callees` → [[parallel.outer_route_selection__route_ranked_candidates|_route_ranked_candidates]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:498-498`
- `callees` → [[parallel.select_outer_route_select_outer_route_bang|select_outer_route!]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:530-530`
<!-- vulcan:connections:end -->

## Limitations
Subtracting the width makes this an optimistic (lower-bound) selection, so a high-variance route can be chosen over a consistently faster one; `regret_s` exposes that cost but the caller only prints it under `tuning.trace`. Success rate influences the choice only on 1e-12 ties. Ranked iteration order silently breaks exact ties.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 477.
