---
id: parallel.outer_route_selection__route_elapsed_stats
label: _route_elapsed_stats
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _route_elapsed_stats
  lines:
  - 203
  - 203
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
  type: NamedTuple{(:mean_s,
  units: n/a
  description: Return value of `_route_elapsed_stats`.
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

# _route_elapsed_stats

## Purpose
Derives the mean and population standard deviation of elapsed wall-clock seconds from the running sums stored in an `OuterRouteStats` record, feeding the confidence-width calculation of the UCB-style route selector.

## Theory & Math
$$\mu = \frac{1}{n}\sum_i t_i,\qquad \sigma = \sqrt{\max\left(0,\ \frac{1}{n}\sum_i t_i^2 - \mu^2\right)}$$ where $t_i$ are recorded elapsed times in seconds and $n = \max(1, \text{samples})$.

## Design & Implementation
`@inline _route_elapsed_stats(stats::OuterRouteStats)` returns the NamedTuple `(mean_s, std_s)`. It uses `samples = max(1, stats.samples)`, `mean_s = elapsed_sum_s / samples`, `mean_sq_s = elapsed_sq_sum_s / samples`, and `variance_s = max(0.0, mean_sq_s - mean_s^2)` before taking `sqrt`. The `max(0.0, ...)` guards against negative variance from floating-point cancellation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `stats` | OuterRouteStats | n/a | yes | Positional argument `stats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NamedTuple{(:mean_s, | n/a | — | Return value of `_route_elapsed_stats`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__outer_route_stats_snapshot_internal|_outer_route_stats_snapshot_internal]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:228-228`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The one-pass sum-of-squares formula loses precision when `mean_s^2` is large relative to the variance (catastrophic cancellation), which is why the clamp exists; a Welford update would be more robust. It computes the population rather than sample standard deviation, biasing `std_s` low for small `n`. With zero samples it returns `(0.0, 0.0)` rather than signalling absence.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 203.
