---
id: parallel.outer_route_selection__candidate_confidence_width
label: _candidate_confidence_width
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _candidate_confidence_width
  lines:
  - 464
  - 464
inputs:
- id: std_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `std_s`.
- id: samples
  type: Int
  units: n/a
  required: true
  description: Positional argument `samples`.
- id: total_samples
  type: Int
  units: n/a
  required: true
  description: Positional argument `total_samples`.
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
  type: Float64
  units: n/a
  description: Return value of `_candidate_confidence_width`.
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

# _candidate_confidence_width

## Purpose
Computes the optimism bonus for one route in the UCB-style selector: a confidence half-width in seconds that shrinks as the route accumulates samples and grows slowly with the total number of observations across all routes.

## Theory & Math
$$w = c\,\sigma\sqrt{\frac{\ln\left(\max(2,\ N+1)\right)}{\max(1,\ n)}}$$ where $c$ is `adaptive_exploration_c`, $\sigma$ is the route's elapsed-time standard deviation in seconds, $n$ its sample count, and $N$ the total samples across candidates.

## Design & Implementation
`@inline _candidate_confidence_width(std_s::Float64, samples::Int, total_samples::Int, exploration_c::Float64)::Float64` returns `0.0` for `samples <= 0`, clamps `exploration = max(0, exploration_c)`, computes `scaled = std_s * sqrt(log(max(2, total_samples + 1)) / max(1, samples))`, and returns `max(0, exploration * scaled)` or `0.0` if that is not finite.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `std_s` | Float64 | n/a | yes | Positional argument `std_s`. |
| in | `samples` | Int | n/a | yes | Positional argument `samples`. |
| in | `total_samples` | Int | n/a | yes | Positional argument `total_samples`. |
| in | `exploration_c` | Float64 | n/a | yes | Positional argument `exploration_c`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_candidate_confidence_width`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__best_candidate_confidence|_best_candidate_confidence]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:506-506`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:472-472`
<!-- vulcan:connections:end -->

## Limitations
Scaling by `std_s` means a route whose few samples happened to be identical (zero variance) receives no exploration bonus, which can prematurely lock out a route. Non-finite `std_s` yields width 0 rather than an error. The formula is a heuristic UCB1-like bound, not a calibrated confidence interval.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 464.
