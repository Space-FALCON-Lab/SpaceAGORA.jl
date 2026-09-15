---
id: simulation.adaptive_routing__record_campaign_route_feedback_bang
label: _record_campaign_route_feedback!
kind: function
source:
  file: src/simulation/campaigns/adaptive_routing.jl
  symbol: _record_campaign_route_feedback!
  lines:
  - 204
  - 204
inputs:
- id: state
  type: OuterRouteState
  units: n/a
  required: true
  description: Positional argument `state`.
- id: features
  type: OuterRouteFeatures
  units: n/a
  required: true
  description: Positional argument `features`.
- id: route
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `route`.
- id: result
  type: MonteCarloResult
  units: n/a
  required: true
  description: Positional argument `result`.
- id: tuning
  type: OuterRouteTuning
  units: n/a
  required: true
  description: Keyword argument `tuning`.
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
  type: Nothing
  units: n/a
  description: Return value of `_record_campaign_route_feedback!`; mutates `state`
    in place.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _record_campaign_route_feedback!

## Purpose
Converts a finished `MonteCarloResult` into per-sample timing statistics and records them into the adaptive `OuterRouteState` under the campaign's feature signature, so later campaigns with the same shape can select the empirically fastest route.

## Design & Implementation
Returns `nothing` immediately when `length(result.samples) <= 0`. Otherwise it counts `successes = length(result.successful)` and `failures = length(result.failed)` and computes `per_sample_s = result.elapsed_s / n_samples`, the amortised campaign wall time per sample. It then calls `record_outer_route_feedback!(state, features; route, successes, failures, elapsed_success_s = per_sample_s * successes, elapsed_success_sq_sum_s = per_sample_s^2 * successes, tuning)`. Using amortised wall time rather than summed per-sample latencies credits a threaded route for throughput instead of penalising it for concurrency. The function mutates `state` and returns `Nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `state` | OuterRouteState | n/a | yes | Positional argument `state`. |
| in | `features` | OuterRouteFeatures | n/a | yes | Positional argument `features`. |
| in | `route` | Symbol | n/a | yes | Positional argument `route`. |
| in | `result` | MonteCarloResult | n/a | yes | Positional argument `result`. |
| in | `tuning` | OuterRouteTuning | n/a | yes | Keyword argument `tuning`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_record_campaign_route_feedback!`; mutates `state` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/campaigns/adaptive_routing.jl`
- [[simulation.adaptive_routing__run_campaign_adaptive|_run_campaign_adaptive]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:247-247`

**Downstream**

- `callees` → [[parcore.outer_route_metrics_record_outer_route_feedback_bang|record_outer_route_feedback!]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:219-219`
<!-- vulcan:connections:end -->

## Limitations
All successful samples are assigned the identical `per_sample_s`, so the recorded variance (`elapsed_success_sq_sum_s`) reflects zero within-campaign spread and understates true dispersion. Failed samples still consume wall time but only successes are credited, biasing the mean upward when failures are cheap. Nothing here is synchronised; concurrent campaigns sharing the global state depend on `record_outer_route_feedback!` for locking.

## Provenance
Mapped from `src/simulation/campaigns/adaptive_routing.jl` line 204.
