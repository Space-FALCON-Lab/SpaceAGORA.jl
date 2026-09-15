---
id: simulation.adaptive_routing__run_campaign_adaptive
label: _run_campaign_adaptive
kind: function
source:
  file: src/simulation/campaigns/adaptive_routing.jl
  symbol: _run_campaign_adaptive
  lines:
  - 232
  - 232
inputs:
- id: f
  type: Any
  units: n/a
  required: true
  description: Positional argument `f`.
- id: seeds
  type: Any
  units: n/a
  required: true
  description: Positional argument `seeds`.
- id: fail_fast
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `fail_fast`.
- id: features
  type: OuterRouteFeatures
  units: n/a
  required: true
  description: Keyword argument `features`.
- id: state
  type: OuterRouteState
  units: n/a
  required: true
  description: Keyword argument `state`.
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
  type: MonteCarloResult
  units: n/a
  description: Return value of `_run_campaign_adaptive`.
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

# _run_campaign_adaptive

## Purpose
Executes a seeded campaign end-to-end under adaptive outer routing: normalise features, plan the route, build a `MonteCarloSpec`, run it, and feed timing results back into the route state.

## Design & Implementation
Seeds are materialised with `collect(seeds)`; an empty list short-circuits to `MonteCarloResult(MonteCarloSampleResult[], 0.0, 0)`. `_campaign_features_for_routing` forces the `"montecarlo"` category with the true sample count, `_campaign_route_plan` yields `(route, threads, inner_thread_budget, record)`, and a `MonteCarloSpec(seeds=seed_values, threads=plan.threads, fail_fast=fail_fast)` is dispatched through `_run_campaign_with_route_env`. When `plan.record` is true, `_record_campaign_route_feedback!` updates `state` with the result. The return type is `MonteCarloResult`. All keyword arguments (`fail_fast`, `features`, `state`, `tuning`) are required.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | Any | n/a | yes | Positional argument `f`. |
| in | `seeds` | Any | n/a | yes | Positional argument `seeds`. |
| in | `fail_fast` | Bool | n/a | yes | Keyword argument `fail_fast`. |
| in | `features` | OuterRouteFeatures | n/a | yes | Keyword argument `features`. |
| in | `state` | OuterRouteState | n/a | yes | Keyword argument `state`. |
| in | `tuning` | OuterRouteTuning | n/a | yes | Keyword argument `tuning`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | MonteCarloResult | n/a | — | Return value of `_run_campaign_adaptive`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.adaptive_routing__run_monte_carlo_adaptive|_run_monte_carlo_adaptive]] · `callees` → `callers` · feedback · `src/simulation/campaigns/adaptive_routing.jl:263-263`
- [[simx.campaigns_constellation_ensemble_run_constellation_ensemble|run_constellation_ensemble]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:168-168`

**Downstream**

- `callees` → [[simulation.adaptive_routing__campaign_features_for_routing|_campaign_features_for_routing]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:242-242`
- `callees` → [[simulation.adaptive_routing__campaign_route_plan|_campaign_route_plan]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:243-243`
- `callees` → [[simulation.adaptive_routing__record_campaign_route_feedback_bang|_record_campaign_route_feedback!]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:247-247`
- `callees` → [[simulation.adaptive_routing__run_campaign_with_route_env|_run_campaign_with_route_env]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:245-245`
- `callees` → [[simulation.monte_carlo_result|MonteCarloResult]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:241-241`
- `callees` → [[simulation.monte_carlo_spec|MonteCarloSpec]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:244-244`
<!-- vulcan:connections:end -->

## Limitations
If `f` throws under `fail_fast`, the exception propagates before feedback is recorded, so a failing route never accumulates negative evidence. `collect(seeds)` materialises lazily generated seed ranges in memory. Feedback is recorded even when `plan.threads` was clamped to 1 by a small sample count, attributing serial timing to the chosen route label. The empty-seeds result reports `0` workers, which callers must not divide by.

## Provenance
Mapped from `src/simulation/campaigns/adaptive_routing.jl` line 232.
