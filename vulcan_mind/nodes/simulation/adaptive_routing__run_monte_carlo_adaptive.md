---
id: simulation.adaptive_routing__run_monte_carlo_adaptive
label: _run_monte_carlo_adaptive
kind: function
source:
  file: src/simulation/campaigns/adaptive_routing.jl
  symbol: _run_monte_carlo_adaptive
  lines:
  - 252
  - 252
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
- id: route_features
  type: Union{Nothing, OuterRouteFeatures}
  units: n/a
  required: true
  description: Keyword argument `route_features`.
- id: route_state
  type: Union{Nothing, OuterRouteState}
  units: n/a
  required: true
  description: Keyword argument `route_state`.
- id: route_tuning
  type: Union{Nothing, OuterRouteTuning}
  units: n/a
  required: true
  description: Keyword argument `route_tuning`.
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
  description: Return value of `_run_monte_carlo_adaptive`.
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

# _run_monte_carlo_adaptive

## Purpose
Entry point used by `run_monte_carlo(threads=:auto)`: fills in default routing features, shared route state, and tuning for any `nothing` arguments, then delegates to `_run_campaign_adaptive`.

## Design & Implementation
Three ternaries resolve the optional inputs: `route_features === nothing` becomes `campaign_route_features(samples=0)` (a minimal montecarlo feature set that `_campaign_features_for_routing` later corrects with the real sample count), `route_state === nothing` becomes the process-global `campaign_outer_route_state()`, and `route_tuning === nothing` becomes `_campaign_route_tuning()`. The resolved values are passed as keywords to `_run_campaign_adaptive` along with `f`, `seeds`, and `fail_fast`. The return type is `MonteCarloResult`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | Any | n/a | yes | Positional argument `f`. |
| in | `seeds` | Any | n/a | yes | Positional argument `seeds`. |
| in | `fail_fast` | Bool | n/a | yes | Keyword argument `fail_fast`. |
| in | `route_features` | Union{Nothing, OuterRouteFeatures} | n/a | yes | Keyword argument `route_features`. |
| in | `route_state` | Union{Nothing, OuterRouteState} | n/a | yes | Keyword argument `route_state`. |
| in | `route_tuning` | Union{Nothing, OuterRouteTuning} | n/a | yes | Keyword argument `route_tuning`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | MonteCarloResult | n/a | — | Return value of `_run_monte_carlo_adaptive`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.run_monte_carlo|run_monte_carlo]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:285-285`

**Downstream**

- `callees` → [[simulation.adaptive_routing__campaign_route_tuning|_campaign_route_tuning]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:262-262`
- `callees` → [[simulation.adaptive_routing__run_campaign_adaptive|_run_campaign_adaptive]] · `callers` · feedback · `src/simulation/campaigns/adaptive_routing.jl:263-263`
- `callees` → [[simulation.campaign_route_features|campaign_route_features]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:260-260`
<!-- vulcan:connections:end -->

## Limitations
The default feature set has `n_sats=1`, `density_family="unknown"`, and `mission_time_s=0.0`, so campaigns that omit `route_features` all share one coarse statistics bucket regardless of workload cost. Using the global state means unrelated campaigns in the same session influence each other's routing. No validation is performed on `seeds` here; it is deferred to `_run_campaign_adaptive` and `MonteCarloSpec`.

## Provenance
Mapped from `src/simulation/campaigns/adaptive_routing.jl` line 252.
