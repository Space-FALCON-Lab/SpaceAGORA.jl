---
id: simulation.adaptive_routing__campaign_route_tuning
label: _campaign_route_tuning
kind: function
source:
  file: src/simulation/campaigns/adaptive_routing.jl
  symbol: _campaign_route_tuning
  lines:
  - 116
  - 116
inputs:
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
  type: OuterRouteTuning
  units: n/a
  description: Return value of `_campaign_route_tuning`.
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

# _campaign_route_tuning

## Purpose
Supplies the `OuterRouteTuning` thresholds used by campaign runners when the caller does not pass explicit `route_tuning`, centralising the decision that campaigns use unmodified defaults.

## Design & Implementation
Returns `OuterRouteTuning()` constructed with all default keyword values. The accompanying comment records the rationale: campaign runners can dispatch to a real `:process` backend (`ParallelProcess`), so the default thresholds (including `process_max_workers`, which defaults to `Sys.CPU_THREADS`) are appropriate without adjustment. The return type is annotated `::OuterRouteTuning`. It is consumed by `_run_monte_carlo_adaptive` and forwarded into `_campaign_route_plan` and `_record_campaign_route_feedback!`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | OuterRouteTuning | n/a | — | Return value of `_campaign_route_tuning`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.adaptive_routing__run_monte_carlo_adaptive|_run_monte_carlo_adaptive]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:262-262`
- [[simx.campaigns_constellation_ensemble_run_constellation_ensemble|run_constellation_ensemble]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:167-167`

**Downstream**

- `callees` → [[parallel.outer_route_state_outerroutetuning|OuterRouteTuning]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:119-119`
<!-- vulcan:connections:end -->

## Limitations
A fresh struct is allocated on every call rather than cached as a constant, which is negligible per campaign but not free. There is no environment-variable override at this layer, so users who want different campaign thresholds must pass `route_tuning` explicitly through `run_monte_carlo` or `run_constellation_ensemble`. Any change to `OuterRouteTuning` defaults silently changes campaign routing behaviour.

## Provenance
Mapped from `src/simulation/campaigns/adaptive_routing.jl` line 116.
