---
id: simulation.campaign_route_features
label: campaign_route_features
kind: function
source:
  file: src/simulation/campaigns/adaptive_routing.jl
  symbol: campaign_route_features
  lines:
  - 64
  - 113
outputs:
- id: features
  type: NamedTuple
  units: n/a
  description: Measured campaign features used by adaptive outer-route selection.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
- parallel
- routing
charts:
- simulation
origin: agent
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
---

# campaign_route_features

## Purpose
`campaign_route_features` extracts the measurable properties that adaptive campaign routing needs before selecting an execution route. It converts campaign configuration and recent runtime context into a stable feature record, keeping route policy independent from the detailed Monte Carlo and simulation result structures.

## Theory & Math
The output is a feature vector `φ` containing discrete configuration values and measured quantities such as worker demand, sample count, and route history. A route policy evaluates `φ` against candidate tuning rules. These features describe execution cost and availability; they are not spacecraft state variables and do not enter the ODE directly.

## Model & Assumptions
Feature fields must use consistent units and names across route candidates. Missing measurements are expected to receive explicit sentinel values or trigger the default route. The function assumes the campaign metadata reflects the work that will actually be executed, especially sample count and worker requirements.

## Design & Implementation
`adaptive_routing.jl` gathers configuration values, campaign dimensions, and current route context into the returned named tuple. The Monte Carlo driver passes the result to `select_outer_route!`, which records the selected route and uses it to determine process-worker reconciliation. Keeping feature assembly separate makes the policy testable with synthetic campaign inputs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `features` | NamedTuple | n/a | — | Measured campaign features used by adaptive outer-route selection. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.adaptive_routing__campaign_density_family|_campaign_density_family]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:41-41`
- [[simulation.adaptive_routing__run_monte_carlo_adaptive|_run_monte_carlo_adaptive]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:260-260`
- [[simx.campaigns_constellation_ensemble_run_constellation_ensemble|run_constellation_ensemble]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:165-165`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:75-75`
- `callees` → [[parallel.outer_route_state_outerroutefeatures|OuterRouteFeatures]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:72-72`
- `callees` → [[simulation.adaptive_routing__campaign_density_family|_campaign_density_family]] · `callers` · feedback · `src/simulation/campaigns/adaptive_routing.jl:98-98`
<!-- vulcan:connections:end -->

## Limitations
Feature extraction cannot predict resource contention created by unrelated processes or later route changes. Historical feedback may be stale, and a feature record can remain structurally valid while describing a workload with different native-library behavior. The adaptive selector must therefore retain a fallback path.

## Provenance
Mapped from `src/simulation/campaigns/adaptive_routing.jl:41-113`.
