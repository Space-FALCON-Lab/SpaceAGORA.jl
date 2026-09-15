---
id: simx.campaigns_simulation_campaigns_simulationcampaigns
label: SimulationCampaigns
kind: struct
source:
  file: src/simulation/campaigns/simulation_campaigns.jl
  symbol: SimulationCampaigns
  lines:
  - 1
  - 22
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: engine_api
  type: Module
  units: n/a
  required: true
  description: RuntimeServices namespace and the SimulationEngine/SimulationCampaigns
    scope it aggregates, supplying the SPICE lock and the SimulationModel types used
    here.
- id: parallel_profiles
  type: Module
  units: n/a
  required: true
  description: ParallelProfiles namespace supplying OuterRouteFeatures, OuterRouteTuning,
    OuterRouteState and the select/record routing calls.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: campaign_api
  type: Module
  units: n/a
  description: Namespace exporting MonteCarloSpec, MonteCarloSampleResult, MonteCarloResult,
    run_monte_carlo, run_constellation_ensemble, campaign_route_features and campaign_outer_route_state.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simx
origin: agent
---
# SimulationCampaigns

## Purpose
`SimulationCampaigns` is the module wrapper for everything that runs many simulations rather than one. It owns no algorithm of its own: the file is an include manifest that pulls in the Monte Carlo driver, the adaptive outer routing policy and the constellation ensemble, then republishes their public names through a single export list.

## Model & Assumptions
The manifest fixes the load order deliberately. `monte_carlo.jl` defines the specification and result types, `adaptive_routing.jl` builds the bandit-style outer route selection on top of those types, and `constellation_ensemble.jl` consumes both. Reordering the includes breaks method resolution at precompile time because the later files reference types defined by the earlier ones at top level.

## Design & Implementation
Imports are explicit rather than blanket: `ParallelProfiles` is brought in whole and again by name for `OuterRouteFeatures`, `OuterRouteTuning`, `OuterRouteState`, `select_outer_route!` and `record_outer_route_feedback!`; `ParallelProcess` contributes `ProcessPool`, `campaign_process_pool` and `ensure_process_workers!` for the distributed route. `Distributed` is imported both as a module and by name for `remotecall_fetch` and `CachingPool`, which the process-pool campaign path uses to ship member work to worker processes. Seven names are exported, which is the entire supported campaign surface.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `parallel_profiles` | Module | n/a | yes | ParallelProfiles namespace supplying OuterRouteFeatures, OuterRouteTuning, OuterRouteState and the select/record routing calls. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `campaign_api` | Module | n/a | — | Namespace exporting MonteCarloSpec, MonteCarloSampleResult, MonteCarloResult, run_monte_carlo, run_constellation_ensemble, campaign_route_features and campaign_outer_route_state. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/campaigns/simulation_campaigns.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the module re-exports into the parent package scope, a name collision between a campaign export and another SpaceAGORA module surfaces as an ambiguity at using-time rather than at the definition site. The manifest carries no version guard, so adding a file without adding its include silently leaves the new code out of the build.

## Provenance
Mapped from `src/simulation/campaigns/simulation_campaigns.jl:1-22`; the three included files live alongside it in `src/simulation/campaigns/`.
