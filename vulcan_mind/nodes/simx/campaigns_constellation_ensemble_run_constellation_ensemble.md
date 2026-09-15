---
id: simx.campaigns_constellation_ensemble_run_constellation_ensemble
label: run_constellation_ensemble
kind: function
source:
  file: src/simulation/campaigns/constellation_ensemble.jl
  symbol: run_constellation_ensemble
  lines:
  - 129
  - 202
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
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Full constellation configuration whose dynamics_model.spacecraft vector
    is split into one single-satellite member configuration per spacecraft.
- id: threads
  type: Union{Integer,Symbol}
  units: count
  required: true
  description: Outer worker count; a positive integer pins the Monte Carlo thread
    pool, and the symbol :auto hands routing to the adaptive bandit.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: ensemble_result
  type: MonteCarloResult
  units: n/a
  description: Aggregated per-member sample results produced by run_monte_carlo or
    by the adaptive campaign router.
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
# run_constellation_ensemble

## Purpose
`run_constellation_ensemble` turns one constellation configuration into an embarrassingly parallel ensemble: each spacecraft is propagated as its own single-satellite simulation, and the per-member results are collected through the Monte Carlo campaign machinery. It is the campaign-level shortcut for constellations whose members do not physically interact, trading the coupled flat-constellation right-hand side for outer parallelism across members.

## Model & Assumptions
The rewrite is only valid when members are genuinely uncoupled, and `_validate_ensemble_uncoupled` enforces that before any work starts: collision handling, inter-satellite interaction effectors and, unless `allow_gnc_effectors` is set, guidance/navigation/control effectors all reject the ensemble path because they read state belonging to other spacecraft. An empty spacecraft vector raises an `ArgumentError` rather than returning an empty result.

## Design & Implementation
Member configurations are built once by `_ensemble_member_configuration` with the tag `sat_<idx>_id_<sc.id>` so artifact paths stay distinct. Three execution routes exist. With `threads=:auto` the function computes route features via `campaign_route_features` (samples equal to the member count, `n_sats=1` because each routed sample propagates one satellite), then defers to `_run_campaign_adaptive`. With an integer thread count above one it builds a `MonteCarloSpec`, deepcopies each member inside the worker closure, sets `isolate_state=false` so isolation is paid exactly once instead of twice, and wraps the call in `withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1")` so the inner parallel policy knows an outer pool is already saturating cores. The single-worker route calls `run_simulation` directly and leaves `run_simulation`'s own `isolate_state` copy in place. Passing `route_state` or `route_tuning` without `threads=:auto` is rejected.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `args` | SimulationConfiguration | n/a | yes | Full constellation configuration whose dynamics_model.spacecraft vector is split into one single-satellite member configuration per spacecraft. |
| in | `threads` | Union{Integer,Symbol} | count | yes | Outer worker count; a positive integer pins the Monte Carlo thread pool, and the symbol :auto hands routing to the adaptive bandit. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `ensemble_result` | MonteCarloResult | n/a | — | Aggregated per-member sample results produced by run_monte_carlo or by the adaptive campaign router. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.constellation_ensemble__validate_ensemble_uncoupled|_validate_ensemble_uncoupled]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:69-69`

**Downstream**

- `callees` → [[simulation.adaptive_routing__campaign_route_tuning|_campaign_route_tuning]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:167-167`
- `callees` → [[simulation.adaptive_routing__run_campaign_adaptive|_run_campaign_adaptive]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:168-168`
- `callees` → [[simulation.campaign_route_features|campaign_route_features]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:165-165`
- `callees` → [[simulation.constellation_ensemble__ensemble_member_configuration|_ensemble_member_configuration]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:145-145`
- `callees` → [[simulation.constellation_ensemble__validate_ensemble_uncoupled|_validate_ensemble_uncoupled]] · `callers` · feedback · `src/simulation/campaigns/constellation_ensemble.jl:142-142`
- `callees` → [[simulation.monte_carlo_spec|MonteCarloSpec]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:182-182`
- `callees` → [[simulation.run_monte_carlo|run_monte_carlo]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:196-196`
- `callees` → [[simulation.run_simulation|run_simulation]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:160-160`
- `callees` → [[simx.engine_execution_run_simulation|run_simulation]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:160-160`
- `callees` → [[spaceagora.run_simulation|run_simulation]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:160-160`
<!-- vulcan:connections:end -->

## Limitations
Member splits alias the parent's environment and GNC model objects, and a solve mutates model state such as `planet.L_PI` in the planet-frame callback, so the per-member `deepcopy` inside the worker is mandatory on every concurrent route; removing it produces silent cross-member corruption rather than an error. Because each member is an independent solve, constellation-wide quantities such as relative geometry or conjunction screening are unavailable from the ensemble result.

## Provenance
Mapped from `src/simulation/campaigns/constellation_ensemble.jl:129-202`, with the uncoupling guard at line 52 and the member-configuration builder at line 32 of the same file.
