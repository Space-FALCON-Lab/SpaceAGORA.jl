---
id: simulation.adaptive_routing__campaign_route_plan
label: _campaign_route_plan
kind: function
source:
  file: src/simulation/campaigns/adaptive_routing.jl
  symbol: _campaign_route_plan
  lines:
  - 136
  - 136
inputs:
- id: features
  type: OuterRouteFeatures
  units: n/a
  required: true
  description: Positional argument `features`.
- id: n_samples
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_samples`.
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
  type: Tuple
  units: n/a
  description: Return value of `_campaign_route_plan`. Returns `(route=:none, threads=1,
    inner_thread_budget=1, record=false)` or `(route=route, threads=workers, inner_thread_budget=inner_thread_budget,
    record=t`.
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

# _campaign_route_plan

## Purpose
Chooses the outer execution route (`:none`, `:threads`, or `:process`) and worker sizing for a campaign of `n_samples` runs, delegating the choice to the adaptive `select_outer_route!` and returning a `NamedTuple` plan that `_run_campaign_with_route_env` executes.

## Design & Implementation
If `ParallelPolicy.outer_parallel_active()` reports that an enclosing campaign already owns the outer split, the function returns `(route=:none, threads=1, inner_thread_budget=1, record=false)` so nested workers do not oversubscribe the pool or poison shared statistics. Otherwise `select_outer_route!(state, features; tuning, machine_class=ParallelProfiles._machine_parallel_class(), threads_available=nthreads()>1, parallel_enabled=true)` picks the route. Worker count is `min(n_samples, tuning.process_max_workers)` for `:process`, `min(n_samples, Threads.nthreads())` for `:threads`, else 1. `inner_thread_budget = max(1, fld(nthreads(), workers))` splits the pool between outer workers and per-sample inner parallelism. `record=true` on this path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `features` | OuterRouteFeatures | n/a | yes | Positional argument `features`. |
| in | `n_samples` | Int | n/a | yes | Positional argument `n_samples`. |
| in | `state` | OuterRouteState | n/a | yes | Keyword argument `state`. |
| in | `tuning` | OuterRouteTuning | n/a | yes | Keyword argument `tuning`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_campaign_route_plan`. Returns `(route=:none, threads=1, inner_thread_budget=1, record=false)` or `(route=route, threads=workers, inner_thread_budget=inner_thread_budget, record=t`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/campaigns/adaptive_routing.jl`
- [[simulation.adaptive_routing__run_campaign_adaptive|_run_campaign_adaptive]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:243-243`

**Downstream**

- `callees` → [[parallel.env_mapping__machine_parallel_class|_machine_parallel_class]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:154-154`
- `callees` → [[parallel.select_outer_route_select_outer_route_bang|select_outer_route!]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:150-150`
<!-- vulcan:connections:end -->

## Limitations
`select_outer_route!` mutates `state` (its `!` contract) even though this function itself has no `!`. Process workers are sized from `tuning.process_max_workers` rather than measured CPU availability, so on a loaded machine they can still oversubscribe. `inner_thread_budget` is computed from `Threads.nthreads()` even for the `:process` route, where workers run single-threaded; the value is then unused for that route. The machine-class query is repeated on each call.

## Provenance
Mapped from `src/simulation/campaigns/adaptive_routing.jl` line 136.
