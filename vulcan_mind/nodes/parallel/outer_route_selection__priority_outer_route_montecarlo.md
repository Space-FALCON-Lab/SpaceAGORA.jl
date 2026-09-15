---
id: parallel.outer_route_selection__priority_outer_route_montecarlo
label: _priority_outer_route_montecarlo
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _priority_outer_route_montecarlo
  lines:
  - 285
  - 285
inputs:
- id: f
  type: OuterRouteFeatures
  units: n/a
  required: true
  description: Positional argument `f`.
- id: t
  type: OuterRouteTuning
  units: n/a
  required: true
  description: Positional argument `t`.
- id: machine_class
  type: Symbol
  units: n/a
  required: true
  description: Keyword argument `machine_class`.
- id: threads_available
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `threads_available`.
- id: parallel_enabled
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `parallel_enabled`.
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
  type: Symbol
  units: n/a
  description: Return value of `_priority_outer_route_montecarlo`.
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

# _priority_outer_route_montecarlo

## Purpose
Routing rule specialised for Monte Carlo campaigns, where the natural parallelism is across independent samples rather than within one simulation, choosing between `:none`, `:threads`, and `:process` based on sample count, mission length, and machine class.

## Design & Implementation
Signature `_priority_outer_route_montecarlo(f, t; machine_class::Symbol, threads_available::Bool, parallel_enabled::Bool)::Symbol`. Returns `:none` if `!parallel_enabled` or `f.montecarlo_samples <= 1`; returns `:process` if `machine_class in (:large, :medium)` and either `f.montecarlo_samples >= t.mc_process_min_samples` or `f.mission_time_s >= t.mc_process_min_mission_s`; otherwise `_threads_or_none(threads_available)`. `default_outer_route` calls it when `f.category` is "montecarlo", and `outer_route_candidates` uses `== :process` to decide process eligibility.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | OuterRouteFeatures | n/a | yes | Positional argument `f`. |
| in | `t` | OuterRouteTuning | n/a | yes | Positional argument `t`. |
| in | `machine_class` | Symbol | n/a | yes | Keyword argument `machine_class`. |
| in | `threads_available` | Bool | n/a | yes | Keyword argument `threads_available`. |
| in | `parallel_enabled` | Bool | n/a | yes | Keyword argument `parallel_enabled`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_priority_outer_route_montecarlo`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.default_outer_route|default_outer_route]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:328-328`
- [[parallel.outer_route_selection_outer_route_candidates|outer_route_candidates]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:390-390`

**Downstream**

- `callees` → [[parallel.default_outer_route|default_outer_route]] · `callers` · feedback · `src/parallel/routing/outer_route_selection.jl:307-307`
- `callees` → [[parallel.outer_route_selection__threads_or_none|_threads_or_none]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:303-303`
<!-- vulcan:connections:end -->

## Limitations
A `:small` machine never receives `:process` regardless of sample count. The two tuning thresholds are OR-ed, so a long mission with only two samples still triggers process mode on a large machine, spawning workers for minimal gain. Category matching is done by the caller via `lowercase(strip(f.category)) == "montecarlo"`, so variants such as "monte_carlo" bypass this rule.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 285.
