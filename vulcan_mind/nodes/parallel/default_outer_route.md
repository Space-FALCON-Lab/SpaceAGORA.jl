---
id: parallel.default_outer_route
label: default_outer_route
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: default_outer_route
  lines:
  - 312
  - 361
outputs:
- id: route
  type: Symbol
  units: n/a
  description: Deterministic fallback route used when adaptive selection has no usable
    evidence.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
- routing
charts:
- parallel
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

# default_outer_route

## Purpose
`default_outer_route` supplies the fallback execution route for a campaign when no adaptive feedback is available or when selection is disabled. It gives the rest of the routing system a stable answer and prevents an empty route state from reaching worker setup.

## Theory & Math
The decision is categorical. A route is selected from the configured route vocabulary, so the function can be viewed as `r₀ = arg default`, where `r₀` is the configured fallback symbol. It does not estimate throughput or alter physical state.

## Model & Assumptions
The returned symbol must be accepted by `outer_route_candidates` and by the process or thread execution layer. The fallback is assumed to be safe for the host’s available resources, but the function does not query the scheduler or launch a worker to prove that assumption.

## Design & Implementation
`outer_route_selection.jl` reads the profile and route-state inputs required by the default policy and returns the configured route. `select_outer_route!` calls this path when adaptive features are absent, invalid, or explicitly bypassed. Keeping the fallback separate makes route selection deterministic in tests and gives the campaign layer a clear policy boundary.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `route` | Symbol | n/a | — | Deterministic fallback route used when adaptive selection has no usable evidence. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.outer_route_selection__priority_outer_route_montecarlo|_priority_outer_route_montecarlo]] · `callees` → `callers` · feedback · `src/parallel/routing/outer_route_selection.jl:307-307`
- [[parallel.select_outer_route_select_outer_route_bang|select_outer_route!]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:543-543`

**Downstream**

- `callees` → [[parallel.outer_route_selection__feature_is_lightweight|_feature_is_lightweight]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:337-337`
- `callees` → [[parallel.outer_route_selection__is_native_gram_point_density|_is_native_gram_point_density]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:323-323`
- `callees` → [[parallel.outer_route_selection__priority_outer_route_montecarlo|_priority_outer_route_montecarlo]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:328-328`
- `callees` → [[parallel.outer_route_selection__threads_or_none|_threads_or_none]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:348-348`
- `callees` → [[parallel.outer_route_state_outerroutetuning|OuterRouteTuning]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:314-314`
<!-- vulcan:connections:end -->

## Limitations
The default route can be conservative or inefficient for a particular machine. It cannot react to queue pressure, memory pressure, or observed task duration. If the route vocabulary changes without updating the default, downstream selection fails at symbol validation rather than silently translating the old value.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl:305-361`.
