---
id: parcore.parallel_profiles_parallelprofiles
label: ParallelProfiles
kind: struct
source:
  file: src/parallel/routing/parallel_profiles.jl
  symbol: ParallelProfiles
  lines:
  - 2
  - 20
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: 'Routing files included by this module: profile definitions, environment
    mapping, outer route state, metrics and selection.'
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: routing_api
  type: Module
  units: n/a
  description: 'Outer routing surface: profile types and parsing, profile configuration
    and environment mapping, route state and statistics, route selection, feedback
    recording and state persistence.'
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parcore
origin: agent
---

# ParallelProfiles

## Purpose
`ParallelProfiles` is the module that owns outer routing. It includes the routing files and exports the profile abstraction, the outer-route state model and the selection and feedback entry points, giving the campaign layer one place to ask how a workload should be parallelised at the top level.

## Model & Assumptions
The module separates three concerns that the export list makes explicit. A profile is a user-facing preset that maps to environment configuration. A route is the runtime decision between no outer parallelism, threads or processes. The state is the accumulated evidence that lets the route decision adapt. Selection reads the state and the tuning, execution reports back through feedback, and persistence carries the evidence across runs.

## Design & Implementation
The file is a twenty-line aggregator whose export lines follow those three groups. Profile exports are `ParallelProfile`, `ParallelProfileConfig`, `parse_parallel_profile`, `parallel_profile_name`, `profile_config`, `profile_env_pairs` and `with_parallel_profile`. State exports are `OuterRouteFeatures`, `OuterRouteTuning`, `OuterRouteState`, `reset_outer_route_state!`, `outer_route_signature` and `outer_route_stats_snapshot`. Decision exports are `default_outer_route`, `outer_route_candidates`, `select_outer_route!` and `record_outer_route_feedback!`, with `load_outer_route_state!` and `save_outer_route_state` covering persistence. This module is the entry point the master chart names for the parallel subsystem.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Routing files included by this module: profile definitions, environment mapping, outer route state, metrics and selection. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `routing_api` | Module | n/a | — | Outer routing surface: profile types and parsing, profile configuration and environment mapping, route state and statistics, route selection, feedback recording and state persistence. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/parallel_profiles.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
As with the policy module, the included files share one flat namespace, so internal helpers are reachable from anywhere inside it and the real dependency order is implicit in the include sequence. The exported surface mixes the stable profile API with routing internals that exist mainly for tests. Nothing here enforces that a caller which selects a route also reports feedback for it, so a code path that forgets the feedback call silently stops contributing evidence.

## Provenance
Mapped from `src/parallel/routing/parallel_profiles.jl:2-20`.
