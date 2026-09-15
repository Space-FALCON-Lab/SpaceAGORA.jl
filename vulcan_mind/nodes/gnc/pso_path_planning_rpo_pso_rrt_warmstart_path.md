---
id: gnc.pso_path_planning_rpo_pso_rrt_warmstart_path
label: rpo_pso_rrt_warmstart_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: rpo_pso_rrt_warmstart_path
  lines:
  - 124
  - 124
inputs:
- id: start
  type: Any
  units: n/a
  required: true
  description: Positional argument `start`.
- id: goal
  type: Any
  units: n/a
  required: true
  description: Positional argument `goal`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: safe_distance_m
  type: Real
  units: n/a
  required: true
  description: Positional argument `safe_distance_m`.
- id: rng
  type: Any
  units: n/a
  required: true
  description: Positional argument `rng`.
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
  type: Any
  units: n/a
  description: 'Return value of `rpo_pso_rrt_warmstart_path`. Returns `plan.path_found
    ? plan.path : nothing, diagnostics`.'
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# rpo_pso_rrt_warmstart_path

## Purpose
Optionally runs RRT-Connect between start and goal to seed the swarm with a feasible polyline before PSO begins.

## Design & Implementation
Returns `nothing` plus empty diagnostics if warm start is disabled. Otherwise it assembles `RPORRTConnectSettings` from the six `rrt_warmstart_*` config fields, calls `rpo_rrt_connect_plan_path` with the safety distance, runtime limit and `rng`, and returns the path only when `plan.path_found`, together with diagnostics carrying the iteration count, costs and point count.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `start` | Any | n/a | yes | Positional argument `start`. |
| in | `goal` | Any | n/a | yes | Positional argument `goal`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Real | n/a | yes | Positional argument `safe_distance_m`. |
| in | `rng` | Any | n/a | yes | Positional argument `rng`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_pso_rrt_warmstart_path`. Returns `plan.path_found ? plan.path : nothing, diagnostics`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:198-198`

**Downstream**

- `callees` → [[gnc.pso_path_planning_rpo_pso_empty_warmstart_diagnostics|rpo_pso_empty_warmstart_diagnostics]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:125-125`
- `callees` → [[gnc.rrt_connect_rporrtconnectsettings|RPORRTConnectSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:126-126`
- `callees` → [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:134-134`
<!-- vulcan:connections:end -->

## Limitations
A failed RRT still consumes up to `rrt_warmstart_runtime_limit_s` of wall time before PSO starts, and the diagnostics record the attempt but the planner has no way to shorten its own iteration budget to compensate.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 124.
