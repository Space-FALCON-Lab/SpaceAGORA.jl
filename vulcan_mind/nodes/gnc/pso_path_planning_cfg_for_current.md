---
id: gnc.pso_path_planning_cfg_for_current
label: cfg_for_current
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: cfg_for_current
  lines:
  - 285
  - 285
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
  type: Any
  units: n/a
  description: Return value of `cfg_for_current`. Returns `rpo_pso_config(cfg; n_waypoints=current_n_waypoints,
    search_margin_m=current_sea`.
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

# cfg_for_current

## Purpose
Rebuilds the working config to reflect the waypoint count and search margin that reexploration may have grown since the run began.

## Design & Implementation
A closure calling `rpo_pso_config(cfg; n_waypoints=current_n_waypoints, search_margin_m=current_search_margin)`. It is used by `reset_swarm!` to size bounds and by the final refinement so the reported `config` matches the geometry actually optimised.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `cfg_for_current`. Returns `rpo_pso_config(cfg; n_waypoints=current_n_waypoints, search_margin_m=current_sea`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:611-611`
- [[gnc.pso_path_planning_reset_swarm_bang|reset_swarm!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:324-324`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:285-285`

**Downstream**

- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:286-286`
<!-- vulcan:connections:end -->

## Limitations
Every call allocates a new config struct; it is called rarely enough that this does not matter.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 285.
