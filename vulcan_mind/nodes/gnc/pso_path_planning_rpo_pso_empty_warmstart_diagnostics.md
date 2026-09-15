---
id: gnc.pso_path_planning_rpo_pso_empty_warmstart_diagnostics
label: rpo_pso_empty_warmstart_diagnostics
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: rpo_pso_empty_warmstart_diagnostics
  lines:
  - 111
  - 111
inputs:
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  description: Return value of `rpo_pso_empty_warmstart_diagnostics`. Returns `(`.
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

# rpo_pso_empty_warmstart_diagnostics

## Purpose
Provides the warm-start diagnostics record for runs where RRT seeding was disabled or skipped, so the planner result always has the same shape.

## Design & Implementation
Returns a named tuple with `enabled` reflecting `cfg.rrt_warmstart_enable`, `attempted` and `path_found` false, zero iterations and points, and infinite `cost` and `raw_cost`. Downstream reporting can therefore read the same fields whether or not RRT ran.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_pso_empty_warmstart_diagnostics`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_rpo_pso_rrt_warmstart_path|rpo_pso_rrt_warmstart_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:125-125`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:199-199`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Infinite costs in a diagnostics record can trip naive aggregation such as averaging warm-start cost across a campaign.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 111.
