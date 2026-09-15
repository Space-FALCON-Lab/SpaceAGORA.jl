---
id: gnc.replanning_rpo_retime_existing_plan
label: rpo_retime_existing_plan
kind: function
source:
  file: src/gnc/guidance/rpo/hypr_planning/replanning.jl
  symbol: rpo_retime_existing_plan
  lines:
  - 263
  - 263
inputs:
- id: plan
  type: RPOPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: current_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `current_rtn`.
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
- id: t
  type: Real
  units: n/a
  required: true
  description: Positional argument `t`.
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
  description: Return value of `rpo_retime_existing_plan`. Returns `rpo_plan_from_path(`.
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

# rpo_retime_existing_plan

## Purpose
Implements the `:retime` replanning action: keeps the geometric route of the current plan but regenerates the time schedule from the chaser's present position, used when clearance or tracking error is degraded but not enough to justify a full path search.

## Design & Implementation
Signature `rpo_retime_existing_plan(plan::RPOPlan, current_rtn, geometry, cfg::RPOPSOConfig, safe_distance_m::Real, t::Real)`. It computes `tail = rpo_remaining_reference_path(plan, current; sample_ds_m=cfg.sample_ds_m)`, overwrites the first column with `current` when the tail has at least two columns (a no-op in practice since the remaining-path routine already prepends `current`), and returns `rpo_plan_from_path(tail, geometry, cfg, safe_distance_m, t; cost=plan.cost, diagnostics=merge(plan.diagnostics, (replanning_action=:retime,)))`. The old plan is not mutated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RPOPlan | n/a | yes | Positional argument `plan`. |
| in | `current_rtn` | Any | n/a | yes | Positional argument `current_rtn`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Real | n/a | yes | Positional argument `safe_distance_m`. |
| in | `t` | Real | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_retime_existing_plan`. Returns `rpo_plan_from_path(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rpo_guidance_hooks_maybe_update_rpo_replanning_bang|maybe_update_rpo_replanning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:115-115`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`

**Downstream**

- `callees` → [[gnc.replanning_rpo_plan_from_path|rpo_plan_from_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:269-269`
- `callees` → [[gnc.replanning_rpo_remaining_reference_path|rpo_remaining_reference_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:265-265`
<!-- vulcan:connections:end -->

## Limitations
The retimed plan inherits `plan.cost` even though its path differs (a truncated tail), so cost is stale after retiming. The `size(tail, 2) >= 2` guard leaves a single-column tail untouched, meaning a chaser at the goal yields a one-point plan whose validity depends entirely on `rpo_reference_from_path` handling that case. It uses `cfg.sample_ds_m` rather than `config.remaining_sample_ds_m`, so the retime resolution can differ from the resolution used in the decision that triggered it.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr_planning/replanning.jl` line 263.
