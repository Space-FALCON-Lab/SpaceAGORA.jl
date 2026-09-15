---
id: gncy.planner_comparison_rpo_run_planner_comparison_batch
label: rpo_run_planner_comparison_batch
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_run_planner_comparison_batch
  lines:
  - 544
  - 646
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: comparison_request
  type: Tuple
  units: n/a
  required: true
  description: Comparison cases, station reference geometry, and the planner comparison
    configuration selecting planners and tracking settings.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: comparison_batch
  type: NamedTuple
  units: n/a
  description: Batch record holding cases, geometry, planner order, per-planner result
    rows, per-planner plans, and the resolved configuration.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncy
origin: agent
---

# rpo_run_planner_comparison_batch

## Purpose
`rpo_run_planner_comparison_batch` runs every configured RPO planner over every comparison case and returns one structured batch that downstream plotting and CSV writers consume. It is the top-level harness used to benchmark the HYPR PSO planner against CHOMP, STOMP, RRT-Connect, and RRT-Star on identical geometry.

## Model & Assumptions
Fairness across planners is enforced by two mechanisms. First, `_rpo_comparison_config_with_fixed_safe_distance` normalises the safe distance so every planner sees the same keep-out inflation. Second, when `cfg.optimizer.match_hypr_runtime` is set, HYPR is reordered to run first and each later planner is given a runtime limit equal to the HYPR plan time recorded for the same case, so the comparison measures solution quality at equal compute rather than at equal iteration count. Every planned path is then flown through the same LQ-MPC tracker so reported fuel, clearance, and final error come from closed-loop behaviour rather than the open-loop path.

## Design & Implementation
Randomness is reproducible through a `MersenneTwister` seeded from `cfg.rng_seed`; a per-run generator is derived by drawing from the master stream and offsetting by the planner index times one thousand plus the case index, so a case can be reproduced without replaying the whole batch. For each run the harness calls `rpo_plan_comparison_path`, then `rpo_track_retimed_path_lqmpc`, and merges planner-side metrics with tracking metrics into a flat row. Planner metrics include compute time, plan time, refinement time, iteration count, and cost. Tracking metrics include success, fuel used and its percentage, total and translational control effort, thrust saturation fraction, minimum clearance, keep-out violations, final position error, and both planned and actual travel duration. Progress is reported through `_rpo_comparison_progress_line!` before and after each run when enabled.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `comparison_request` | Tuple | n/a | yes | Comparison cases, station reference geometry, and the planner comparison configuration selecting planners and tracking settings. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `comparison_batch` | NamedTuple | n/a | — | Batch record holding cases, geometry, planner order, per-planner result rows, per-planner plans, and the resolved configuration. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[gnc.planner_comparison__rpo_comparison_config_with_fixed_safe_distance|_rpo_comparison_config_with_fixed_safe_distance]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:545-545`
- `callees` → [[gnc.planner_comparison__rpo_comparison_progress_line_bang|_rpo_comparison_progress_line!]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:557-557`
- `callees` → [[gnc.planner_comparison_rpo_plan_comparison_path|rpo_plan_comparison_path]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:575-575`
- `callees` → [[gnc.planner_comparison_rpo_track_retimed_path_lqmpc|rpo_track_retimed_path_lqmpc]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:584-584`
- `callees` → [[gnc.planner_comparison_rpoplannercomparisonconfig|RPOPlannerComparisonConfig]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:544-544`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:624-624`
- `callees` → [[gnc.trajectory_optimizers_normalize_rpo_comparison_planner_type|normalize_rpo_comparison_planner_type]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:547-547`
- `callees` → [[gnc.trajectory_optimizers_rpo_comparison_planner_label|rpo_comparison_planner_label]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:594-594`
<!-- vulcan:connections:end -->

## Limitations
The runtime-matching path requires HYPR to be among the selected planners; without it, every planner falls back to the configured runtime limit and the equal-compute comparison is lost. Per-run seeds are formed by unsigned addition, so two different planner and case index pairs can in principle collide on the same offset. Results are held entirely in memory as named tuples, so very large batches grow the working set linearly with cases times planners.

## Provenance
Mapped from planner_comparison.jl lines 544-646; include site observed at guidance_hooks.jl line 76.
