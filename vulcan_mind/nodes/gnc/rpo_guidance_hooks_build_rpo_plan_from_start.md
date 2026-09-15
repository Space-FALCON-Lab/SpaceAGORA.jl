---
id: gnc.rpo_guidance_hooks_build_rpo_plan_from_start
label: build_rpo_plan_from_start
kind: function
source:
  file: src/gnc/guidance/rpo/rpo_guidance_hooks.jl
  symbol: build_rpo_plan_from_start
  lines:
  - 16
  - 16
inputs:
- id: model
  type: RPOGuidanceModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: start_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `start_rtn`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: safe_distance_override
  type: Any
  units: n/a
  required: false
  description: Keyword argument `safe_distance_override` (default `nothing`).
- id: force_rrt_warmstart
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `force_rrt_warmstart` (default `false`).
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
  type: RPOPlan
  units: n/a
  description: Return value of `build_rpo_plan_from_start`. Returns `RPOPlan(`.
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

# build_rpo_plan_from_start

## Purpose
Core plan-construction routine: runs the particle-swarm path planner from an explicit RTN start position to the model's goal against a supplied geometry snapshot, then converts the resulting waypoint path into a time-tagged position and velocity reference wrapped in an `RPOPlan`.

## Design & Implementation
Signature is `(model, start_rtn, geometry, t::Float64; safe_distance_override=nothing, force_rrt_warmstart::Bool=false)`. The start is normalised to `SVector{3, Float64}`. The base configuration is `model.pso_config` or a default `RPOPSOConfig()` when that is `nothing`; `force_rrt_warmstart` rebuilds it through `rpo_pso_config(base_cfg; rrt_warmstart_enable=true)`. The safe distance is resolved by a three-level precedence: an explicit `safe_distance_override` wins, then `model.safe_distance_m` when strictly positive, then `base_cfg.safe_distance_m`. `rpo_pso_plan_path` produces the waypoint path and cost; `rpo_reference_from_path` turns it into `(t_ref, r_ref, v_ref)`. The returned `RPOPlan` is built with `valid=true` unconditionally and carries a rich `diagnostics` named tuple recording the cost components, adaptive-sampling record, whether refinement improved the solution, early-stop and iteration-timeout state with the iteration index and phase, the warmstart record, and `planned_at_s = t`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | RPOGuidanceModel | n/a | yes | Positional argument `model`. |
| in | `start_rtn` | Any | n/a | yes | Positional argument `start_rtn`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `safe_distance_override` | Any | n/a | no | Keyword argument `safe_distance_override` (default `nothing`). |
| in | `force_rrt_warmstart` | Bool | n/a | no | Keyword argument `force_rrt_warmstart` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPlan | n/a | — | Return value of `build_rpo_plan_from_start`. Returns `RPOPlan(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rpo_guidance_hooks_build_rpo_plan|build_rpo_plan]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:12-12`
- [[gnc.rpo_guidance_hooks_maybe_update_rpo_replanning_bang|maybe_update_rpo_replanning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:128-128`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:22-22`
- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:19-19`
- `callees` → [[gnc.rpo_plan_buffer_rpoplan|RPOPlan]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:36-36`
- `callees` → [[gncy.pso_parameters_rpopsoconfig|RPOPSOConfig]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:18-18`
- `callees` → [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:23-23`
- `callees` → [[gncz.rpo_reference_trajectory_rpo_reference_from_path|rpo_reference_from_path]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:30-30`
<!-- vulcan:connections:end -->

## Limitations
`valid=true` is hard-coded, so a planner run that terminated on an iteration timeout or that returned a path still intersecting the keep-out volume is still published as a valid plan; the caller must inspect `diagnostics.iteration_timed_out` to notice. The zero-valued sentinel used for safe distance means a deliberately configured `safe_distance_m = 0.0` on the model is indistinguishable from unset and silently falls through to the configuration default. Both the planner and the reference generator are called with no wall-clock budget enforced here, and any exception they raise propagates to the caller uncaught.

## Provenance
Mapped from `src/gnc/guidance/rpo/rpo_guidance_hooks.jl` line 16.
