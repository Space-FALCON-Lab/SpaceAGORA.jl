---
id: gnc.rpo_guidance_hooks_maybe_update_rpo_replanning_bang
label: maybe_update_rpo_replanning!
kind: function
source:
  file: src/gnc/guidance/rpo/rpo_guidance_hooks.jl
  symbol: maybe_update_rpo_replanning!
  lines:
  - 83
  - 83
inputs:
- id: model
  type: RPOGuidanceModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: t
  type: Float64
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
  type: Bool
  units: n/a
  description: Return value of `maybe_update_rpo_replanning!`; mutates `model` in
    place. Returns `false` or `true`.
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

# maybe_update_rpo_replanning!

## Purpose
The in-flight replanning gate: each guidance step it re-evaluates the active plan against the current relative state and obstacle field, and either leaves the plan alone, retimes it, or rebuilds it from the present position.

## Design & Implementation
Four early exits return `false` before any work: no replanning configuration, `config.enabled` false, or an invalid `model.plan_buffer`. Otherwise it recomputes the RTN relative state through `_rpo_state_pos_vel` and `inertial_to_rtn_relative_state`, takes `start = x_rel[1:3]`, and asks `rpo_replanning_decision` for an action. A `rpo_replanning_signature` of `decision.spheres` identifies the obstacle configuration. On `decision.action == :none` the signature is stored and `replanning_persistence_count` is reset to zero. Otherwise the counter increments when the signature repeats and resets to one when it changes, implementing hysteresis: the action is suppressed until `replanning_persistence_count >= config.hysteresis_samples`, and again suppressed while `t - model.last_replanning_time_s < config.min_replan_interval_s`. Safe distance resolves in the order `config.safe_distance_m`, `model.safe_distance_m`, `base_cfg.safe_distance_m`, each gated on being strictly positive. The `:retime` branch calls `rpo_retime_existing_plan`, bumps `retime_count`, and records a `:retime` event. The `:replan` branch temporarily overwrites `model.goal_rtn` with `config.desired_goal_rtn` when present, calls `build_rpo_plan_from_start` with `force_rrt_warmstart=true`, and on success updates the buffer, bumps `replan_count`, and records `:replan`. Its `catch` restores the original goal when it was changed, bumps `replan_failure_count`, records `:replan_failed`, and warns rather than propagating, so a planner failure never aborts the simulation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | RPOGuidanceModel | n/a | yes | Positional argument `model`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `maybe_update_rpo_replanning!`; mutates `model` in place. Returns `false` or `true`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.rpo_guidance_hooks_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:158-158`

**Downstream**

- `callees` → [[core.reference_system_inertial_to_rtn_relative_state|inertial_to_rtn_relative_state]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:91-91`
- `callees` → [[gnc.replanning_rpo_replanning_signature|rpo_replanning_signature]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:94-94`
- `callees` → [[gnc.replanning_rpo_retime_existing_plan|rpo_retime_existing_plan]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:115-115`
- `callees` → [[gnc.rpo_guidance_hooks__rpo_record_replanning_event_bang|_rpo_record_replanning_event!]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:118-118`
- `callees` → [[gnc.rpo_guidance_hooks__rpo_replanning_config|_rpo_replanning_config]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:84-84`
- `callees` → [[gnc.rpo_guidance_hooks__rpo_state_pos_vel|_rpo_state_pos_vel]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:89-89`
- `callees` → [[gnc.rpo_guidance_hooks_build_rpo_plan_from_start|build_rpo_plan_from_start]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:128-128`
- `callees` → [[gncy.pso_parameters_rpopsoconfig|RPOPSOConfig]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:110-110`
- `callees` → [[gncy.replanning_rpo_replanning_decision|rpo_replanning_decision]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:93-93`
- `callees` → [[gncz.rpo_plan_buffer_update_rpo_plan_buffer_bang|update_rpo_plan_buffer!]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:116-116`
<!-- vulcan:connections:end -->

## Limitations
The goal restoration in the `catch` block only runs on failure: a successful replan leaves `model.goal_rtn` permanently overwritten with `config.desired_goal_rtn`, which is a silent state change that persists for every subsequent plan. The `try` has no `finally`, so an error thrown by `update_rpo_plan_buffer!` after the plan was built also skips restoration paths in the non-failure sense. Hysteresis counting is driven by call frequency rather than elapsed time, so `hysteresis_samples` means something different at different integrator step sizes. The retime branch has no failure handling at all, unlike the replan branch, so an exception there escapes to the integrator. The function mutates six model fields and is not thread-safe.

## Provenance
Mapped from `src/gnc/guidance/rpo/rpo_guidance_hooks.jl` line 83.
