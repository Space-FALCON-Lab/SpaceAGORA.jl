---
id: gnc.bridge_helpers__bridge_get_cnf
label: _bridge_get_cnf
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _bridge_get_cnf
  lines:
  - 209
  - 209
inputs:
- id: args
  type: Any
  units: n/a
  required: false
  description: Positional argument `args` (default `nothing`).
- id: cnf
  type: Any
  units: n/a
  required: false
  description: Keyword argument `cnf` (default `nothing`).
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
  description: Return value of `_bridge_get_cnf`. Returns `cnf` or `getproperty(args,
    :cnf)`.
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

# _bridge_get_cnf

## Purpose
Retrieves the aerobraking control-state object `cnf`, preferring an explicit keyword and falling back to a typed `args.cnf` field.

## Design & Implementation
Signature `_bridge_get_cnf(args=nothing; cnf=nothing)`. Returns `cnf` unchanged when it is not `nothing`; otherwise returns `getproperty(args, :cnf)` when `args` has that property. If neither is available it throws `ArgumentError("Control state `cnf` not found. Pass `cnf=` or typed args.cnf field.")`. The function is `@inline` and guarded by `isdefined` to tolerate multiple inclusion.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | no | Positional argument `args` (default `nothing`). |
| in | `cnf` | Any | n/a | no | Keyword argument `cnf` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_bridge_get_cnf`. Returns `cnf` or `getproperty(args, :cnf)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.second_switch_solver_second_time_switch_recalc_with_integration|second_time_switch_recalc_with_integration]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:4-4`
- [[gnc.targeting_solver__target_planning_impl|_target_planning_impl]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:13-13`
- [[gnc.targeting_solver_func_targeting_heatload|func_targeting_heatload]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:149-149`
- [[gnc.targeting_solver_func_targeting_num_int|func_targeting_num_int]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:121-121`
- [[gnc.tracking_executor__control_solarpanels_heatload_impl|_control_solarpanels_heatload_impl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:202-202`
- [[gnc.tracking_executor__control_solarpanels_openloop_impl|_control_solarpanels_openloop_impl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:270-270`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:11-11`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:9-9`
- [[gncx.e_edg_strategy_compute_e_edg_guidance_window_bang|compute_e_edg_guidance_window!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:12-12`
- [[gncx.energy_profile_solver_security_mode|security_mode]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/energy_profile_solver.jl:2-2`
- [[gncx.second_switch_solver_second_time_switch_recalc|second_time_switch_recalc]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:61-61`
- [[gncx.tracking_executor_control_solarpanels_heatrate|control_solarpanels_heatrate]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:85-85`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:12-12`
- [[gncy.t_edg_strategy_compute_t_edg_guidance_window_bang|compute_t_edg_guidance_window!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/t_edg_strategy.jl:3-3`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:12-12`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The dictionary-style `getindex` path used by `_bridge_optional_field` is not supported here, so a `Dict` containing `:cnf` still raises `ArgumentError`. No type check is applied to the returned object, so any value stored under `cnf` passes through. A `cnf` explicitly set to `nothing` in `args` is returned as `nothing` without error.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 209.
