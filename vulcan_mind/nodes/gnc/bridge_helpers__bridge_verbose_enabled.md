---
id: gnc.bridge_helpers__bridge_verbose_enabled
label: _bridge_verbose_enabled
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _bridge_verbose_enabled
  lines:
  - 6
  - 6
inputs:
- id: args
  type: Any
  units: n/a
  required: false
  description: Positional argument `args` (default `nothing`).
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
  description: Return value of `_bridge_verbose_enabled`.
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

# _bridge_verbose_enabled

## Purpose
Decides whether the legacy aerobraking control bridge should emit verbose diagnostics, merging an environment-variable override with whatever verbosity flag the runtime `args` object carries.

## Design & Implementation
Returns `true` immediately when `ENV["SPACEAGORA_DEBUG_LEGACY_CONTROL"]` equals the string `"1"`. Otherwise it probes `args` (default `nothing`) with `hasproperty`: first `args.simulation_settings.verbose`, then a top-level `args.verbose`, converting whichever is found with `Bool(...)`. When `args` is `nothing` or has neither field the result is `false`. The function is `@inline` and wrapped in an `isdefined` guard so re-including the file does not redefine it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | no | Positional argument `args` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_bridge_verbose_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.eom_predictor_shooting_residual_bang|shooting_residual!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1037-1037`
- [[gnc.targeting_solver__target_planning_impl|_target_planning_impl]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:22-22`
- [[gnc.targeting_solver_control_solarpanels_targeting_heatload|control_solarpanels_targeting_heatload]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:140-140`
- [[gnc.targeting_solver_control_solarpanels_targeting_num_int|control_solarpanels_targeting_num_int]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:117-117`
- [[gnc.tracking_executor__control_exception_fallback|_control_exception_fallback]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:16-16`
- [[gnc.tracking_executor_df|df]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:134-134`
- [[gnc.tracking_executor_f|f]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:59-59`
- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:203-203`
- [[gncx.tracking_executor_control_solarpanels_heatrate|control_solarpanels_heatrate]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:134-134`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1037-1037`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:203-203`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`Bool(getproperty(...))` throws `InexactError` if the stored verbose flag is a non-0/1 number, and `MethodError` if it is a string. The environment variable is re-read from `ENV` on every call, which is a dictionary lookup on a hot path when callers invoke this per step. Only the exact string `"1"` enables the override; `"true"` or `"yes"` are ignored.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 6.
