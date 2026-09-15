---
id: gncy.switch_window_solver_switch_calculation_with_integration
label: switch_calculation_with_integration
kind: function
source:
  file: src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl
  symbol: switch_calculation_with_integration
  lines:
  - 3
  - 62
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: guidance_call
  type: Tuple
  units: n/a
  required: true
  description: Aerobraking guidance invocation carrying the input profile, mission
    model, pass position, argument dictionary, and time.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: time_switch
  type: Vector
  units: s
  description: Pair of drag-passage switch times produced by the heat-load root solve,
    or a saturated sentinel pair when no root brackets.
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

# switch_calculation_with_integration

## Purpose
Solves for the E-EDG drag-passage switch window by driving the integrated heat load of a candidate pass to the vehicle heat-load limit. It wraps a scalar residual around a full controlled aerobraking integration and finds the control gain that makes the accumulated heat load match the aerodynamic limit exactly.

## Theory & Math
The residual is $g(k)=Q_{\mathrm{final}}(k)-Q_{\lim}$ where $Q_{\mathrm{final}}$ is the heat load accumulated across the pass under gain $k=k_{cf}/100$. Bisection converges on the sign change of $g$ over $[k_{lo},k_{hi}]$ when $g(k_{hi})g(k_{lo})<0$.

## Model & Assumptions
The inner residual `func(k_cf)` rescales its argument by one hundredth, runs `_control_asim_ctrl` with the GRAM atmosphere enabled, reads the final accumulated heat value from the last element of the returned state history, and subtracts `m.aerodynamics.heat_load_limit`. The residual is assumed monotone enough over the bracket that bisection converges. Switch times are communicated out of the residual through globals rather than a return value, so the last successful integration determines the answer.

## Design & Implementation
Bracket limits are chosen per planet by name, with Mars starting at 3.3, Venus and Earth at 1.0, and Titan at 1.0e-5, all capped at 10.0. The function first probes the residual at 10 and 0.1. When the two probes straddle zero it calls `find_zero` with `Bisection()` at a relative tolerance of 1e-5. When both probes fall on the same side the routine returns a saturated pair chosen from `args[:heat_load_sol]`, mapping to either a zero-length or a one-thousand-second nominal window depending on which heat-load solution mode the mission requested.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `guidance_call` | Tuple | n/a | yes | Aerobraking guidance invocation carrying the input profile, mission model, pass position, argument dictionary, and time. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `time_switch` | Vector | s | — | Pair of drag-passage switch times produced by the heat-load root solve, or a saturated sentinel pair when no root brackets. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.e_edg_strategy_compute_e_edg_guidance_window_bang|compute_e_edg_guidance_window!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:24-24`

**Downstream**

- `callees` → [[gnc.guidance_hooks__control_asim_ctrl|_control_asim_ctrl]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:13-13`
- `callees` → [[gnc.heat_rate_models_func|func]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:6-6`
- `callees` → [[gnc.second_switch_solver_func|func]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:6-6`
- `callees` → [[gnc.switch_window_solver_func|func]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:6-6`
<!-- vulcan:connections:end -->

## Limitations
Planet dispatch is by string name, so an unmodelled planet leaves the bracket limits undefined and the call fails at the comparison. The saturated branches test `args[:heat_load_sol] == 0 && args[:heat_load_sol] == 3`, a condition that can never hold, so those particular fallbacks are unreachable. Reliance on module globals for `time_switch` makes the routine unsafe to call concurrently for two vehicles. Bisection cost scales with the cost of a full aerobraking pass integration, so each residual evaluation is expensive.

## Provenance
Mapped from switch_window_solver.jl lines 3-62; include site observed at guidance_hooks.jl line 86.
