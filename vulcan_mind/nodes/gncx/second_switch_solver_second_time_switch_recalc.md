---
id: gncx.second_switch_solver_second_time_switch_recalc
label: second_time_switch_recalc
kind: function
source:
  file: src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl
  symbol: second_time_switch_recalc
  lines:
  - 60
  - 171
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: GuidanceHooks namespace supplying the second-switch recalculation with
    the input parameters, mission, drag-passage indicator, argument dictionary, current
    time, and re-evaluation mode.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: switch_pair
  type: Tuple{Float64,Float64}
  units: s
  description: Unchanged first switch time and recomputed second switch time that
    brings the predicted passage heat load back onto the limit.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncx
origin: agent
---
# second_time_switch_recalc

## Purpose
`second_time_switch_recalc` re-solves the second angle-of-attack switch time mid-passage using the closed-form propagator. The first switch is treated as already committed, so the routine adjusts only the moment at which the panels return to the high-drag attitude, which is the one remaining degree of freedom once a passage is underway.

## Model & Assumptions
The prediction is assembled in two segments. A first closed-form call is made purely to obtain the trajectory time grid; a bang profile is then built on that grid by masking the angle to zero between the committed first switch and the candidate second switch and to the nominal mission angle elsewhere; and a second closed-form call propagates the passage with that profile so altitude, flight-path angle, and velocity are consistent with the control being evaluated. Density along the resulting altitude history comes from `density_exp`, and the molecular speed ratio is formed from the closed-form velocity and the thermal speed at the fixed planet temperature.

## Design & Implementation
Heat rate is then integrated separately over the two future segments. The interval from the current time to the candidate switch is flown at zero angle, and the interval beyond the switch at the mission angle, each evaluated through `heat_rate_calc` with the accommodation factor scaled by `args[:multiplicative_factor_heatload]`. Each segment is integrated as its mean rate times its span, which tolerates the non-uniform sample counts produced by the masks, and a degenerate case where both segments are empty returns the already-accumulated maximum load so the residual stays finite. When `args[:control_mode]` is three the predicted rate on the remaining segment is additionally capped at `args[:max_heat_rate]` by a masked blend, which models the fact that the heat-rate controller would itself intervene rather than letting the rate run free. The resulting predicted load is the residual driven to the mission limit, and the committed first switch is returned alongside the recomputed second.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | GuidanceHooks namespace supplying the second-switch recalculation with the input parameters, mission, drag-passage indicator, argument dictionary, current time, and re-evaluation mode. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `switch_pair` | Tuple{Float64,Float64} | s | — | Unchanged first switch time and recomputed second switch time that brings the predicted passage heat load back onto the limit. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.e_edg_strategy_compute_e_edg_guidance_window_bang|compute_e_edg_guidance_window!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:49-49`

**Downstream**

- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:61-61`
- `callees` → [[gnc.heat_rate_models_func|func]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:69-69`
- `callees` → [[gnc.second_switch_solver_func|func]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:69-69`
- `callees` → [[gnc.switch_window_solver_func|func]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:69-69`
- `callees` → [[gnc.tracking_executor_heat_rate_calc|heat_rate_calc]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:103-103`
- `callees` → [[gncx.closed_form_solution_closed_form|closed_form]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:72-72`
<!-- vulcan:connections:end -->

## Limitations
Both segments use a mean-rate approximation, which is coarse when the rate varies sharply across a segment. The routine trusts the committed first switch and therefore cannot recover from a poor initial schedule. `density_exp` is an exponential fit and does not reflect a sampled or GRAM atmosphere, so the integration-based sibling in the same file exists precisely for cases where that fidelity matters.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:60-171`.
