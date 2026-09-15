---
id: gnc.switch_window_solver_func
label: func
kind: function
source:
  file: src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl
  symbol: func
  lines:
  - 6
  - 6
inputs:
- id: k_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `k_cf`.
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
  description: Return value of `func`. Returns `Q - m.aerodynamics.heat_load_limit`.
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

# func

## Purpose
Inner residual function of `switch_calculation_with_integration`. For a candidate control gain it runs the full controlled drag-passage integration and returns how far the accumulated heat load overshoots or undershoots the vehicle's heat-load limit, giving the bracketing root solver a scalar to drive to zero.

## Design & Implementation
`func(k_cf)` rescales its argument as `k = k_cf/100`, then calls `_control_asim_ctrl(ip, m, t, position, args, k, heat_rate_control, true, gram_atmosphere; cnf=cnf)`, capturing the trajectory array `y` and the resulting `time_switch` pair. The accumulated heat load is read from the final element `y[end,end]`, and the return value is `Q - m.aerodynamics.heat_load_limit`. The `time_switch` produced alongside is published through a `global time_switch, k` declaration so the enclosing function can recover the switch times belonging to the accepted root; the residual itself returns only the scalar the solver needs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `k_cf` | Any | n/a | yes | Positional argument `k_cf`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `func`. Returns `Q - m.aerodynamics.heat_load_limit`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.switch_window_solver_switch_calculation|switch_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:87-87`
- [[gncx.second_switch_solver_second_time_switch_recalc|second_time_switch_recalc]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:69-69`
- [[gncy.switch_window_solver_switch_calculation_with_integration|switch_calculation_with_integration]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:6-6`

**Downstream**

- `callees` → [[gnc.guidance_hooks__control_asim_ctrl|_control_asim_ctrl]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:13-13`
<!-- vulcan:connections:end -->

## Limitations
Using `global time_switch, k` writes into module-level bindings, so two satellites solving concurrently overwrite each other's switch times and the value left behind after `find_zero` belongs to whatever the last evaluation happened to be, not necessarily the accepted root. Each call performs a complete trajectory integration, making the bisection expensive. The `y[end,end]` indexing assumes heat load is the last stored column, an invariant nothing in this file checks.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl` line 6.
