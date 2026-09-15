---
id: gnc.second_switch_solver_func
label: func
kind: function
source:
  file: src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl
  symbol: func
  lines:
  - 7
  - 7
inputs:
- id: t_s
  type: Any
  units: n/a
  required: true
  description: Positional argument `t_s`.
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
Inner residual closure defined inside `second_time_switch_recalc_with_integration`. Given a candidate second bank-reversal time `t_s` in seconds, it returns the signed margin between the predicted end-of-pass integrated heat load and the vehicle's allowable `m.aerodynamics.heat_load_limit`, so that a Brent root finder can drive that margin to zero.

## Theory & Math
Residual $f(t_s) = Q(t_s) - Q_{\max}$, where $Q(t_s)$ is the integrated heat load $\int \dot q\,dt$ over the pass when the second angle-of-attack switch occurs at time $t_s$, and $Q_{\max} = $ `m.aerodynamics.heat_load_limit`. The outer solver seeks $t_s^\star$ with $f(t_s^\star) = 0$.

## Design & Implementation
The closure calls `_control_asim_ctrl(ip, m, t, current_position, args, 0, heat_rate_control, false, gram_atmosphere, t_s, reevaluation_mode; cnf=cnf_state)`, which runs a full forward integration of the atmospheric pass with the switch time set to `t_s`. It takes `Q = y[end, end]`, the last column of the last row of the returned state history, as the accumulated heat load, and returns `Q - m.aerodynamics.heat_load_limit`. It captures `ip`, `m`, `t`, `current_position`, `args`, `heat_rate_control`, `gram_atmosphere`, `reevaluation_mode` and `cnf_state` from the enclosing scope, so each evaluation is a complete trajectory propagation rather than an algebraic expression.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t_s` | Any | n/a | yes | Positional argument `t_s`. |
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

- `callees` → [[gnc.guidance_hooks__control_asim_ctrl|_control_asim_ctrl]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:9-9`
<!-- vulcan:connections:end -->

## Limitations
Every evaluation costs one full numerical propagation, so root finding here is expensive and the enclosing `fzero` call is wrapped in a bare `try`/`catch` that silently discards any failure. The closure assumes `y` is a matrix whose final element is the cumulative heat load; a propagation that returns a shorter or differently shaped history will index incorrectly rather than error clearly. It also assumes `Q` is monotone enough in `t_s` for a bracketed Brent search on `[t, b]` to have a sign change.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl` line 7.
