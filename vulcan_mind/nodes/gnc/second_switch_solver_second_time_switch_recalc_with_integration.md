---
id: gnc.second_switch_solver_second_time_switch_recalc_with_integration
label: second_time_switch_recalc_with_integration
kind: function
source:
  file: src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl
  symbol: second_time_switch_recalc_with_integration
  lines:
  - 3
  - 3
inputs:
- id: ip
  type: Any
  units: n/a
  required: true
  description: Positional argument `ip`.
- id: m
  type: Any
  units: n/a
  required: true
  description: Positional argument `m`.
- id: position
  type: Any
  units: n/a
  required: true
  description: Positional argument `position`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
- id: heat_rate_control
  type: Any
  units: n/a
  required: true
  description: Positional argument `heat_rate_control`.
- id: reevaluation_mode
  type: Any
  units: n/a
  required: true
  description: Positional argument `reevaluation_mode`.
- id: gram_atmosphere
  type: Any
  units: n/a
  required: false
  description: Positional argument `gram_atmosphere` (default `nothing`).
- id: current_position
  type: Any
  units: n/a
  required: false
  description: Positional argument `current_position` (default `0`).
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
  description: Return value of `second_time_switch_recalc_with_integration`. Returns
    `Q - m.aerodynamics.heat_load_limit` or `cnf_state.time_switch_1, cnf_state.time_switch_2`
    or `cnf_state.time_switch_1, t` or `cnf_state.time_switch_1, time_switch`.
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

# second_time_switch_recalc_with_integration

## Purpose
Recomputes the second angle-of-attack switch time for an E-EDG aerobraking pass by re-integrating the trajectory, so the accumulated heat load lands on `m.aerodynamics.heat_load_limit`. Returns the tuple `(time_switch_1, time_switch_2)` in seconds, leaving the first switch untouched.

## Theory & Math
Solve $f(t_s)=Q(t_s)-Q_{\max}=0$ for $t_s \in [t,\,b]$ by Brent's method, where $Q(t_s)$ is the pass heat load produced by integrating with the second switch at $t_s$, $Q_{\max}$ is the heat-load limit, and $b = t + 200\,\mathrm{s}$ or $t + 1500\,\mathrm{s}$ depending on the selected heat-load solution branch.

## Design & Implementation
It fetches guidance state with `_bridge_get_cnf(args; cnf=cnf)` and seeds `time_switch` from `cnf_state.time_switch_2`. The nested `func` evaluates the heat-load residual by propagation. It samples the residual at the stored switch time and at the current time `t`, then short-circuits: if `abs(delta_Q_switch) <= 0.01` the existing switch is kept (with an extra guard that rejects the keep when the two residuals differ by less than 0.5 while `abs(t - time_switch) > 20` seconds), and if the current-time residual is already negative under `heat_load_sol` 0 or 3 it returns `t` as the new switch. Otherwise it brackets on `[t, b]`, where `b = t + 200` for `heat_load_sol` 0 or 1 and `t + 1500` otherwise, and calls `fzero(ts -> func(ts), [t, b], Roots.Brent())`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ip` | Any | n/a | yes | Positional argument `ip`. |
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `position` | Any | n/a | yes | Positional argument `position`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `heat_rate_control` | Any | n/a | yes | Positional argument `heat_rate_control`. |
| in | `reevaluation_mode` | Any | n/a | yes | Positional argument `reevaluation_mode`. |
| in | `gram_atmosphere` | Any | n/a | no | Positional argument `gram_atmosphere` (default `nothing`). |
| in | `current_position` | Any | n/a | no | Positional argument `current_position` (default `0`). |
| in | `cnf` | Any | n/a | no | Keyword argument `cnf` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `second_time_switch_recalc_with_integration`. Returns `Q - m.aerodynamics.heat_load_limit` or `cnf_state.time_switch_1, cnf_state.time_switch_2` or `cnf_state.time_switch_1, t` or `cnf_state.time_switch_1, time_switch`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.e_edg_strategy_compute_e_edg_guidance_window_bang|compute_e_edg_guidance_window!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:28-28`

**Downstream**

- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl:4-4`
<!-- vulcan:connections:end -->

## Limitations
The `fzero` call sits in a `try`/`catch` that swallows all errors, including a non-bracketing interval, so a failed solve silently returns the unchanged incoming `time_switch_2` with no diagnostic. The local `x_tol` is assigned in both branches of the `reevaluation_mode` test as 0.01 and is then never used. The 0.01, 0.5, 20-second, 200-second and 1500-second thresholds are hard-coded rather than configured. `args[:heat_load_sol]` must be present or the lookup throws a `KeyError`.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl` line 3.
