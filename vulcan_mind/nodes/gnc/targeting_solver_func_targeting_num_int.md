---
id: gnc.targeting_solver_func_targeting_num_int
label: func_targeting_num_int
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl
  symbol: func_targeting_num_int
  lines:
  - 119
  - 119
inputs:
- id: t_switch
  type: Any
  units: n/a
  required: true
  description: Positional argument `t_switch`.
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
  description: Return value of `func_targeting_num_int`. Returns `(energy_fin - energy_f)
    / 1e6`.
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

# func_targeting_num_int

## Purpose
Residual closure defined inside `control_solarpanels_targeting_num_int`: for a trial switch time `t_switch` it integrates the aerobraking pass and returns the scaled difference between the achieved final specific energy and the requested `energy_f`. Brent's method in the enclosing function drives this residual to zero.

## Design & Implementation
Captures `param`, `time_0`, `in_cond`, `energy_f` and `log_enabled` from the enclosing scope. It calls `asim_ctrl_targeting(t_switch, param, time_0, in_cond; cnf=_bridge_get_cnf(param))` to obtain a solution array `sol`, reads `m = param.mission`, and evaluates `energy_fin = norm(sol[4:6,end])^2/2 - m.planet.μ / norm(sol[1:3,end])` in J/kg. When logging is enabled it prints the `(t_switch, energy_fin)` pair. The return value is `(energy_fin - energy_f) / 1e6`, dividing by one million so the root finder works with numbers of order unity.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t_switch` | Any | n/a | yes | Positional argument `t_switch`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `func_targeting_num_int`. Returns `(energy_fin - energy_f) / 1e6`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:128-128`
- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:121-121`
- `callees` → [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:121-121`
<!-- vulcan:connections:end -->

## Limitations
Every evaluation performs a full ODE propagation, so the cost of the outer root search scales with Brent iteration count. It assumes `sol` is indexable with rows 1:3 as position and 4:6 as velocity in SI units. The closure has no guard for a failed or truncated integration (`sol` shorter than expected simply raises a `BoundsError`). The `1e6` scale factor is fixed regardless of the energy magnitude for the planet in use.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl` line 119.
