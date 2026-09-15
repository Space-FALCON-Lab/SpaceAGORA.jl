---
id: gnc.targeting_solver_func_targeting_heatload
label: func_targeting_heatload
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl
  symbol: func_targeting_heatload
  lines:
  - 142
  - 142
inputs:
- id: v_E
  type: Any
  units: n/a
  required: true
  description: Positional argument `v_E`.
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
  description: Return value of `func_targeting_heatload`. Returns `(energy_fin - energy_f)
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

# func_targeting_heatload

## Purpose
Residual closure defined inside `control_solarpanels_targeting_heatload`: for a trial value of the control parameter `v_E` it propagates a heat-rate-constrained aerobraking pass and returns the scaled energy miss against `energy_f`. It is the function whose root `find_zero` locates on the bracket `[1, 1000]`.

## Design & Implementation
Captures `param`, `OE`, `energy_f` and `log_enabled`. On each call it unpacks `m = param.mission`, `ip = param.ip`, `time_0 = param.time_0`, `args = param.args` and `gram_atmosphere = param.gram_atmosphere`, then calls `_control_asim_ctrl_rf(ip, m, time_0, OE, args, v_E, 1.0, false, gram_atmosphere; cnf=_bridge_get_cnf(param))`, keeping only the first return value `sol`. The final specific energy is `norm(sol[4:6,end])^2/2 - m.planet.μ / norm(sol[1:3,end])`, optionally printed alongside `v_E`, and the residual `(energy_fin - energy_f) / 1e6` is returned.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `v_E` | Any | n/a | yes | Positional argument `v_E`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `func_targeting_heatload`. Returns `(energy_fin - energy_f) / 1e6`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:154-154`
- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:149-149`
- `callees` → [[gnc.guidance_hooks__control_asim_ctrl_rf|_control_asim_ctrl_rf]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:149-149`
<!-- vulcan:connections:end -->

## Limitations
The literal positional arguments `1.0` and `false` passed to `_control_asim_ctrl_rf` are undocumented at the call site. The closure re-reads all fields from `param` on every call, which is harmless but means a mutated `param` mid-search changes the residual. No exception handling exists for a propagation that terminates early. The residual assumes SI energy units through the `1e6` divisor.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl` line 142.
