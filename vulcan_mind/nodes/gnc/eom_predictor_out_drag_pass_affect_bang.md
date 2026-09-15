---
id: gnc.eom_predictor_out_drag_pass_affect_bang
label: out_drag_pass_affect!
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl
  symbol: out_drag_pass_affect!
  lines:
  - 728
  - 728
inputs:
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
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
  description: Return value of `out_drag_pass_affect!`; mutates `integrator` in place.
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

# out_drag_pass_affect!

## Purpose
Affect function paired with `out_drag_pass_condition` in the `out_drag_pass` `ContinuousCallback`. When the integrated trajectory crosses the exit interface it records the crossing time and stops the integration, so each shooting iteration and the final propagation in `asim_ctrl_targeting_plot` end exactly at atmospheric exit rather than at the 1500 s hard limit.

## Design & Implementation
Receives the DiffEq `integrator`, writes `integrator.t` (seconds since `time_0`) into `cnf_state.t_out_drag_passage` where `cnf_state` is the configuration object captured from the enclosing `asim_ctrl_targeting_plot` scope, then calls `terminate!(integrator)`. Termination makes `sol[:, end]` the state at exit, which is what `shooting_residual!` reads to form its boundary residuals.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `out_drag_pass_affect!`; mutates `integrator` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:278-278`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:298-298`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:728-728`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The function mutates `cnf_state` captured by closure, so every Newton iteration of `nlsolve` overwrites `t_out_drag_passage`; only the last propagation's value survives. It cannot be reused outside the enclosing function. Because it terminates unconditionally, a trajectory that starts above the interface and re-enters would end at the first crossing. There is no guard against the callback never firing, in which case the integration silently runs to `final_time` and `t_out_drag_passage` retains its previous value.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl` line 728.
