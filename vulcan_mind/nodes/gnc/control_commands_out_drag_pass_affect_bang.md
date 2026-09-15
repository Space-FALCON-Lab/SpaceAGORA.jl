---
id: gnc.control_commands_out_drag_pass_affect_bang
label: out_drag_pass_affect!
kind: function
source:
  file: src/gnc/control/aerobraking/control_commands.jl
  symbol: out_drag_pass_affect!
  lines:
  - 298
  - 298
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
Callback action executed when the vehicle crosses the atmospheric exit interface. It records the exit epoch into the shared per-satellite state so downstream guidance knows when the drag pass ended, and stops the integration immediately.

## Design & Implementation
`out_drag_pass_affect!(integrator)` performs two mutations: it assigns `cnf_state.t_out_drag_passage = integrator.t`, writing the crossing time in seconds into the closure-captured configuration state obtained earlier from `_bridge_get_cnf(args; cnf=cnf)`, then calls `terminate!(integrator)` so `solve` returns with the pass truncated at the interface. It is registered as the downcrossing action of `ContinuousCallback(out_drag_pass_condition, out_drag_pass_affect!, nothing)`.

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
The captured `cnf_state` is shared mutable state: if the same configuration object backs more than one concurrent integration, the last writer wins and `t_out_drag_passage` no longer corresponds to the pass the caller is reading. `terminate!` ends the solve unconditionally, so an early spurious crossing — for instance from a trajectory that grazes the interface without entering the atmosphere — truncates the pass with no diagnostic. Nothing records that the callback fired, so a caller cannot distinguish a terminated solve from one that simply reached its final time.

## Provenance
Mapped from `src/gnc/control/aerobraking/control_commands.jl` line 298.
