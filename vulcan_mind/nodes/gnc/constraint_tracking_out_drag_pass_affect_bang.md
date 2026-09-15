---
id: gnc.constraint_tracking_out_drag_pass_affect_bang
label: out_drag_pass_affect!
kind: function
source:
  file: src/gnc/control/aerobraking/constraint_tracking.jl
  symbol: out_drag_pass_affect!
  lines:
  - 278
  - 278
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
Callback effect that records when the vehicle exits the drag passage and stops the integration at that point.

## Design & Implementation
On a root of `out_drag_pass_condition`, it writes `integrator.t` into `cnf_state.t_out_drag_passage` and calls `terminate!(integrator)`. The `cnf_state` binding is captured from the enclosing `asim_ctrl_plot` scope, so the exit time is visible to the caller after `solve` returns, and terminating means the returned solution always ends exactly at the interface crossing rather than at the nominal final time of `time_0 + 1500` seconds.

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
The field is overwritten rather than appended, so with a callback registered across the repeated forward and reverse solves of the shooting loop only the most recent crossing survives, and after a reverse-time solve the stored value is a time from that backward pass. Mutating captured configuration state from inside a solver callback makes the routine non-reentrant: two passes sharing one `cnf_state` object, or running on separate threads, race on the same field.

## Provenance
Mapped from `src/gnc/control/aerobraking/constraint_tracking.jl` line 278.
