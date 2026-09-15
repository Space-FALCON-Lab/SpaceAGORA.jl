---
id: gnc.control_commands_time_switch_func_affect_bang
label: time_switch_func_affect!
kind: function
source:
  file: src/gnc/control/aerobraking/control_commands.jl
  symbol: time_switch_func_affect!
  lines:
  - 320
  - 320
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
  description: Return value of `time_switch_func_affect!`; mutates `integrator` in
    place.
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

# time_switch_func_affect!

## Purpose
Callback action recording each detected switching-surface crossing. It appends the crossing epoch to a growing list in the shared configuration state, from which `asim_ctrl` later reads the first and last entries to form the pair of angle-of-attack switch times returned to the guidance law.

## Design & Implementation
`time_switch_func_affect!(integrator)` calls `append!(cnf_state.t_time_switch_func, integrator.t)`, pushing the crossing time in seconds onto the captured vector, and returns `nothing` so the integration continues rather than terminating. It is the sole action of `ContinuousCallback(time_switch_func_condition, time_switch_func_affect!)`, registered without a separate upcrossing handler so both crossing directions are recorded. `asim_ctrl` reads `cnf_state.t_time_switch_func` after the solve and resets it to `[]` before the next pass.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `time_switch_func_affect!`; mutates `integrator` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:302-302`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:320-320`

**Downstream**

- `callees` → [[gnc.bridge_helpers__make_aerobraking_runtime_context|_make_aerobraking_runtime_context]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:365-365`
- `callees` → [[gnc.bridge_helpers__with_control_gain|_with_control_gain]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:364-364`
- `callees` → [[gncx.closed_form_solution_closed_form|closed_form]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:331-331`
<!-- vulcan:connections:end -->

## Limitations
Appending unconditionally means every crossing is recorded, including chattering near a shallow tangency of the switching surface, so the vector can hold more than the two entries the caller expects and gives no way to tell a genuine switch from numerical noise. The reset to `[]` happens in the enclosing function, so a solve that throws leaves stale entries behind that contaminate the next pass. The vector lives in shared `cnf_state`, so concurrent integrations interleave their crossings into one list with no owner tag.

## Provenance
Mapped from `src/gnc/control/aerobraking/control_commands.jl` line 320.
