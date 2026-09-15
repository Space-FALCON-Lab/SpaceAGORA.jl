---
id: gnc.constraint_tracking_time_switch_func_affect_bang
label: time_switch_func_affect!
kind: function
source:
  file: src/gnc/control/aerobraking/constraint_tracking.jl
  symbol: time_switch_func_affect!
  lines:
  - 302
  - 302
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
Records the time of each detected bang-bang control switch so the caller can report and plot the switching structure of the optimal panel-angle profile.

## Design & Implementation
On each root of `time_switch_func_condition` it calls `append!(cnf_state.t_time_switch_func, integrator.t)` and returns `nothing` without modifying the integrator state, so the pass continues uninterrupted. `cnf_state` is captured from the enclosing scope. The commented-out branch shows the alternative of terminating once two switches have been collected, which is left disabled so the full trajectory is produced.

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

- `callees` → [[gnc.bridge_helpers__make_aerobraking_runtime_context|_make_aerobraking_runtime_context]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:353-353`
- `callees` → [[gnc.bridge_helpers__with_control_gain|_with_control_gain]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:352-352`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:462-462`
<!-- vulcan:connections:end -->

## Limitations
The vector grows monotonically and is never cleared, so the repeated forward and reverse solves of the shooting loop and the two final refinement solves all append into the same list, interleaving times from opposite integration directions with no marker of which solve produced them. Since it mutates captured state during a callback, concurrent passes sharing a configuration object corrupt the list. Duplicate or near-duplicate entries appear when the residual grazes zero.

## Provenance
Mapped from `src/gnc/control/aerobraking/constraint_tracking.jl` line 302.
