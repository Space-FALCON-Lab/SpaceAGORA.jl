---
id: gnc.trajectory_predictor_out_drag_passage_affect_bang
label: out_drag_passage_affect!
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl
  symbol: out_drag_passage_affect!
  lines:
  - 392
  - 392
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
  description: Return value of `out_drag_passage_affect!`; mutates `integrator` in
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

# out_drag_passage_affect!

## Purpose
Terminating action of the T-EDG targeting propagation. When the vehicle reaches the atmospheric exit interface, this stops the integration so the shooting method evaluates its targeting residual at the end of the drag passage rather than at an arbitrary final time.

## Design & Implementation
`out_drag_passage_affect!(integrator)` has a single statement, `terminate!(integrator)`, which sets the solver's return code and unwinds the integration loop. It is paired with `out_drag_passage_condition` in `ContinuousCallback(out_drag_passage_condition, out_drag_passage_affect!, nothing)`; passing `nothing` as the third argument suppresses the upcrossing action, so only the crossing in the decreasing-residual direction terminates. The nominal final time supplied to the `ODEProblem` is `(time_0 + 1e8)/cnf_state.TU`, roughly three years of canonical time, which exists purely as a backstop for a pass that never exits.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `out_drag_passage_affect!`; mutates `integrator` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:392-392`

**Downstream**

- `callees` → [[gnc.bridge_helpers__with_time_switch|_with_time_switch]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:400-400`
<!-- vulcan:connections:end -->

## Limitations
Unlike the corresponding handler in the control-command path, this records nothing: no exit epoch, no state snapshot, so the caller must recover the crossing from `sol.t[end]` and cannot distinguish a genuine termination from the solver failing for another reason. If the interface is never reached the solve grinds toward the `1e8` second horizon at `abstol`/`reltol` of `1e-9`, which is effectively a hang rather than an error. No exception is raised and no flag is set, so an upstream targeting iteration can silently consume a meaningless trajectory.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl` line 392.
