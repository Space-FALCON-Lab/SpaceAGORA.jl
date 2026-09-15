---
id: gnc.trajectory_predictor_out_drag_passage_condition
label: out_drag_passage_condition
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl
  symbol: out_drag_passage_condition
  lines:
  - 378
  - 378
inputs:
- id: y
  type: Any
  units: n/a
  required: true
  description: Positional argument `y`.
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
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
  description: Return value of `out_drag_passage_condition`. Returns `norm(y[1:3])
    * cnf_state.DU - mission.planet.Rp_e - settings.exit_interface_m`.
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

# out_drag_passage_condition

## Purpose
Exit-interface root function for the T-EDG targeting propagation, expressed in canonical length units. Besides locating the crossing for the terminating callback, it latches the post-passage angle-of-attack command so the vehicle leaves the atmosphere in the attitude the selected heat-load solution requires.

## Design & Implementation
`out_drag_passage_condition(y, t, integrator)` reads `mission` and `settings` from `integrator.p` and returns `norm(y[1:3]) * cnf_state.DU - mission.planet.Rp_e - settings.exit_interface_m`, converting the non-dimensional position back to metres with the captured distance unit `DU`. Before returning it checks whether the same expression is within `1e-5` m of zero, and if so writes `cnf_state.α`: `mission.aerodynamics.α` when `settings.heat_load_solution` is `0` or `2`, and `0.0` when it is `1` or `3`. The function is registered as the condition of `ContinuousCallback(out_drag_passage_condition, out_drag_passage_affect!, nothing)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `y` | Any | n/a | yes | Positional argument `y`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `out_drag_passage_condition`. Returns `norm(y[1:3]) * cnf_state.DU - mission.planet.Rp_e - settings.exit_interface_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:378-378`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Performing a state mutation inside a root-finding condition is unsound: the callback machinery evaluates the condition at trial points, during rejected steps and during bracket refinement, so `cnf_state.α` can be latched from a probe that is never accepted. The `1e-5` m tolerance on a radius of several million metres is far below what the root solver needs to hit, so in practice the latch may never fire at all. `Rp_e` is the equatorial radius, making the interface spherical rather than altitude-constant, and `heat_load_solution` values outside `0` through `3` leave the attitude untouched with no warning.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl` line 378.
