---
id: gncx.struct_load_control__edg_structural_alpha
label: _edg_structural_alpha
kind: function
source:
  file: src/gnc/control/struct_load_control.jl
  symbol: _edg_structural_alpha
  lines:
  - 97
  - 107
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace supplying the structural-load constraint projector
    with the energy-depletion configuration, ODE parameters, environment, spacecraft,
    controlled panel links, and baseline angle.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: alpha_limited
  type: Float64
  units: rad
  description: Baseline angle of attack reduced as needed so the predicted panel drag
    load stays within the structural limit.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncx
origin: agent
---
# _edg_structural_alpha

## Purpose
`_edg_structural_alpha` projects a desired angle of attack onto the structurally feasible set. It is the third constraint in the energy-depletion stack, alongside the instantaneous heat-rate limit and the integrated heat-load limit, and it protects the deployed solar panels from aerodynamic loading during a drag passage.

## Model & Assumptions
The routine is gated on flight phase: `_edg_in_drag_passage` decides whether the vehicle is inside the atmospheric segment at all, and outside that segment the baseline angle is returned unchanged so the constraint costs nothing during the coast arc. Inside the passage the work is delegated to `_energy_depletion_struct_load_root_alpha`, which searches for the angle at which the predicted drag load on the controlled panel links equals the structural limit.

## Design & Implementation
The underlying load model is built by `_energy_depletion_struct_drag_area`, which sums the projected drag-producing area contributed by the controlled panel links at a candidate angle, so the constraint depends on which links are declared controllable rather than on the whole vehicle. The residual is that predicted load minus the limit, and the solve uses `Roots.find_zero` with an explicit `Roots.Bisection()` over the bracket from the configured minimum angle to the baseline angle. Bisection is chosen over a derivative method here because the projected-area model is piecewise in the link geometry and not reliably smooth. A `catch` around the solve falls back to the minimum angle, the safest end of the bracket, and the returned value is clamped into the bracket regardless of which path produced it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace supplying the structural-load constraint projector with the energy-depletion configuration, ODE parameters, environment, spacecraft, controlled panel links, and baseline angle. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `alpha_limited` | Float64 | rad | — | Baseline angle of attack reduced as needed so the predicted panel drag load stays within the structural limit. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_command_alpha_bang|_edg_command_alpha!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:199-199`

**Downstream**

- `callees` → [[gnc.struct_load_control__energy_depletion_struct_load_root_alpha|_energy_depletion_struct_load_root_alpha]] · `callers` · call · `src/gnc/control/struct_load_control.jl:106-106`
- `callees` → [[gnc.targeting_control__edg_in_drag_passage|_edg_in_drag_passage]] · `callers` · call · `src/gnc/control/struct_load_control.jl:105-105`
<!-- vulcan:connections:end -->

## Limitations
The load model accounts only for the aerodynamic drag on the controlled panels and ignores bending moments, hinge and boom loads, thermal-structural coupling, and any dynamic amplification. Bisection needs the residual to change sign across the bracket; when it does not, the fallback to the minimum angle is conservative but discards the intended control authority. The constraint is instantaneous, so it does not bound cumulative fatigue over repeated passages.

## Provenance
Mapped from `src/gnc/control/struct_load_control.jl:97-107`, with the root solve at lines 27-95.
