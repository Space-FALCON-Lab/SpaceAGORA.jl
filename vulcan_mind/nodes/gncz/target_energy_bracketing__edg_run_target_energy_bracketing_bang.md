---
id: gncz.target_energy_bracketing__edg_run_target_energy_bracketing_bang
label: _edg_run_target_energy_bracketing!
kind: function
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: _edg_run_target_energy_bracketing!
  lines:
  - 261
  - 334
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: GNC guidance namespace exposing the energy depletion guidance model,
    its configuration, and its per-spacecraft state.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: guidance_state
  type: AerobrakingEnergyDepletionState
  units: mixed
  description: Per-spacecraft bracket energies, target energy, reachability flag,
    and the selected guidance mode for the current drag passage.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncz
origin: agent
---

# _edg_run_target_energy_bracketing!

## Purpose
`_edg_run_target_energy_bracketing!` decides, once per drag passage and per spacecraft, whether the configured target apoapsis is reachable by modulating drag between the two extreme flight profiles. It is the deciding step of the aerobraking energy depletion guidance model and it writes the mode that the control layer subsequently executes.

## Theory & Math
Two endpoint profiles are evaluated for the current pass, a safe low-drag profile and a maximum energy depletion profile, each producing an exit specific energy, a periapsis radius, and an apoapsis radius. Periapsis is treated as affine in exit energy between the endpoints, so for a candidate energy $E$ the interpolated periapsis is $r_p(E) = r_p^{min} + \frac{E - E_{min}}{E_{max} - E_{min}}(r_p^{max} - r_p^{min})$. The target energy solves the fixed-point residual $E - E_{des}(r_a^{target}, r_p(E)) = 0$, where $E_{des}$ is the specific energy of an orbit with the requested apoapsis and the interpolated periapsis.

## Model & Assumptions
The bracket is evaluated at most once per passage, guarded by an evaluated flag, and only while the vehicle is inside a drag passage. Reachability requires the solved energy to lie inside the endpoint energy bracket and the requested apoapsis radius to lie inside the endpoint apoapsis bracket, each with a relative tolerance scaled by the magnitudes involved.

## Design & Implementation
State is pulled through accessors that tolerate both component-array spacecraft states and bare vectors. Heat load over the controlled panel links, heat rate control, and structural load control are forwarded to the endpoint evaluation as configuration-driven switches. On success the routine records the bracket bounds, increments a counter, and selects targeting; otherwise it falls back to maximum energy depletion when that mode is configured and to safe low drag when it is not.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | GNC guidance namespace exposing the energy depletion guidance model, its configuration, and its per-spacecraft state. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `guidance_state` | AerobrakingEnergyDepletionState | mixed | — | Per-spacecraft bracket energies, target energy, reachability flag, and the selected guidance mode for the current drag passage. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.target_energy_bracketing_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:203-203`

**Downstream**

- `callees` → [[gnc.guidance_hooks__control_module|_control_module]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:274-274`
- `callees` → [[gnc.target_energy_bracketing__edg_pos_vel_mass|_edg_pos_vel_mass]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:282-282`
- `callees` → [[gnc.target_energy_bracketing__edg_sat_state|_edg_sat_state]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:281-281`
- `callees` → [[gnc.target_energy_bracketing__edg_set_targeting_fallback_bang|_edg_set_targeting_fallback!]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:277-277`
- `callees` → [[gnc.target_energy_bracketing__edg_target_energy_from_reachable_bracket|_edg_target_energy_from_reachable_bracket]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:309-309`
<!-- vulcan:connections:end -->

## Limitations
Linear interpolation of periapsis against energy is only valid when the two endpoints straddle a mildly nonlinear response, and the reachability test inherits that error. The bracket is not re-evaluated later in the same passage, so a mid-pass atmosphere change cannot reopen a decision.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl:191-334`.
