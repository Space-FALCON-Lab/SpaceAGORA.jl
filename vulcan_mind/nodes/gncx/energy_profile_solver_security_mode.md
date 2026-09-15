---
id: gncx.energy_profile_solver_security_mode
label: security_mode
kind: function
source:
  file: src/gnc/guidance/aerobraking/e_edg/energy_profile_solver.jl
  symbol: security_mode
  lines:
  - 1
  - 28
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: GuidanceHooks namespace supplying the security-mode feasibility check
    with the input parameters, mission, drag-passage indicator, argument dictionary,
    and current time.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: switch_window
  type: Vector{Float64}
  units: s
  description: Either the unchanged switch pair or a saturated window that commands
    zero angle of attack for the remainder of the passage.
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
# security_mode

## Purpose
`security_mode` answers a single question: even if the vehicle flies the rest of this passage at the lowest possible heat rate, will the accumulated heat load still exceed the mission limit? If the answer is yes, the guidance abandons energy depletion for the remainder of the passage and commands the minimum-drag attitude.

## Model & Assumptions
The check is built on the best case rather than the current plan. A closed-form passage is propagated from the current position, and a zero angle-of-attack profile is constructed for its entire length, since zero angle is the minimum-heating orientation for the panel geometry. The free-molecular heat rate is evaluated along that profile through `heat_rate_calc`, with the thermal accommodation factor scaled by `args[:multiplicative_factor_heatload]` so the check inherits whatever conservatism the run is configured with. Temperature is held at the planet value and the molecular speed ratio is formed from the closed-form velocity and the thermal speed.

## Design & Implementation
Only the portion of the profile still ahead of the vehicle contributes: the mask `t_cf .> T` zeroes past samples before the heat rate is summed. The integral is approximated as the sum of the rate samples times the uniform trajectory step taken from the first two time entries, and the already-accumulated `heat_load_past` is added to it. When that total exceeds the mission heat-load limit, the routine latches `security_mode` on the configuration state and returns the pair zero and the final closed-form time plus ten thousand seconds, a window that begins immediately and extends far past the end of the passage so the commanded angle stays at zero throughout. Otherwise the existing switch pair is returned unchanged and nothing is latched.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | GuidanceHooks namespace supplying the security-mode feasibility check with the input parameters, mission, drag-passage indicator, argument dictionary, and current time. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `switch_window` | Vector{Float64} | s | — | Either the unchanged switch pair or a saturated window that commands zero angle of attack for the remainder of the passage. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.e_edg_strategy_compute_e_edg_guidance_window_bang|compute_e_edg_guidance_window!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:56-56`

**Downstream**

- `callees` → [[environment.density_models_density_polyfit|density_polyfit]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/energy_profile_solver.jl:11-11`
- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/energy_profile_solver.jl:2-2`
- `callees` → [[gnc.tracking_executor_heat_rate_calc|heat_rate_calc]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/energy_profile_solver.jl:15-15`
- `callees` → [[gncx.closed_form_solution_closed_form|closed_form]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/energy_profile_solver.jl:5-5`
<!-- vulcan:connections:end -->

## Limitations
The rectangular sum assumes a uniform time grid, which the closed-form propagator provides but which is not verified. The check is one-sided and latching, so it can only ever be more conservative and never recovers within a passage. Density comes from a polynomial fit rather than a sampled atmosphere, so an atmosphere significantly denser than the fit can defeat the very margin the check is meant to protect.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/e_edg/energy_profile_solver.jl:1-28`.
