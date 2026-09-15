---
id: gncx.e_edg_strategy_compute_e_edg_guidance_window_bang
label: compute_e_edg_guidance_window!
kind: function
source:
  file: src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl
  symbol: compute_e_edg_guidance_window!
  lines:
  - 1
  - 64
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: GuidanceHooks namespace supplying the energy-depletion guidance window
    computation with the aerobraking guidance input record.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: guidance_window
  type: AerobrakingGuidanceOutput
  units: s
  description: First and second angle-of-attack switch times together with the security-mode
    flag.
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
# compute_e_edg_guidance_window!

## Purpose
`compute_e_edg_guidance_window!` is the energy-depletion guidance strategy. It decides, once per guidance call, whether to compute the angle-of-attack switch window for the first time, to re-evaluate an existing window, or to fall back to a security mode that abandons energy depletion in favour of thermal survival.

## Model & Assumptions
The three branches are mutually exclusive and selected by state rather than by argument. On the first call of a passage, when `evaluate_switch_heat_load` is not yet set, the window is computed from scratch; whether it is computed by the closed-form `switch_calculation` or by the integration-based variants is chosen by `args[:flash2_through_integration]` and by the `args[:heat_load_sol]` selector, which routes values zero and one to the switch calculation and values two and three to the second-switch recalculation. The flag is then set so later calls take a different branch. Re-evaluation is heavily rate-limited: it requires either the ascending phase or proximity to the second switch, and then a timer condition that tightens as the switch approaches, allowing ten seconds between updates far out, three seconds within fifty seconds of the switch, and eight tenths of a second within three seconds of it.

## Design & Implementation
The security branch is the safety net. It triggers when the recorded heat load has reached ninety-eight percent of the mission limit, when the load is still rising by more than two units per sample, when security mode is enabled by `args[:security_mode]` and not already latched, and when the drag-passage indicator says the vehicle is inside the passage. It then calls `security_mode`, which can return a window that holds the angle at zero for the remainder of the passage. Configuration state is resolved once through `_bridge_get_cnf` and mutated in place, and enabling integration-based reevaluation clears `args[:security_mode]` at the top of the function. The result is packaged as an `AerobrakingGuidanceOutput`, and `compute_aerobraking_guidance(::EEdgStrategy, input)` at the end of the file is the single-line dispatch target that makes this the E-EDG strategy.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | GuidanceHooks namespace supplying the energy-depletion guidance window computation with the aerobraking guidance input record. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `guidance_window` | AerobrakingGuidanceOutput | s | — | First and second angle-of-attack switch times together with the security-mode flag. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.e_edg_strategy_compute_aerobraking_guidance|compute_aerobraking_guidance]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:67-67`

**Downstream**

- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:12-12`
- `callees` → [[gnc.interfaces_aerobrakingguidanceoutput|AerobrakingGuidanceOutput]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:59-59`
- `callees` → [[gnc.second_switch_solver_second_time_switch_recalc_with_integration|second_time_switch_recalc_with_integration]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:28-28`
- `callees` → [[gnc.switch_window_solver_switch_calculation|switch_calculation]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:32-32`
- `callees` → [[gncx.energy_profile_solver_security_mode|security_mode]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:56-56`
- `callees` → [[gncx.second_switch_solver_second_time_switch_recalc|second_time_switch_recalc]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:49-49`
- `callees` → [[gncy.switch_window_solver_switch_calculation_with_integration|switch_calculation_with_integration]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:24-24`
<!-- vulcan:connections:end -->

## Limitations
The branch conditions mix several thresholds written as literals rather than configuration, so tuning the re-evaluation cadence requires editing the source. Because the security flag latches in the configuration state, a passage that enters security mode cannot return to nominal energy depletion within that passage. The function mutates its input configuration in place, so calling it speculatively to preview a window is not safe.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:1-64`.
