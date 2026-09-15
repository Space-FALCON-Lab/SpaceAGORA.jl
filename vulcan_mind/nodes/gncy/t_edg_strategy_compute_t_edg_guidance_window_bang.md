---
id: gncy.t_edg_strategy_compute_t_edg_guidance_window_bang
label: compute_t_edg_guidance_window!
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/t_edg_strategy.jl
  symbol: compute_t_edg_guidance_window!
  lines:
  - 1
  - 13
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: guidance_input
  type: AerobrakingGuidanceInput
  units: n/a
  required: true
  description: Aerobraking guidance input record supplying the argument dictionary
    and the configuration handle for the current pass.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: guidance_output
  type: AerobrakingGuidanceOutput
  units: n/a
  description: Switch-time pair and security-mode flag carried forward from the persistent
    configuration state.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncy
origin: agent
---

# compute_t_edg_guidance_window!

## Purpose
`compute_t_edg_guidance_window!` is the runtime control-hook entry point for the targeted energy-depletion guidance strategy. It produces the switch window that the aerobraking controller consumes during a pass, and it is the method that `compute_aerobraking_guidance` dispatches to for `TEdgStrategy`.

## Model & Assumptions
The T-EDG targeting problem itself is solved offline by the targeting solver and equations-of-motion predictor paths rather than inside this hook. This function therefore treats the persistent configuration record as the authoritative source of the current switch schedule, and it holds that schedule steady unless mission policy has supplied a dedicated targeting schedule. The assumption is that the schedule was already computed and stored before the controller asks for it.

## Design & Implementation
The body resolves the shared configuration state through `_bridge_get_cnf`, passing both the argument dictionary from the input record and the explicit configuration handle so an override can be threaded through without touching the dictionary. It then constructs an `AerobrakingGuidanceOutput` from `cnf_state.time_switch_1`, `cnf_state.time_switch_2`, and `cnf_state.security_mode`. The companion method on line 15 gives `TEdgStrategy` its dispatch entry and simply forwards the input record here, keeping the strategy tag layer free of logic.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `guidance_input` | AerobrakingGuidanceInput | n/a | yes | Aerobraking guidance input record supplying the argument dictionary and the configuration handle for the current pass. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `guidance_output` | AerobrakingGuidanceOutput | n/a | — | Switch-time pair and security-mode flag carried forward from the persistent configuration state. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.t_edg_strategy_compute_aerobraking_guidance|compute_aerobraking_guidance]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/t_edg_strategy.jl:16-16`

**Downstream**

- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/t_edg_strategy.jl:3-3`
- `callees` → [[gnc.interfaces_aerobrakingguidanceoutput|AerobrakingGuidanceOutput]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/t_edg_strategy.jl:8-8`
<!-- vulcan:connections:end -->

## Limitations
Despite the exclamation-mark naming convention the function does not currently mutate the configuration state; it only reads three fields from it. Because the schedule is passed through unchanged, a stale or never-populated configuration silently yields a zero switch window rather than an error. The security-mode flag is forwarded without being re-evaluated against the current pass conditions, so a security condition raised earlier persists until something else clears it.

## Provenance
Mapped from t_edg_strategy.jl lines 1-17; include site observed at guidance_hooks.jl line 90.
