---
id: gnc.targeting_control_evaluate_candidate
label: evaluate_candidate
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: evaluate_candidate
  lines:
  - 973
  - 973
inputs:
- id: t_switch
  type: Any
  units: n/a
  required: true
  description: Positional argument `t_switch`.
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
  description: Return value of `evaluate_candidate`. Returns `_edg_targeting_outcome_with_heat_load(         config,         p,         spacec`.
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

# evaluate_candidate

## Purpose
The closure that evaluates one candidate switch time to its full outcome with heat load, used by certification and both residual solvers.

## Design & Implementation
Captures the config, parameters, spacecraft, state, time and accumulated heat load, and calls `_edg_targeting_outcome_with_heat_load` with the candidate switch time. Defined inside `_edg_solve_targeting_switch`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t_switch` | Any | n/a | yes | Positional argument `t_switch`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `evaluate_candidate`. Returns `_edg_targeting_outcome_with_heat_load(         config,         p,         spacec`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_certify_targeting_candidates|_edg_certify_targeting_candidates]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:898-898`
- [[gnc.targeting_control_apoapsis_residual|apoapsis_residual]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:1033-1033`
- [[gnc.targeting_control_energy_residual|energy_residual]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:1028-1028`
- [[gnc.targeting_control_solve_apoapsis_switch|solve_apoapsis_switch]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:1065-1065`

**Downstream**

- `callees` → [[gnc.targeting_control__edg_certify_targeting_candidates|_edg_certify_targeting_candidates]] · `callers` · feedback · `src/gnc/control/targeting_control.jl:986-986`
- `callees` → [[gnc.targeting_control__edg_disable_uncertified_targeting_bang|_edg_disable_uncertified_targeting!]] · `callers` · call · `src/gnc/control/targeting_control.jl:994-994`
- `callees` → [[gnc.targeting_control__edg_targeting_outcome_with_heat_load|_edg_targeting_outcome_with_heat_load]] · `callers` · call · `src/gnc/control/targeting_control.jl:973-973`
<!-- vulcan:connections:end -->

## Limitations
No memoisation: the certification sweep, the Brent solve and the apoapsis refinement each call it afresh, so the same switch time may be predicted more than once.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 973.
