---
id: gnc.targeting_control__edg_predict_targeting_outcome
label: _edg_predict_targeting_outcome
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_predict_targeting_outcome
  lines:
  - 670
  - 670
inputs:
- id: config
  type: AerobrakingEnergyDepletionConfig
  units: n/a
  required: true
  description: Positional argument `config`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: pos
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos`.
- id: vel
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel`.
- id: mass_state
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass_state`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: switch_time_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `switch_time_s`.
- id: heat_rate_control
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `heat_rate_control`.
- id: structural_control
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `structural_control`.
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
  description: Return value of `_edg_predict_targeting_outcome`. Returns `(`.
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

# _edg_predict_targeting_outcome

## Purpose
Predicts the end-of-pass orbit for one candidate switch time: energy, periapsis and apoapsis, together with the full track.

## Design & Implementation
Estimates the mass through `_edg_predict_mass`, the passage duration through `_edg_drag_passage_duration`, builds the time grid, integrates the targeting trajectory, and evaluates orbit metrics at the final state. Returns a named tuple carrying the switch time, duration, energy, apsides, track and angle profile.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | AerobrakingEnergyDepletionConfig | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `mass_state` | Float64 | n/a | yes | Positional argument `mass_state`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `switch_time_s` | Float64 | n/a | yes | Positional argument `switch_time_s`. |
| in | `heat_rate_control` | Bool | n/a | yes | Keyword argument `heat_rate_control`. |
| in | `structural_control` | Bool | n/a | yes | Keyword argument `structural_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_predict_targeting_outcome`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_targeting_bracket_outcomes|_edg_targeting_bracket_outcomes]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:819-819`
- [[gnc.targeting_control__edg_targeting_outcome_with_heat_load|_edg_targeting_outcome_with_heat_load]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:861-861`
- [[gnc.targeting_control__edg_targeting_switch_outcomes|_edg_targeting_switch_outcomes]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:780-780`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_drag_passage_duration|_edg_drag_passage_duration]] · `callers` · call · `src/gnc/control/targeting_control.jl:683-683`
- `callees` → [[gnc.heat_load_control__edg_predict_mass|_edg_predict_mass]] · `callers` · call · `src/gnc/control/targeting_control.jl:682-682`
- `callees` → [[gnc.targeting_control__edg_integrated_targeting_trajectory|_edg_integrated_targeting_trajectory]] · `callers` · call · `src/gnc/control/targeting_control.jl:685-685`
- `callees` → [[gnc.targeting_control__edg_orbit_metrics_from_rv|_edg_orbit_metrics_from_rv]] · `callers` · call · `src/gnc/control/targeting_control.jl:699-699`
- `callees` → [[gnc.targeting_control__edg_targeting_prediction_time_grid|_edg_targeting_prediction_time_grid]] · `callers` · call · `src/gnc/control/targeting_control.jl:684-684`
<!-- vulcan:connections:end -->

## Limitations
The duration estimate is made once from the initial state, so if the constrained profile changes how long the vehicle spends below the interface, the grid may end before exit and the final metrics are taken mid-passage.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 670.
