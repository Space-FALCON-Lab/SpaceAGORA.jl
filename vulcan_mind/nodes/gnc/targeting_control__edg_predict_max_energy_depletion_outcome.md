---
id: gnc.targeting_control__edg_predict_max_energy_depletion_outcome
label: _edg_predict_max_energy_depletion_outcome
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_predict_max_energy_depletion_outcome
  lines:
  - 711
  - 711
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
- id: env
  type: Any
  units: n/a
  required: true
  description: Positional argument `env`.
- id: heat_load_j_cm2
  type: Float64
  units: n/a
  required: true
  description: Positional argument `heat_load_j_cm2`.
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
  description: Return value of `_edg_predict_max_energy_depletion_outcome`. Returns
    `(`.
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

# _edg_predict_max_energy_depletion_outcome

## Purpose
Predicts the end-of-pass orbit under the max-energy-depletion profile, solving the heat-load window first if that sub-mode is active.

## Design & Implementation
Estimates mass and duration, builds the grid, solves `_edg_solve_heat_load_switches` when the heat-load sub-mode is on and otherwise uses an infinite window, integrates the max-energy trajectory and evaluates final orbit metrics. Returns the same tuple shape as the targeting prediction with `switch_time_s` set to `NaN` and the window added.

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
| in | `env` | Any | n/a | yes | Positional argument `env`. |
| in | `heat_load_j_cm2` | Float64 | n/a | yes | Positional argument `heat_load_j_cm2`. |
| in | `heat_rate_control` | Bool | n/a | yes | Keyword argument `heat_rate_control`. |
| in | `structural_control` | Bool | n/a | yes | Keyword argument `structural_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_predict_max_energy_depletion_outcome`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_targeting_bracket_outcomes|_edg_targeting_bracket_outcomes]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:832-832`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_drag_passage_duration|_edg_drag_passage_duration]] · `callers` · call · `src/gnc/control/targeting_control.jl:725-725`
- `callees` → [[gnc.heat_load_control__edg_predict_mass|_edg_predict_mass]] · `callers` · call · `src/gnc/control/targeting_control.jl:724-724`
- `callees` → [[gnc.targeting_control__edg_integrated_max_energy_depletion_trajectory|_edg_integrated_max_energy_depletion_trajectory]] · `callers` · call · `src/gnc/control/targeting_control.jl:742-742`
- `callees` → [[gnc.targeting_control__edg_orbit_metrics_from_rv|_edg_orbit_metrics_from_rv]] · `callers` · call · `src/gnc/control/targeting_control.jl:756-756`
- `callees` → [[gnc.targeting_control__edg_targeting_prediction_time_grid|_edg_targeting_prediction_time_grid]] · `callers` · call · `src/gnc/control/targeting_control.jl:726-726`
- `callees` → [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callers` · call · `src/gnc/control/targeting_control.jl:728-728`
<!-- vulcan:connections:end -->

## Limitations
The heat-load solve is itself a prediction, so this function nests two passes of integration and can be the single most expensive call in the controller.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 711.
