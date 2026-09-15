---
id: gnc.targeting_control__edg_targeting_bracket_outcomes
label: _edg_targeting_bracket_outcomes
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_targeting_bracket_outcomes
  lines:
  - 807
  - 807
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
- id: mass
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: heat_load_j_cm2
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `heat_load_j_cm2` (default `0.0`).
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
  description: Return value of `_edg_targeting_bracket_outcomes`. Returns `low_drag,
    max_energy_depletion`.
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

# _edg_targeting_bracket_outcomes

## Purpose
Brackets the achievable end-of-pass energy between the immediate-switch outcome and the max-energy-depletion outcome, respecting the heat-load window.

## Design & Implementation
Predicts the low-drag outcome with the switch at `t`, samples the environment, and predicts the max-energy-depletion outcome with the accumulated heat load. Returns the pair.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | AerobrakingEnergyDepletionConfig | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `heat_load_j_cm2` | Float64 | n/a | no | Keyword argument `heat_load_j_cm2` (default `0.0`). |
| in | `heat_rate_control` | Bool | n/a | yes | Keyword argument `heat_rate_control`. |
| in | `structural_control` | Bool | n/a | yes | Keyword argument `structural_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_targeting_bracket_outcomes`. Returns `low_drag, max_energy_depletion`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.targeting_control__edg_predict_max_energy_depletion_outcome|_edg_predict_max_energy_depletion_outcome]] · `callers` · call · `src/gnc/control/targeting_control.jl:832-832`
- `callees` → [[gnc.targeting_control__edg_predict_targeting_outcome|_edg_predict_targeting_outcome]] · `callers` · call · `src/gnc/control/targeting_control.jl:819-819`
- `callees` → [[gnc.targeting_control__edg_targeting_prediction_environment|_edg_targeting_prediction_environment]] · `callers` · call · `src/gnc/control/targeting_control.jl:831-831`
<!-- vulcan:connections:end -->

## Limitations
Shares the double-prediction cost of the switch-outcomes variant; because the max-energy profile can enter the heat-load window, its energy is not guaranteed to be the true minimum achievable.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 807.
