---
id: gnc.targeting_control__edg_targeting_outcome_with_heat_load
label: _edg_targeting_outcome_with_heat_load
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_targeting_outcome_with_heat_load
  lines:
  - 848
  - 848
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
- id: switch_time_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `switch_time_s`.
- id: accumulated_heat_load_j_cm2
  type: Float64
  units: n/a
  required: true
  description: Positional argument `accumulated_heat_load_j_cm2`.
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
  description: Return value of `_edg_targeting_outcome_with_heat_load`. Returns `merge(outcome,
    (heat_load_j_cm2=accumulated_heat_load_j_cm2 + future_heat_load,)`.
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

# _edg_targeting_outcome_with_heat_load

## Purpose
Extends a targeting prediction with the heat load the predicted profile would accumulate, so candidates can be certified against the heat-load limit.

## Design & Implementation
Runs `_edg_predict_targeting_outcome`, integrates the heat rate along the returned track and angle profile with `_edg_profile_heat_load` (heat-rate control off, since the profile already reflects it), and merges `heat_load_j_cm2` as the accumulated load plus the future load into the outcome tuple.

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
| in | `switch_time_s` | Float64 | n/a | yes | Positional argument `switch_time_s`. |
| in | `accumulated_heat_load_j_cm2` | Float64 | n/a | yes | Positional argument `accumulated_heat_load_j_cm2`. |
| in | `heat_rate_control` | Bool | n/a | yes | Keyword argument `heat_rate_control`. |
| in | `structural_control` | Bool | n/a | yes | Keyword argument `structural_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_targeting_outcome_with_heat_load`. Returns `merge(outcome, (heat_load_j_cm2=accumulated_heat_load_j_cm2 + future_heat_load,)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control_evaluate_candidate|evaluate_candidate]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:973-973`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_profile_heat_load|_edg_profile_heat_load]] · `callers` · call · `src/gnc/control/targeting_control.jl:873-873`
- `callees` → [[gnc.targeting_control__edg_predict_targeting_outcome|_edg_predict_targeting_outcome]] · `callers` · call · `src/gnc/control/targeting_control.jl:861-861`
<!-- vulcan:connections:end -->

## Limitations
The heat-load integration re-evaluates the Maxwellian heat rate at every track sample, adding a third pass over the prediction.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 848.
