---
id: gnc.targeting_control__edg_certify_targeting_candidates
label: _edg_certify_targeting_candidates
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_certify_targeting_candidates
  lines:
  - 883
  - 883
inputs:
- id: candidate_times
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `candidate_times`.
- id: evaluate_candidate
  type: Any
  units: n/a
  required: true
  description: Positional argument `evaluate_candidate`.
- id: heat_load_limit_j_cm2
  type: Real
  units: n/a
  required: true
  description: Keyword argument `heat_load_limit_j_cm2`.
- id: energy_order_tolerance_jkg
  type: Real
  units: n/a
  required: true
  description: Keyword argument `energy_order_tolerance_jkg`.
- id: heat_load_tolerance_j_cm2
  type: Real
  units: n/a
  required: true
  description: Keyword argument `heat_load_tolerance_j_cm2`.
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
  description: Return value of `_edg_certify_targeting_candidates`. Returns `(times=certified_times,
    outcomes=certified_outcomes, failure=:nonfinite_energy, ` or `(times=certified_times,
    outcomes=certified_outcomes, failure=:heat_load, failure` or `(times=certified_times,
    outcomes=certified_outcomes, failure=:energy_order, fail` or `(times=certified_times,
    outcomes=certified_outcomes, failure=:none, failure_time`.
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

# _edg_certify_targeting_candidates

## Purpose
Walks the candidate switch times in order and keeps only the prefix along which predicted energy decreases monotonically and heat load stays within limit, so the root solve is confined to a well-behaved bracket.

## Design & Implementation
Validates a non-empty candidate list and non-negative tolerances. It evaluates the first candidate and fails with `:nonfinite_energy` or `:heat_load` if it is unusable. For each later candidate it requires finite energy, heat load at or below the limit plus tolerance, and energy strictly below the previous by more than `energy_order_tolerance_jkg`, failing with `:energy_order` otherwise. Returns the certified times, outcomes, the failure kind and the failing time.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `candidate_times` | AbstractVector{<:Real} | n/a | yes | Positional argument `candidate_times`. |
| in | `evaluate_candidate` | Any | n/a | yes | Positional argument `evaluate_candidate`. |
| in | `heat_load_limit_j_cm2` | Real | n/a | yes | Keyword argument `heat_load_limit_j_cm2`. |
| in | `energy_order_tolerance_jkg` | Real | n/a | yes | Keyword argument `energy_order_tolerance_jkg`. |
| in | `heat_load_tolerance_j_cm2` | Real | n/a | yes | Keyword argument `heat_load_tolerance_j_cm2`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_certify_targeting_candidates`. Returns `(times=certified_times, outcomes=certified_outcomes, failure=:nonfinite_energy, ` or `(times=certified_times, outcomes=certified_outcomes, failure=:heat_load, failure` or `(times=certified_times, outcomes=certified_outcomes, failure=:energy_order, fail` or `(times=certified_times, outcomes=certified_outcomes, failure=:none, failure_time`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control_evaluate_candidate|evaluate_candidate]] · `callees` → `callers` · feedback · `src/gnc/control/targeting_control.jl:986-986`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/targeting_control.jl:891-891`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/control/targeting_control.jl:909-909`
- `callees` → [[gnc.targeting_control_evaluate_candidate|evaluate_candidate]] · `callers` · call · `src/gnc/control/targeting_control.jl:898-898`
<!-- vulcan:connections:end -->

## Limitations
Certification stops at the first violation, so a single noisy prediction early in the sweep truncates the bracket to a few points even if later candidates would have been fine.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 883.
