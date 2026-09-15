---
id: gnc.heat_load_control_residual
label: residual
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: residual
  lines:
  - 648
  - 648
inputs:
- id: k
  type: Any
  units: n/a
  required: true
  description: Positional argument `k`.
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
  description: Return value of `residual`. Returns `f`.
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

# residual

## Purpose
Nested closure inside `_edg_solve_heat_load_switches` that maps a heating weight `k` to the signed difference between the predicted total heat load and the target, serving as the root-finding objective.

## Design & Implementation
Captures `config`, `p`, `spacecraft`, `pos`, `vel`, `mass`, `t`, `env`, `coeffs`, the control flags, and five `Ref` cells. For each `k` it calls `_edg_heat_load_profile_for_k`, builds a constrained all-high profile, projects the raw profile to two switches with `_edg_first_two_switch_alpha_profile`, and stores the track and profile in `last_track[]`/`last_profile[]`. The predicted load is `heat_load_j_cm2 + _edg_profile_heat_load(config, p, track, profile; heat_rate_control = false)`, and `f = predicted_load - target_load`. Whenever `f <= 0` and `|f|` beats `best_under_residual[]`, the track and a copy of the profile are saved in the `best_under_*` refs. Returns `f`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `k` | Any | n/a | yes | Positional argument `k`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `residual`. Returns `f`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:648-648`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_balanced_tpbvp_heat_load_window|_edg_balanced_tpbvp_heat_load_window]] · `callers` · call · `src/gnc/control/heat_load_control.jl:725-725`
- `callees` → [[gnc.heat_load_control__edg_constrained_heat_load_alpha_profile|_edg_constrained_heat_load_alpha_profile]] · `callers` · call · `src/gnc/control/heat_load_control.jl:664-664`
- `callees` → [[gnc.heat_load_control__edg_drag_passage_duration|_edg_drag_passage_duration]] · `callers` · call · `src/gnc/control/heat_load_control.jl:750-750`
- `callees` → [[gnc.heat_load_control__edg_first_two_switch_alpha_profile|_edg_first_two_switch_alpha_profile]] · `callers` · call · `src/gnc/control/heat_load_control.jl:674-674`
- `callees` → [[gnc.heat_load_control__edg_heat_load_profile_for_k|_edg_heat_load_profile_for_k]] · `callers` · call · `src/gnc/control/heat_load_control.jl:649-649`
- `callees` → [[gnc.heat_load_control__edg_low_alpha_switch_window|_edg_low_alpha_switch_window]] · `callers` · call · `src/gnc/control/heat_load_control.jl:741-741`
- `callees` → [[gnc.heat_load_control__edg_padded_heat_load_window|_edg_padded_heat_load_window]] · `callers` · call · `src/gnc/control/heat_load_control.jl:743-743`
- `callees` → [[gnc.heat_load_control__edg_profile_heat_load|_edg_profile_heat_load]] · `callers` · call · `src/gnc/control/heat_load_control.jl:677-677`
<!-- vulcan:connections:end -->

## Limitations
The function is not monotonic in `k` in general because the two-switch projection can change discontinuously, so Brent may fail (the caller falls back to bisection). Side effects through the `Ref` cells mean the closure is not reentrant and must not be evaluated concurrently. Each evaluation is a full trajectory prediction.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 648.
