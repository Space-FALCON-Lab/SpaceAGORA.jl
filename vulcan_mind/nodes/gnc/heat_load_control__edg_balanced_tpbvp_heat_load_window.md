---
id: gnc.heat_load_control__edg_balanced_tpbvp_heat_load_window
label: _edg_balanced_tpbvp_heat_load_window
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_balanced_tpbvp_heat_load_window
  lines:
  - 460
  - 460
inputs:
- id: config
  type: Any
  units: n/a
  required: true
  description: Positional argument `config`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: track
  type: Any
  units: n/a
  required: true
  description: Positional argument `track`.
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: controlled_panel_links
  type: Tuple{Vararg{Int}}
  units: n/a
  required: true
  description: Positional argument `controlled_panel_links`.
- id: heat_load_j_cm2
  type: Float64
  units: n/a
  required: true
  description: Positional argument `heat_load_j_cm2`.
- id: target_load
  type: Float64
  units: n/a
  required: true
  description: Positional argument `target_load`.
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
  type: Tuple
  units: n/a
  description: Return value of `_edg_balanced_tpbvp_heat_load_window`. Returns `(t
    + first(track.time), t + last(track.time))` or `(Inf, Inf)` or `(t + track.time[best_start],
    t + track.time[best_stop])`.
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

# _edg_balanced_tpbvp_heat_load_window

## Purpose
Chooses a low-drag window directly from heat-rate savings so that the predicted heat load lands just under the target, used as the preferred answer when the TPBVP integration solver is selected.

## Design & Implementation
Signature `(config, p, t, track, spacecraft, controlled_panel_links, heat_load_j_cm2, target_load; heat_rate_control, structural_control)`; returns `(Inf, Inf)` for fewer than three nodes. It builds a constrained all-high profile and an all-low profile, evaluates both heat-rate histories with `heat_rate_control = false`, and computes `required_saving = heat_load_j_cm2 + ∫qdot_high - target_load`. If no saving is needed or the per-node saving `qdot_high - qdot_low` is never positive it returns `(Inf, Inf)`. It then forms the trapezoidal cumulative saving, and if the total is insufficient returns the whole pass. Otherwise it scans every `start_idx`, uses `searchsortedfirst` on `cumulative` to find the earliest `stop_idx` that achieves the saving, and keeps the window with the smallest overshoot (tolerance 0.02 J/cm^2), preferring later stops on ties. Returns absolute times `(t + time[start], t + time[stop])`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `track` | Any | n/a | yes | Positional argument `track`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `controlled_panel_links` | Tuple{Vararg{Int}} | n/a | yes | Positional argument `controlled_panel_links`. |
| in | `heat_load_j_cm2` | Float64 | n/a | yes | Positional argument `heat_load_j_cm2`. |
| in | `target_load` | Float64 | n/a | yes | Positional argument `target_load`. |
| in | `heat_rate_control` | Bool | n/a | yes | Keyword argument `heat_rate_control`. |
| in | `structural_control` | Bool | n/a | yes | Keyword argument `structural_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_edg_balanced_tpbvp_heat_load_window`. Returns `(t + first(track.time), t + last(track.time))` or `(Inf, Inf)` or `(t + track.time[best_start], t + track.time[best_stop])`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control_residual|residual]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:725-725`
- [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:725-725`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_constrained_heat_load_alpha_profile|_edg_constrained_heat_load_alpha_profile]] · `callers` · call · `src/gnc/control/heat_load_control.jl:475-475`
- `callees` → [[gnc.heat_load_control__edg_integrate_series|_edg_integrate_series]] · `callers` · call · `src/gnc/control/heat_load_control.jl:488-488`
- `callees` → [[gnc.heat_load_control__edg_profile_heat_rates|_edg_profile_heat_rates]] · `callers` · call · `src/gnc/control/heat_load_control.jl:486-486`
<!-- vulcan:connections:end -->

## Limitations
Per-node savings are computed assuming the two profiles do not change the trajectory, so drag feedback on the track is ignored. The O(n log n) scan is cheap but only considers a single contiguous window. The 0.02 J/cm^2 tie tolerance is hard-coded.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 460.
