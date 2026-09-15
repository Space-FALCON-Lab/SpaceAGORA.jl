---
id: gnc.heat_load_control__edg_low_alpha_switch_window
label: _edg_low_alpha_switch_window
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_low_alpha_switch_window
  lines:
  - 449
  - 449
inputs:
- id: config
  type: Any
  units: n/a
  required: true
  description: Positional argument `config`.
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
- id: alpha_profile
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `alpha_profile`.
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
  description: Return value of `_edg_low_alpha_switch_window`. Returns `(Inf, Inf)`
    or `(t + track.time[first_low], t + track.time[last_low])`.
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

# _edg_low_alpha_switch_window

## Purpose
Converts the first low-alpha interval of a profile into absolute switch-on and switch-off times for the low-drag window.

## Design & Implementation
Signature `(config, t::Float64, track, alpha_profile)`. Uses `_edg_first_low_alpha_interval_indices`; returns `(Inf, Inf)` if no interval exists or if it is a single node (`first_low == last_low`). Otherwise returns `(t + track.time[first_low], t + track.time[last_low])` in seconds of absolute simulation time.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `track` | Any | n/a | yes | Positional argument `track`. |
| in | `alpha_profile` | Vector{Float64} | n/a | yes | Positional argument `alpha_profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_edg_low_alpha_switch_window`. Returns `(Inf, Inf)` or `(t + track.time[first_low], t + track.time[last_low])`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control_residual|residual]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:741-741`
- [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:741-741`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_first_low_alpha_interval_indices|_edg_first_low_alpha_interval_indices]] · `callers` · call · `src/gnc/control/heat_load_control.jl:450-450`
<!-- vulcan:connections:end -->

## Limitations
A one-node window is treated as no window, so very short optimal dips are dropped. The window ends at the last low node rather than at the following high node, shortening it by one grid step.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 449.
