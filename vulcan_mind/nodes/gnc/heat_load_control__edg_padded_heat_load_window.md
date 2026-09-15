---
id: gnc.heat_load_control__edg_padded_heat_load_window
label: _edg_padded_heat_load_window
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_padded_heat_load_window
  lines:
  - 526
  - 526
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
- id: window
  type: Tuple{Float64, Float64}
  units: n/a
  required: true
  description: Positional argument `window`.
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
  description: Return value of `_edg_padded_heat_load_window`. Returns `(`.
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

# _edg_padded_heat_load_window

## Purpose
Widens or clips the solved low-drag window to the prediction horizon, adding safety margin around it for the closed-form solver where trajectory error is larger.

## Design & Implementation
`@inline` function `(config, t, track, window::Tuple{Float64,Float64})`. Non-finite windows are returned unchanged. When `config.heat_load_switch_solver == :tpbvp_integration` the window is simply clipped to `[t + first(time), t + last(time)]`. Otherwise the start is moved 10 s earlier and the end is extended by `max(10, 0.15 * duration)` seconds, both then clipped to the track span.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `track` | Any | n/a | yes | Positional argument `track`. |
| in | `window` | Tuple{Float64, Float64} | n/a | yes | Positional argument `window`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_padded_heat_load_window`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control_residual|residual]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:743-743`
- [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:743-743`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The 10 s and 15 percent pads are hard-coded and asymmetric (larger at the end), reflecting an assumption that the closed-form track exits the atmosphere too early. Padding always lengthens the low-drag window, which reduces energy depletion relative to the optimum.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 526.
