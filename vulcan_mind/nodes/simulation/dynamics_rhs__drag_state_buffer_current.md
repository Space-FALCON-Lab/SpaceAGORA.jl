---
id: simulation.dynamics_rhs__drag_state_buffer_current
label: _drag_state_buffer_current
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _drag_state_buffer_current
  lines:
  - 1937
  - 1937
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: Bool
  units: n/a
  description: Return value of `_drag_state_buffer_current`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _drag_state_buffer_current

## Purpose
Tests whether a satellite's in-atmosphere flag was staged by the drag-state callback at exactly the current RHS time, so the implicit-atmosphere short-circuit can trust it.

## Design & Implementation
Returns true when `sat_idx` is within `in_atmosphere_sample_t` and that entry equals `t`. The timestamp is `NaN` until the callback has fired, so the test is false on the first evaluation. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_drag_state_buffer_current`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__spacecraft_outside_atmosphere_for_current_state|_spacecraft_outside_atmosphere_for_current_state]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1951-1951`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Exact floating-point equality, so RK stage times other than the accepted-step time miss and fall back to a fresh altitude computation.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1937.
