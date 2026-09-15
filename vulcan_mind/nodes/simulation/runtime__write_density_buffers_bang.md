---
id: simulation.runtime__write_density_buffers_bang
label: _write_density_buffers!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/runtime.jl
  symbol: _write_density_buffers!
  lines:
  - 12
  - 12
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
- id: rho
  type: Float64
  units: n/a
  required: true
  description: Positional argument `rho`.
- id: T
  type: Float64
  units: n/a
  required: true
  description: Positional argument `T`.
- id: wind_vec
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `wind_vec`.
- id: t
  type: Float64
  units: n/a
  required: false
  description: Positional argument `t` (default `p.shared_buffers.current_time[]`).
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
  type: Nothing
  units: n/a
  description: Return value of `_write_density_buffers!`; mutates `p` in place. Returns
    `nothing`.
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

# _write_density_buffers!

## Purpose
Stores one satellite's freshly sampled atmosphere into the shared buffers the right-hand side reads, together with the sample time.

## Design & Implementation
Writes `rho`, `T`, `wind_vec` and `t` into `p.shared_buffers.densities`, `temperatures`, `winds` and `density_sample_t` at `sat_idx`, each write guarded by a length check so a buffer sized for fewer satellites is skipped rather than overrun. The time defaults to `p.shared_buffers.current_time[]` when the caller does not pass it explicitly. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `rho` | Float64 | n/a | yes | Positional argument `rho`. |
| in | `T` | Float64 | n/a | yes | Positional argument `T`. |
| in | `wind_vec` | SVector{3, Float64} | n/a | yes | Positional argument `wind_vec`. |
| in | `t` | Float64 | n/a | no | Positional argument `t` (default `p.shared_buffers.current_time[]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_write_density_buffers!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime_update_density_sat_bang|update_density_sat!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:243-243`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:243-243`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The four guards are independent, so mismatched buffer lengths produce a partially updated set — a density with a stale temperature — with no indication that any write was skipped.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/runtime.jl` line 12.
