---
id: simulation.thermal_callbacks__heat_rate_buffer_for_sat_bang
label: _heat_rate_buffer_for_sat!
kind: function
source:
  file: src/simulation/callbacks/thermal_callbacks.jl
  symbol: _heat_rate_buffer_for_sat!
  lines:
  - 1
  - 1
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
  description: Return value of `_heat_rate_buffer_for_sat!`; mutates `p` in place.
    Returns `heat_rates`.
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

# _heat_rate_buffer_for_sat!

## Purpose
Returns the zeroed per-link heat-rate scratch vector for spacecraft `sat_idx`, resizing it when the spacecraft's link count has changed. It mutates `p.shared_buffers.heat_rates[sat_idx]` in place and returns that same vector.

## Design & Implementation
Marked `@inline`. It reads `p.args.dynamics_model.spacecraft[sat_idx].links`, compares `length(heat_rates)` against `n_links`, and calls `resize!` only on mismatch, so the steady-state path performs no allocation. It then unconditionally `fill!`s the buffer with `0.0`, guaranteeing that link entries skipped later by the caller read as zero heat rate rather than as a stale value from the previous callback firing.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_heat_rate_buffer_for_sat!`; mutates `p` in place. Returns `heat_rates`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/thermal_callbacks.jl`
- [[simulation.thermal_callbacks__compute_stage_heat_rates_bang|_compute_stage_heat_rates!]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:20-20`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The buffer is shared mutable state indexed by spacecraft, so two tasks working on the same `sat_idx` would corrupt each other; safety depends entirely on the caller partitioning spacecraft across threads. No bounds check guards `sat_idx`, so an out-of-range index throws from the array access rather than producing a clear message. Growing the buffer with `resize!` leaves the added slots undefined until the `fill!` on the following line.

## Provenance
Mapped from `src/simulation/callbacks/thermal_callbacks.jl` line 1.
