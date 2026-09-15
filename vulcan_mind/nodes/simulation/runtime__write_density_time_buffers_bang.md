---
id: simulation.runtime__write_density_time_buffers_bang
label: _write_density_time_buffers!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/runtime.jl
  symbol: _write_density_time_buffers!
  lines:
  - 35
  - 35
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
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
  type: Nothing
  units: n/a
  description: Return value of `_write_density_time_buffers!`; mutates `p` in place.
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

# _write_density_time_buffers!

## Purpose
Stamps the sample time for every satellite after a batch density evaluation, which writes densities directly and would otherwise leave the time buffer stale.

## Design & Implementation
Takes `min(num_sats, length(times))` as the loop bound and fills `density_sample_t[1:limit]` with `t` under `@inbounds`, safe because the bound was clamped. Separating this from the per-satellite writer is what lets the batch path fill densities through `getDensityBatch!` and still keep the time buffer consistent.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_write_density_time_buffers!`; mutates `p` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:334-334`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:334-334`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It assumes every satellite was sampled at exactly `t`, which holds for the batch path but would be wrong if reused after a partial update.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/runtime.jl` line 35.
