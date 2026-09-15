---
id: simulation.save_fields__save_heat_rate
label: _save_heat_rate
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_heat_rate
  lines:
  - 127
  - 127
inputs:
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
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
  description: Return value of `_save_heat_rate`. Returns `heat_rates`.
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

# _save_heat_rate

## Purpose
Save-time getter for the peak convective heat rate on each spacecraft, reported as the maximum over that spacecraft's links so a single scalar per spacecraft captures the worst-case heating.

## Design & Implementation
Marked `@inline`, with three source paths per spacecraft. When the state exposes `u.sc`, the rates are recomputed fresh by `_compute_stage_heat_rates!(integrator.p, u.sc[i], i, Float64(t); use_buffered_density=false)`, so the save reflects a new atmosphere sample rather than the buffered one. Otherwise, if `i <= length(shared_heat_rates)`, the buffered per-link vector is reused. Failing both, an empty `Float64[]` is used. Each case ends with `heat_rates[i] = !isempty(rates) ? maximum(rates) : 0.0`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_save_heat_rate`. Returns `heat_rates`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:184-184`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:132-132`
- `callees` → [[simulation.thermal_callbacks__compute_stage_heat_rates_bang|_compute_stage_heat_rates!]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:132-132`
<!-- vulcan:connections:end -->

## Limitations
The primary path recomputes the atmosphere at every save point, which for a dense save schedule is a substantial extra cost and can differ from the atmosphere the dynamics actually integrated against. Because `_compute_stage_heat_rates!` writes into `p.shared_buffers.heat_rates`, saving mutates the same buffer the thermal callback owns, which is a hazard if saving ever runs concurrently with that callback. An empty link set and a genuinely zero heat rate both save as `0.0`.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 127.
