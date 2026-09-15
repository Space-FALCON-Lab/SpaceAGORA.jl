---
id: simulation.thermal_callbacks_update_thermal_sat_bang
label: update_thermal_sat!
kind: function
source:
  file: src/simulation/callbacks/thermal_callbacks.jl
  symbol: update_thermal_sat!
  lines:
  - 61
  - 61
inputs:
- id: i
  type: Int
  units: n/a
  required: true
  description: Positional argument `i`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
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
  description: Return value of `update_thermal_sat!`; mutates `i` in place. Returns
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

# update_thermal_sat!

## Purpose
Per-spacecraft body of the thermal discrete callback: refreshes the heat rates for spacecraft `i` from the current integrator state at time `t`.

## Design & Implementation
Defined as a closure inside `get_thermal_callback` so it captures nothing beyond its arguments and can be invoked directly by both the serial loop and the threaded fan-out. It forwards `u.sc[i]`, the per-spacecraft slice of the composed state, to `_compute_stage_heat_rates!(p, u.sc[i], i, t; use_buffered_density=true)` and returns `nothing`. The buffered density flag matters: the callback runs after the density callbacks have populated the shared atmosphere buffers, so it reuses those samples instead of paying for another atmosphere evaluation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `update_thermal_sat!`; mutates `i` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.thermal_callbacks_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:76-76`
- [[simulation_a.thermal_callbacks_get_thermal_callback|get_thermal_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:61-61`

**Downstream**

- `callees` → [[simulation.thermal_callbacks__compute_stage_heat_rates_bang|_compute_stage_heat_rates!]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:62-62`
<!-- vulcan:connections:end -->

## Limitations
Using buffered density couples this callback to callback ordering; if the density buffers were not refreshed at this step the heat rates are computed from the previous step's atmosphere. The state slice `u.sc[i]` must exist, so a spacecraft index beyond the composed state throws. The return value carries nothing, so callers cannot tell whether a zeroed buffer means no heating or a rejected atmosphere sample.

## Provenance
Mapped from `src/simulation/callbacks/thermal_callbacks.jl` line 61.
