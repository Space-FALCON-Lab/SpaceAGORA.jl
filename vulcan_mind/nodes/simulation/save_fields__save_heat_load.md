---
id: simulation.save_fields__save_heat_load
label: _save_heat_load
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_heat_load
  lines:
  - 143
  - 143
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
  description: Return value of `_save_heat_load`. Returns `heat_loads`.
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

# _save_heat_load

## Purpose
Save-time getter for accumulated heat load per spacecraft, reported as the maximum across that spacecraft's links, which is the quantity the aerobraking guidance limits against.

## Design & Implementation
Marked `@inline`. Unlike heat rate, heat load is an integrated state rather than an instantaneous evaluation, so it is read straight out of the solver state by `_simulation_engine_module()._state_heat_loads(u, integrator.p.args, i)`, with `args` needed to locate the heat-load block within the composed state layout. The per-link vector is reduced with `maximum`, and an empty vector, meaning a spacecraft with no thermal links, yields `0.0`. The result is a `Vector{Float64}` of length `num_sats`.

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
| out | `result` | Any | n/a | — | Return value of `_save_heat_load`. Returns `heat_loads`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:185-185`

**Downstream**

- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:146-146`
<!-- vulcan:connections:end -->

## Limitations
Reducing with `maximum` discards the per-link distribution, so a configuration where a different link dominates at different times shows a saved curve that switches sources without any marker. An empty link set is indistinguishable in the output from a spacecraft that has genuinely accumulated no heat. The value is only as good as the heat-rate derivative that fed the integration, including any zeroing applied when an atmosphere sample was rejected.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 143.
