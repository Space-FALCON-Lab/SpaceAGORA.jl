---
id: simulation.save_fields__save_positions
label: _save_positions
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_positions
  lines:
  - 17
  - 17
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
  description: Return value of `_save_positions`. Returns `positions`.
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

# _save_positions

## Purpose
Save-time getter for spacecraft inertial positions. It returns a `Vector{SVector{3, Float64}}` of length `num_sats` holding each spacecraft's position in the inertial frame, in metres.

## Design & Implementation
Marked `@inline`. It allocates the result with `Vector{SVector{3, Float64}}(undef, num_sats)` and fills it in an `@inbounds` loop calling `_simulation_engine_module()._state_position_ii(u, i)`, which decodes the composed solver state `u` for spacecraft `i` without copying the underlying arrays. Going through `_simulation_engine_module()` keeps the callback module free of a hard compile-time dependency on the engine. The `t` and `integrator` arguments are accepted to match the uniform getter signature but are unused, because position is read directly from the state vector.

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
| out | `result` | Any | n/a | — | Return value of `_save_positions`. Returns `positions`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:173-173`

**Downstream**

- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:20-20`
<!-- vulcan:connections:end -->

## Limitations
A fresh vector is allocated on every save point, so long runs with dense saving generate significant garbage. The `@inbounds` loop trusts that `u` actually contains `num_sats` spacecraft; a state shorter than the configured constellation reads out of bounds rather than erroring. The values are whatever frame `_state_position_ii` returns, with no conversion or validation applied here.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 17.
