---
id: simulation.save_fields__save_quaternion
label: _save_quaternion
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_quaternion
  lines:
  - 160
  - 160
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
  description: Return value of `_save_quaternion`. Returns `quaternions`.
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

# _save_quaternion

## Purpose
Save-time getter for spacecraft attitude, returning one `SVector{4, Float64}` attitude quaternion per spacecraft. It is only installed when the mission configuration enables orientation simulation.

## Design & Implementation
Marked `@inline`. For each spacecraft it calls `_simulation_engine_module()._state_quaternion(u, i)`, which returns `nothing` when the composed state carries no orientation block. That case is rejected immediately with `throw(ArgumentError("Quaternion save field requires orientation state."))` rather than being papered over with an identity rotation, so a misconfigured run fails loudly at the first save point. `default_save_fields` only pushes this field when `args.mission_configuration.orientation_sim` is true, and gives it the column prefix `"q"`.

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
| out | `result` | Any | n/a | — | Return value of `_save_quaternion`. Returns `quaternions`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:188-188`

**Downstream**

- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:163-163`
<!-- vulcan:connections:end -->

## Limitations
The throw occurs at save time, not at configuration time, so an explicitly supplied `save_fields` list containing `:quaternion` for a translation-only run fails partway through the integration after work is already done. The quaternion is saved exactly as the state holds it, with no renormalisation and no sign-continuity fixing, so accumulated norm drift and hemisphere flips appear in the output and can confuse downstream interpolation.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 160.
