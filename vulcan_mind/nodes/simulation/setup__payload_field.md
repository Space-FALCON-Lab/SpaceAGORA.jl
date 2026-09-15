---
id: simulation.setup__payload_field
label: _payload_field
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _payload_field
  lines:
  - 1626
  - 1626
inputs:
- id: payload
  type: Any
  units: n/a
  required: true
  description: Positional argument `payload`.
- id: field
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `field`.
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
  description: Return value of `_payload_field`. Returns `getproperty(payload, field)`.
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

# _payload_field

## Purpose
Reads one required field from a deserialised N-body cache payload, producing an `ArgumentError` that names the missing field instead of an opaque property error.

## Design & Implementation
Tests `hasproperty(payload, field)` and raises with the field name in the message if absent, otherwise returns `getproperty`. Declared `@inline`. Every field read in `_cache_from_nbody_ephemeris_payload` goes through it, so a truncated or hand-edited file reports exactly which key it lacks.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `payload` | Any | n/a | yes | Positional argument `payload`. |
| in | `field` | Symbol | n/a | yes | Positional argument `field`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_payload_field`. Returns `getproperty(payload, field)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__cache_from_nbody_ephemeris_payload|_cache_from_nbody_ephemeris_payload]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1634-1634`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It does not check the field's type; the caller's `String`, `Float64` or `Matrix` conversion may still raise a less specific error for a field that is present but mistyped.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1626.
