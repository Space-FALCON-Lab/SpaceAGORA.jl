---
id: simulation.setup__cache_from_nbody_ephemeris_payload
label: _cache_from_nbody_ephemeris_payload
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _cache_from_nbody_ephemeris_payload
  lines:
  - 1631
  - 1631
inputs:
- id: payload
  type: Any
  units: n/a
  required: true
  description: Positional argument `payload`.
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
  description: Return value of `_cache_from_nbody_ephemeris_payload`. Returns `(`.
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

# _cache_from_nbody_ephemeris_payload

## Purpose
Reconstructs an N-body ephemeris cache and its build parameters from a deserialized payload, validating every field so a corrupt or stale file cannot silently poison a run.

## Design & Implementation
Requires a `NamedTuple`, checks `schema_version` against `NBODY_EPHEMERIS_CACHE_SCHEMA_VERSION`, and reads the primary name, body names, `et_start`, `mission_end_s`, `dt_s`, the time vector and the position matrix through `_payload_field`. It then validates finiteness and positivity of the scalars, non-empty bodies, at least two samples, and that the position matrix dimensions match samples by bodies. Returns a named tuple of the rebuilt cache plus the three build parameters.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `payload` | Any | n/a | yes | Positional argument `payload`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_cache_from_nbody_ephemeris_payload`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__load_nbody_ephemeris_cache_bang|_load_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1668-1668`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/setup.jl:1641-1641`
- `callees` → [[simulation.setup__nbody_ephemeris_cache_from_samples|_nbody_ephemeris_cache_from_samples]] · `callers` · call · `src/simulation/engine/setup.jl:1655-1655`
- `callees` → [[simulation.setup__payload_field|_payload_field]] · `callers` · call · `src/simulation/engine/setup.jl:1634-1634`
<!-- vulcan:connections:end -->

## Limitations
It does not verify that the times are uniformly spaced at `dt_s` or that they begin at `et_start`, so a file whose header disagrees with its samples is accepted.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1631.
