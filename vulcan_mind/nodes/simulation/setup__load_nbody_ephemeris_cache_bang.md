---
id: simulation.setup__load_nbody_ephemeris_cache_bang
label: _load_nbody_ephemeris_cache!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _load_nbody_ephemeris_cache!
  lines:
  - 1664
  - 1664
inputs:
- id: path
  type: String
  units: n/a
  required: true
  description: Positional argument `path`.
- id: replace
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `replace` (default `true`).
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
  type: SimulationModel.NBodyEphemerisCache
  units: n/a
  description: Return value of `_load_nbody_ephemeris_cache!`; mutates `path` in place.
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

# _load_nbody_ephemeris_cache!

## Purpose
Loads a serialized N-body ephemeris cache from disk and registers it as prewarmed so subsequent runs with matching parameters pick it up without rebuilding.

## Design & Implementation
Deserializes the file, validates it through `_cache_from_nbody_ephemeris_payload`, and registers the cache under its key via `_register_prewarmed_nbody_ephemeris_cache!` with the `replace` flag forwarded. Returns the cache.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | String | n/a | yes | Positional argument `path`. |
| in | `replace` | Bool | n/a | no | Keyword argument `replace` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.NBodyEphemerisCache | n/a | — | Return value of `_load_nbody_ephemeris_cache!`; mutates `path` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.public_api_load_nbody_ephemeris_cache_bang|load_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/public_api.jl:54-54`

**Downstream**

- `callees` → [[simulation.setup__cache_from_nbody_ephemeris_payload|_cache_from_nbody_ephemeris_payload]] · `callers` · call · `src/simulation/engine/setup.jl:1668-1668`
- `callees` → [[simulation.setup__register_prewarmed_nbody_ephemeris_cache_bang|_register_prewarmed_nbody_ephemeris_cache!]] · `callers` · call · `src/simulation/engine/setup.jl:1669-1669`
<!-- vulcan:connections:end -->

## Limitations
The file is trusted after schema validation; there is no digest check, and `deserialize` on an untrusted file can execute arbitrary code.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1664.
