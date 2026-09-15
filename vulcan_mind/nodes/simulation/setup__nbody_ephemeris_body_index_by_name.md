---
id: simulation.setup__nbody_ephemeris_body_index_by_name
label: _nbody_ephemeris_body_index_by_name
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _nbody_ephemeris_body_index_by_name
  lines:
  - 1517
  - 1517
inputs:
- id: body_query_names
  type: Vector{String}
  units: n/a
  required: true
  description: Positional argument `body_query_names`.
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
  type: Dict{String,
  units: n/a
  description: Return value of `_nbody_ephemeris_body_index_by_name`.
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

# _nbody_ephemeris_body_index_by_name

## Purpose
Builds the name-to-column dictionary an `NBodyEphemerisCache` uses to find which column of its position matrix holds a given third body.

## Design & Implementation
Iterates `pairs(body_query_names)` and inserts each SPICE query name with its one-based index into a fresh `Dict{String,Int}`. The dictionary is stored alongside the name vector in the cache so a query can go from body name to column without a linear search.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `body_query_names` | Vector{String} | n/a | yes | Positional argument `body_query_names`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Dict{String, | n/a | — | Return value of `_nbody_ephemeris_body_index_by_name`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__nbody_ephemeris_cache_from_samples|_nbody_ephemeris_cache_from_samples]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1531-1531`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A duplicated name maps to its last index silently; the collector upstream deduplicates names so this does not arise from configuration, but a hand-built cache could trigger it.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1517.
