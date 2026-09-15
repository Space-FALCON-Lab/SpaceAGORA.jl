---
id: simulation.setup__body_query_names_reuse_key
label: _body_query_names_reuse_key
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _body_query_names_reuse_key
  lines:
  - 258
  - 258
inputs:
- id: body_query_names
  type: AbstractVector{String}
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
  type: String
  units: n/a
  description: Return value of `_body_query_names_reuse_key`.
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

# _body_query_names_reuse_key

## Purpose
Collapses the ordered list of third-body names queried by an N-body effector into a single string suitable for use inside a tuple dictionary key.

## Design & Implementation
`_body_query_names_reuse_key(body_query_names::AbstractVector{String})::String` returns `join(body_query_names, '\0')`. The NUL separator cannot appear in a valid SPICE body name, so the encoding is unambiguous and order-preserving. The result becomes the second element of `NBodyEphemerisReuseKey`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `body_query_names` | AbstractVector{String} | n/a | yes | Positional argument `body_query_names`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_body_query_names_reuse_key`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__nbody_ephemeris_reuse_key|_nbody_ephemeris_reuse_key]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:282-282`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Order matters: `["MOON","SUN"]` and `["SUN","MOON"]` produce different keys and therefore separate caches even though the sampled data is the same up to permutation. Case is not normalised, so `"Moon"` and `"MOON"` also differ. An empty vector yields an empty string, which is a valid but degenerate key.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 258.
