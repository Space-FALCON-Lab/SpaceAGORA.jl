---
id: simulation.vacuum_predicted_gram__vacuum_gram_cache_npoints
label: _vacuum_gram_cache_npoints
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl
  symbol: _vacuum_gram_cache_npoints
  lines:
  - 22
  - 22
inputs:
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
  type: Int
  units: n/a
  description: Return value of `_vacuum_gram_cache_npoints`.
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

# _vacuum_gram_cache_npoints

## Purpose
Returns the number of spline knots the vacuum-predicted GRAM cache should sample across its look-ahead horizon, read from the environment so users can trade cache build cost (one density-model call per knot) against interpolation accuracy without changing code.

## Design & Implementation
An `@inline` function returning `Int`. It reads `get(ENV, "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS", "20")`, strips whitespace, and parses with `parse(Int, raw)` inside a `try` block; a parse failure is rethrown as `ArgumentError("SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS must be an integer, got '<raw>'")`. The parsed value is clamped from below with `max(4, parsed)` so that `_natural_cubic_spline_build!` always has at least two interior unknowns. The result feeds `_build_vacuum_gram_cache!` where knot spacing is `h = horizon_s / (n_pts - 1)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_vacuum_gram_cache_npoints`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:185-185`

**Downstream**

- `callees` → [[simulation.config__parse_float_env|_parse_float_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:33-33`
<!-- vulcan:connections:end -->

## Limitations
There is no upper bound, so an absurdly large value produces that many synchronous density-model calls on every cache rebuild. Values below 4 are silently raised to 4 rather than rejected, which can surprise a user who deliberately requested a coarser cache. The environment lookup happens on each call; callers should cache the result outside the integration loop. Negative or zero inputs are accepted and clamped rather than reported.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl` line 22.
