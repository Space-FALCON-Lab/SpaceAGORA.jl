---
id: simulation.setup__parse_nonnegative_int_env
label: _parse_nonnegative_int_env
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _parse_nonnegative_int_env
  lines:
  - 181
  - 181
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: Int
  units: n/a
  required: true
  description: Positional argument `default`.
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
  description: Return value of `_parse_nonnegative_int_env`.
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

# _parse_nonnegative_int_env

## Purpose
Parses an integer engine parameter that may be zero (for example a cache entry limit where 0 disables retention) but never negative, throwing rather than clamping so misconfiguration is visible.

## Design & Implementation
`_parse_nonnegative_int_env(name::String, default::Int)::Int` reads `_engine_env_get(name, string(default))`, strips, and calls `parse(Int, raw)` inside `try`/`catch`, converting parse failures to `ArgumentError("<name> must be an integer value, got '<raw>'")`. It then requires `parsed >= 0` or throws `ArgumentError("<name> must be >= 0, got <parsed>")`. Unlike `ParallelPolicy.parse_nonnegative_int_env`, which silently clamps with `max(0, ...)`, this variant rejects negatives. Used by `_ephemeris_reuse_max_entries`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | Int | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_parse_nonnegative_int_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__ephemeris_reuse_max_entries|_ephemeris_reuse_max_entries]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:237-237`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/setup.jl:182-182`
<!-- vulcan:connections:end -->

## Limitations
Float-formatted integers such as `"32.0"` fail to parse. Very large digit strings raise `OverflowError` which is caught and reported as a generic integer error. No upper bound is enforced.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 181.
