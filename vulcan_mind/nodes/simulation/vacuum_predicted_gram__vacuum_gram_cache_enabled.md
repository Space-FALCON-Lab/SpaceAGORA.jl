---
id: simulation.vacuum_predicted_gram__vacuum_gram_cache_enabled
label: _vacuum_gram_cache_enabled
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl
  symbol: _vacuum_gram_cache_enabled
  lines:
  - 18
  - 18
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
  type: Bool
  units: n/a
  description: Return value of `_vacuum_gram_cache_enabled`.
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

# _vacuum_gram_cache_enabled

## Purpose
Feature switch for the vacuum-predicted GRAM density cache. The density callback consults it to decide whether to route GRAM queries through the spline cache or call the density model directly; the cache is opt-in because it trades exactness for speed.

## Design & Implementation
An `@inline` function with no arguments returning `Bool`. It delegates to `_parse_bool_env("SPACEAGORA_VACUUM_GRAM_CACHE", false)`, so the cache is disabled unless the environment variable is set to a truthy value recognised by the shared parser. The sibling accessors follow the same pattern: `_vacuum_gram_cache_npoints()` (env `..._NPOINTS`, default 20, floored at 4), `_vacuum_gram_cache_horizon_s()` (default 600 s, floored at 10 s) and `_vacuum_gram_cache_deviation_m()` (default 5000 m, floored at 100 m).

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_vacuum_gram_cache_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:184-184`

**Downstream**

- `callees` → [[simulation.config__parse_bool_env|_parse_bool_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:19-19`
<!-- vulcan:connections:end -->

## Limitations
The environment is read on every call rather than once at configuration time, so toggling the variable mid-process changes behaviour between RHS evaluations and adds an `ENV` lookup to the hot path unless the caller hoists it. The truthy-string rules are whatever `_parse_bool_env` implements and are not documented here. There is no per-configuration override; the switch is global to the Julia process, which matters for concurrent ensembles with different needs.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl` line 18.
