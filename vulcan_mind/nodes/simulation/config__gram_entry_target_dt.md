---
id: simulation.config__gram_entry_target_dt
label: _gram_entry_target_dt
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _gram_entry_target_dt
  lines:
  - 35
  - 35
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
  type: Any
  units: n/a
  description: Return value of `_gram_entry_target_dt`. Returns `max(0.05, _parse_float_env("SPACEAGORA_GRAM_ENTRY_TARGET_DT_S",
    0.5))`.
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

# _gram_entry_target_dt

## Purpose
Supplies the fixed propagation time step, in seconds, for the analytic entry-targeting reference trajectory.

## Design & Implementation
Evaluates `_parse_float_env("SPACEAGORA_GRAM_ENTRY_TARGET_DT_S", 0.5)` and applies `max(0.05, ...)`, so the step is never smaller than fifty milliseconds. Combined with `_gram_entry_target_max_steps` this bounds the wall-clock cost of building the reference arc.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gram_entry_target_dt`. Returns `max(0.05, _parse_float_env("SPACEAGORA_GRAM_ENTRY_TARGET_DT_S", 0.5))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.targeting_gram_entry_target_allen_eggers|_gram_entry_target_allen_eggers]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:260-260`

**Downstream**

- `callees` → [[simulation.config__parse_float_env|_parse_float_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:35-35`
<!-- vulcan:connections:end -->

## Limitations
Because the floor is one-sided, a large override yields a coarse arc that can step over the density peak of an entry pass entirely. The step is uniform in time rather than in altitude or dynamic pressure, so resolution is poorest exactly where the atmosphere varies fastest.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 35.
