---
id: simulation.config__gram_entry_target_cd
label: _gram_entry_target_cd
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _gram_entry_target_cd
  lines:
  - 34
  - 34
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
  description: Return value of `_gram_entry_target_cd`. Returns `max(0.05, _parse_float_env("SPACEAGORA_GRAM_ENTRY_TARGET_CD",
    1.5))`.
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

# _gram_entry_target_cd

## Purpose
Supplies the drag coefficient used by the analytic entry-targeting trajectory that seeds GRAM atmospheric queries.

## Design & Implementation
Evaluates `_parse_float_env("SPACEAGORA_GRAM_ENTRY_TARGET_CD", 1.5)` and clamps the result from below with `max(0.05, ...)`. The default of 1.5 is a conventional blunt-body ballistic value; the 0.05 floor exists so a zero or negative override cannot produce a division by zero or a non-decelerating reference.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gram_entry_target_cd`. Returns `max(0.05, _parse_float_env("SPACEAGORA_GRAM_ENTRY_TARGET_CD", 1.5))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:115-115`
- [[simulation_a.targeting_gram_entry_target_allen_eggers|_gram_entry_target_allen_eggers]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:231-231`

**Downstream**

- `callees` → [[simulation.config__parse_float_env|_parse_float_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:34-34`
<!-- vulcan:connections:end -->

## Limitations
The clamp is one-sided: absurdly large drag coefficients pass through unmodified. The value is a single scalar applied to the whole reference arc, so Mach and Knudsen dependence of the real drag coefficient is not represented, and it need not match the drag coefficient used by the actual simulated spacecraft.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 34.
