---
id: simulation.config__gram_entry_target_max_steps
label: _gram_entry_target_max_steps
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _gram_entry_target_max_steps
  lines:
  - 36
  - 36
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
  description: Return value of `_gram_entry_target_max_steps`.
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

# _gram_entry_target_max_steps

## Purpose
Caps the number of integration steps taken when building the analytic entry-targeting reference trajectory.

## Design & Implementation
Reads and strips `SPACEAGORA_GRAM_ENTRY_TARGET_MAX_STEPS`, defaulting to the string `"400"`, and parses it with `parse(Int, raw)`. A parse failure is caught and rethrown as an `ArgumentError` quoting the raw string. The parsed value is then floored with `max(8, parsed)` so at least eight steps always run.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_gram_entry_target_max_steps`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.targeting_gram_entry_target_allen_eggers|_gram_entry_target_allen_eggers]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:261-261`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Floating-point spellings such as `"400.0"` fail to parse and abort the run rather than being truncated. The cap is on step count alone, so with the default half-second step the reference arc spans at most two hundred seconds of flight and a longer pass is silently truncated rather than reported as incomplete.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 36.
