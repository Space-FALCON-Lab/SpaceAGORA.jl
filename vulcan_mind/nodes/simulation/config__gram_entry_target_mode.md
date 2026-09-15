---
id: simulation.config__gram_entry_target_mode
label: _gram_entry_target_mode
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _gram_entry_target_mode
  lines:
  - 24
  - 24
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
  type: Symbol
  units: n/a
  description: Return value of `_gram_entry_target_mode`.
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

# _gram_entry_target_mode

## Purpose
Chooses the analytic model used to pick an entry-targeting reference trajectory for GRAM atmospheric sampling during entry passes.

## Design & Implementation
Reads `SPACEAGORA_GRAM_ENTRY_TARGET_MODE`, defaulting to `"allen_eggers"`, then lowercases and strips it. Maps `off`, `none`, `0`, `false`, `no` to `:off`; maps `allen_eggers`, `allen-eggers`, `allen`, `ae`, `on`, `1`, `true`, `yes`, `auto` to `:allen_eggers`. Any other token throws `ArgumentError` listing the two supported modes.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_gram_entry_target_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:99-99`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only one analytic entry model is implemented, so `:allen_eggers` is the sole non-trivial branch and the Allen-Eggers assumptions (exponential atmosphere, steep ballistic entry, negligible lift) are inherited without check. The alias `auto` maps to the same value rather than selecting a model from flight conditions.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 24.
