---
id: simulation.config__gram_track_cache_target_use_j2
label: _gram_track_cache_target_use_j2
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _gram_track_cache_target_use_j2
  lines:
  - 15
  - 15
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
  description: Return value of `_gram_track_cache_target_use_j2`.
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

# _gram_track_cache_target_use_j2

## Purpose
Selects whether the target trajectory that seeds the GRAM track cache is propagated with J2 oblateness included rather than pure two-body motion.

## Design & Implementation
Delegates to `_parse_bool_env("SPACEAGORA_GRAM_TRACK_CACHE_TARGET_USE_J2", true)` and returns a `Bool`. Like the other knobs in this file it is resolved once per run into `CallbackEnvConfig` so that per-step cache lookups do not re-parse the environment.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_gram_track_cache_target_use_j2`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:182-182`

**Downstream**

- `callees` → [[simulation.config__parse_bool_env|_parse_bool_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:16-16`
<!-- vulcan:connections:end -->

## Limitations
The flag only governs the cache seeding propagation; it does not have to agree with the gravity model actually used by the integrator, so a run can seed the cache with J2 dynamics while integrating two-body, or the reverse, with no warning. Invalid values throw from the underlying parser.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 15.
