---
id: parallel.env_config__telemetry_bucket
label: _telemetry_bucket
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: _telemetry_bucket
  lines:
  - 166
  - 166
inputs:
- id: source
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `source`.
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
  description: Return value of `_telemetry_bucket`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# _telemetry_bucket

## Purpose
Collapses fine-grained policy source identifiers into the coarse buckets under which timing telemetry and persisted hints are aggregated.

## Design & Implementation
Pure mapping on `source::Symbol`: `:density_callback` and `:density_callback_lockfree` both become `:density`; `:control_callback` becomes `:control`; `:multibody` stays `:multibody`; everything else (including `:dynamic_effectors` and `:thermal_callback`) becomes `:other`. Returns a `Symbol` and never throws.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `source` | Symbol | n/a | yes | Positional argument `source`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_telemetry_bucket`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/env_config.jl`
- [[parallel.policy_telemetry__record_policy_decision_bang|_record_policy_decision!]] · `callees` → `callers` · call · `src/parallel/policy/policy_telemetry.jl:36-36`

**Downstream**

- `callees` → [[parcore.env_config_snapshot_policy_decision_env|snapshot_policy_decision_env]] · `callers` · call · `src/parallel/policy/env_config.jl:178-178`
<!-- vulcan:connections:end -->

## Limitations
Locked and lock-free density models share the `:density` bucket even though `auto_thread_min_budget` treats them differently, so persisted hints learned with one model type are applied to the other. Thermal callbacks fall into `:other` alongside unrelated sources, blending their statistics.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 166.
