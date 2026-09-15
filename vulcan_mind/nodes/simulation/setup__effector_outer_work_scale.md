---
id: simulation.setup__effector_outer_work_scale
label: _effector_outer_work_scale
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_outer_work_scale
  lines:
  - 465
  - 465
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
  type: Float64
  units: n/a
  description: Return value of `_effector_outer_work_scale`.
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

# _effector_outer_work_scale

## Purpose
Multiplier applied to the per-worker work threshold when an outer parallel layer is active, raising the bar for inner threading because the pool is already busy.

## Design & Implementation
Returns `_parse_positive_float_env("SPACEAGORA_EFFECTOR_OUTER_WORK_SCALE", 1.5)`. Used only in `_dynamic_effector_thread_decision` as `target_ns = env.effector_work_ns_per_worker_threshold * (outer_active ? env.effector_outer_work_scale : 1.0)`. Captured into `RhsPlanEnvConfig.effector_outer_work_scale`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_effector_outer_work_scale`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:864-864`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:466-466`
- `callees` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:466-466`
<!-- vulcan:connections:end -->

## Limitations
Values below 1.0 are accepted and would make inner threading easier under outer parallelism, the opposite of the intent, with no warning. The scale is only consulted when `allow_with_outer` also permits threading; otherwise the policy refuses regardless.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 465.
