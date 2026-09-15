---
id: simulation.setup__effector_work_ns_per_worker_threshold
label: _effector_work_ns_per_worker_threshold
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_work_ns_per_worker_threshold
  lines:
  - 461
  - 461
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
  description: Return value of `_effector_work_ns_per_worker_threshold`.
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

# _effector_work_ns_per_worker_threshold

## Purpose
Nanoseconds of estimated effector work per thread that must be available before the inner effector loop is considered heavy enough to thread in `:auto` mode.

## Design & Implementation
Returns `_parse_positive_float_env("SPACEAGORA_EFFECTOR_WORK_NS_PER_WORKER_THRESHOLD", 4.0e4)`, defaulting to 40 μs. In `_dynamic_effector_thread_decision`, `work_per_worker_ns = per_effector_cost_ns * n_effectors / max_allotment` is compared against `target_ns = threshold × (outer_active ? effector_outer_work_scale : 1.0)` to set `heavy_work`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_effector_work_ns_per_worker_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:863-863`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:462-462`
- `callees` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:462-462`
<!-- vulcan:connections:end -->

## Limitations
The 40 μs figure encodes an assumption about task-spawn overhead on the machines used for tuning (12-core Mac and 4-vCPU CI per source comments); other hardware may need retuning. `Inf` is accepted and disables inner threading in `:auto` mode.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 461.
