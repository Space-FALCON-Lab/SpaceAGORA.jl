---
id: simulation.setup__rhs_effector_cost_min_samples
label: _rhs_effector_cost_min_samples
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_effector_cost_min_samples
  lines:
  - 457
  - 457
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
  description: Return value of `_rhs_effector_cost_min_samples`.
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

# _rhs_effector_cost_min_samples

## Purpose
Minimum sample count per effector index before `_rhs_effector_observed_cost_ns` returns the measured per-effector cost instead of the caller's fallback.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_RHS_EFFECTOR_COST_MIN_SAMPLES", 4)`, clamped to at least 1. Captured into `RhsPlanEnvConfig.rhs_effector_cost_min_samples` and compared against `shared_buffers.rhs_effector_cost_samples[][eff_idx]`. Distinct from `_effector_cost_min_samples`, which governs the aggregate per-item estimator.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_rhs_effector_cost_min_samples`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:875-875`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:458-458`
<!-- vulcan:connections:end -->

## Limitations
Having two independently tuned minimum-sample knobs for two estimators that share the same α is easy to misconfigure. Because samples are counted per effector index, adding an effector mid-run (not supported) would leave its counter at zero.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 457.
