---
id: simulation.setup__effector_cost_ema_alpha
label: _effector_cost_ema_alpha
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_cost_ema_alpha
  lines:
  - 425
  - 425
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
  description: Return value of `_effector_cost_ema_alpha`.
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

# _effector_cost_ema_alpha

## Purpose
Smoothing weight of the exponential moving average that tracks observed effector evaluation cost, balancing responsiveness to load changes against noise from timer jitter.

## Theory & Math
The update is $c_{k+1} = (1-\alpha)\,c_k + \alpha\,s_k$ where $c_k$ is the current cost estimate in ns, $s_k$ the new sample in ns, and $\alpha \in (0,1]$ this weight.

## Design & Implementation
Returns `_parse_unit_float_env("SPACEAGORA_EFFECTOR_COST_EMA_ALPHA", 0.2)`, constrained to (0, 1]. Both `_update_effector_cost_model!` and `_update_rhs_effector_cost_model!` read it via `_rhs_env_config_from_buffers(shared_buffers).effector_cost_ema_alpha` and apply `new = (1 - α) * previous + α * sample`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_effector_cost_ema_alpha`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:862-862`

**Downstream**

- `callees` → [[simulation.setup__parse_unit_float_env|_parse_unit_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:426-426`
<!-- vulcan:connections:end -->

## Limitations
One α governs both the per-item aggregate estimator and the per-effector RHS estimator, though their sample rates differ. Values near 1 make the estimate track the last sample, amplifying jitter; α = 1 exactly is permitted. Values below about 0.05 make the estimator lag for hundreds of steps.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 425.
