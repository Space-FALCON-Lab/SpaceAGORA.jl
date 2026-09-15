---
id: parallel.env_config_adaptive_rho
label: adaptive_rho
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: adaptive_rho
  lines:
  - 223
  - 223
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
  description: Return value of `adaptive_rho`.
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

# adaptive_rho

## Purpose
Supplies the growth factor ρ that bounds how far the adaptive policy may 'desire' more threads than the pool actually holds, enforcing ρ > 1.

## Design & Implementation
Reads `SPACEAGORA_PARALLEL_POLICY_RHO` through `parse_float_env` with default `1.5`. Throws `ArgumentError("SPACEAGORA_PARALLEL_POLICY_RHO must satisfy ρ > 1, got '<ρ>'")` when `!(ρ > 1.0)`. Returns `Float64`; consumed by `_adaptive_desire_cap` to compute the maximum desired thread count.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `adaptive_rho`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callees` → `callers` · call · `src/parallel/policy/adaptive_decision.jl:45-45`
- [[parcore.env_config_snapshot_policy_decision_env|snapshot_policy_decision_env]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:223-223`
- [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callees` → `callers` · call · `src/parallel/policy/observation_tracking.jl:69-69`

**Downstream**

- `callees` → [[parallel.env_config_parse_float_env|parse_float_env]] · `callers` · call · `src/parallel/policy/env_config.jl:224-224`
<!-- vulcan:connections:end -->

## Limitations
`Inf` satisfies `ρ > 1.0` and would make the desire cap overflow in `ceil(Int, ...)`, raising an `InexactError` later rather than here. `NaN` is rejected only because the comparison is false. No upper bound is enforced.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 223.
