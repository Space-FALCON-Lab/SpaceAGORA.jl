---
id: parallel.env_config__adaptive_desire_cap
label: _adaptive_desire_cap
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: _adaptive_desire_cap
  lines:
  - 231
  - 231
inputs:
- id: pool_size
  type: Int
  units: n/a
  required: true
  description: Positional argument `pool_size`.
- id: rho
  type: Float64
  units: n/a
  required: true
  description: Positional argument `ρ`.
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
  description: Return value of `_adaptive_desire_cap`.
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

# _adaptive_desire_cap

## Purpose
Computes the maximum thread count the adaptive controller is allowed to request, scaling the pool size by ρ so the controller can probe slightly beyond the pool without running away.

## Theory & Math
The cap is $c = \max\!\left(1, \left\lceil \rho \cdot \max(1, p) \right\rceil\right)$ where $p$ is the thread pool size (`pool_size`) and $\rho > 1$ is the growth factor from `adaptive_rho`.

## Design & Implementation
`_adaptive_desire_cap(pool_size::Int, ρ::Float64)::Int` returns `max(1, ceil(Int, ρ * max(1, pool_size)))`. The inner `max(1, pool_size)` guards against zero or negative pool sizes; the outer `max(1, ...)` guarantees at least one thread even if ρ were tiny. `ceil` rounds up so fractional products never truncate below the pool size when ρ ≥ 1.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pool_size` | Int | n/a | yes | Positional argument `pool_size`. |
| in | `rho` | Float64 | n/a | yes | Positional argument `ρ`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_adaptive_desire_cap`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callees` → `callers` · call · `src/parallel/policy/adaptive_decision.jl:46-46`
- [[parcore.env_config_snapshot_policy_decision_env|snapshot_policy_decision_env]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:231-231`
- [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callees` → `callers` · call · `src/parallel/policy/observation_tracking.jl:74-74`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
With `ρ = Inf` (accepted by `adaptive_rho`) `ceil(Int, Inf)` throws `InexactError`. Large `pool_size × ρ` products can exceed `typemax(Int)` on 32-bit builds. The cap is not clipped to the real pool size, so the caller must still clamp before spawning tasks.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 231.
