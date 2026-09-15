---
id: simulation.setup__rhs_flat_packet_heterogeneity_threshold
label: _rhs_flat_packet_heterogeneity_threshold
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_packet_heterogeneity_threshold
  lines:
  - 445
  - 445
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
  description: Return value of `_rhs_flat_packet_heterogeneity_threshold`.
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

# _rhs_flat_packet_heterogeneity_threshold

## Purpose
Ratio of most-expensive to least-expensive effector cost above which the flat queue treats the workload as heterogeneous enough to benefit from cost-aware packets rather than uniform striding.

## Design & Implementation
Returns `_parse_positive_float_env("SPACEAGORA_RHS_FLAT_PACKET_HETEROGENEITY_THRESHOLD", 3.0)`. The planner computes `max_cost / min_cost` over per-effector estimates and compares; a ratio at or above 3 (default) means one effector is at least three times as costly as another. Captured into `RhsPlanEnvConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_rhs_flat_packet_heterogeneity_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:880-880`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:446-446`
- `callees` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:446-446`
<!-- vulcan:connections:end -->

## Limitations
A ratio is undefined when the minimum estimate is zero; the estimator guards this by rejecting non-positive estimates, but a fallback prior for one effector and a measured value for another can produce a spurious ratio. Values below 1.0 are accepted and make every workload look heterogeneous.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 445.
