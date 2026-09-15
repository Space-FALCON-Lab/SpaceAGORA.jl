---
id: simulation.setup__rhs_flat_packet_work_ns_threshold
label: _rhs_flat_packet_work_ns_threshold
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_packet_work_ns_threshold
  lines:
  - 441
  - 441
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
  description: Return value of `_rhs_flat_packet_work_ns_threshold`.
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

# _rhs_flat_packet_work_ns_threshold

## Purpose
Total estimated work, in nanoseconds per RHS call, above which the flat queue's `:auto` scheduler mode considers cost-aware packet scheduling worthwhile.

## Design & Implementation
Returns `_parse_positive_float_env("SPACEAGORA_RHS_FLAT_PACKET_WORK_NS_THRESHOLD", 5.0e6)`, defaulting to 5 ms. The planner sums per-effector cost estimates from `_rhs_effector_observed_cost_ns` across all items and compares with this value. Stored in `RhsPlanEnvConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_rhs_flat_packet_work_ns_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:879-879`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:442-442`
- `callees` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:442-442`
<!-- vulcan:connections:end -->

## Limitations
Because the sum uses EMA estimates, early-run decisions rely on the `_effector_cost_ns_per_item_default` prior and may misclassify. `Inf` is accepted and disables packet scheduling in `:auto` without an explicit off switch.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 441.
