---
id: simulation.setup__effector_cost_ns_per_item_default
label: _effector_cost_ns_per_item_default
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_cost_ns_per_item_default
  lines:
  - 417
  - 417
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
  description: Return value of `_effector_cost_ns_per_item_default`.
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

# _effector_cost_ns_per_item_default

## Purpose
Prior estimate, in nanoseconds, of the wall time to evaluate one effector for one satellite, used by the effector threading decision until enough runtime samples have been collected.

## Design & Implementation
Returns `_parse_positive_float_env("SPACEAGORA_EFFECTOR_COST_NS_PER_ITEM_DEFAULT", 2.5e4)`, so the default assumes 25 μs per effector evaluation. `_effector_observed_cost_ns_per_item` returns this value whenever `shared_buffers` is absent, lacks the cost fields, has fewer than `effector_cost_min_samples`, or holds a non-finite or non-positive estimate.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_effector_cost_ns_per_item_default`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.dynamics_rhs__rhs_effector_static_cost_ns|_rhs_effector_static_cost_ns]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:352-352`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:860-860`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:418-418`
- `callees` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:418-418`
<!-- vulcan:connections:end -->

## Limitations
A single scalar prior cannot represent a mix of cheap inverse-square and expensive high-degree harmonics effectors; the decision may be wrong for the first few steps. Zero or negative values throw `ArgumentError`; `Inf` is accepted and would make every workload appear heavy.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 417.
