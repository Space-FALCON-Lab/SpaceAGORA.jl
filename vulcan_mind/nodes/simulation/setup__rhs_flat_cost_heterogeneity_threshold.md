---
id: simulation.setup__rhs_flat_cost_heterogeneity_threshold
label: _rhs_flat_cost_heterogeneity_threshold
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_cost_heterogeneity_threshold
  lines:
  - 792
  - 792
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
  description: Return value of `_rhs_flat_cost_heterogeneity_threshold`.
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

# _rhs_flat_cost_heterogeneity_threshold

## Purpose
Reads the ratio of most to least expensive effector above which the effector set counts as heterogeneous enough to benefit from the flat constellation effector queue's dynamic scheduling.

## Design & Implementation
Parses `SPACEAGORA_EFFECTOR_FLAT_COST_HETEROGENEITY_THRESHOLD` as a positive float, defaulting to 3.0. Declared `@inline`. The value is captured into `RhsPlanEnvConfig.flat_cost_heterogeneity_threshold` at setup and consumed by `_rhs_effectors_have_heavy_or_heterogeneous_cost`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_rhs_flat_cost_heterogeneity_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:869-869`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:793-793`
- `callees` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:793-793`
<!-- vulcan:connections:end -->

## Limitations
It is compared against static per-type cost estimates rather than measured costs, so the threshold's meaning is only as good as that table.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 792.
