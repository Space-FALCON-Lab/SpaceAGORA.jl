---
id: simulation.setup__rhs_flat_work_ns_threshold
label: _rhs_flat_work_ns_threshold
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_work_ns_threshold
  lines:
  - 784
  - 784
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
  description: Return value of `_rhs_flat_work_ns_threshold`.
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

# _rhs_flat_work_ns_threshold

## Purpose
Reads the estimated total effector work per RHS call, in nanoseconds, below which the flat queue is not worth its overhead.

## Design & Implementation
Parses `SPACEAGORA_EFFECTOR_FLAT_WORK_NS_THRESHOLD` as a positive float with default 200,000 ns. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_rhs_flat_work_ns_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:867-867`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:785-785`
- `callees` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:785-785`
<!-- vulcan:connections:end -->

## Limitations
The estimate it is compared against comes from the EMA of observed per-item cost, which is `NaN` until enough samples exist and then falls back to a default.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 784.
