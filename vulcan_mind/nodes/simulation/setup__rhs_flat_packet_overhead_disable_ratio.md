---
id: simulation.setup__rhs_flat_packet_overhead_disable_ratio
label: _rhs_flat_packet_overhead_disable_ratio
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_packet_overhead_disable_ratio
  lines:
  - 449
  - 449
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
  description: Return value of `_rhs_flat_packet_overhead_disable_ratio`.
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

# _rhs_flat_packet_overhead_disable_ratio

## Purpose
Fraction of measured packet-scheduling overhead relative to useful work above which the `:auto` scheduler mode gives up on packets for the rest of the run.

## Design & Implementation
Returns `_parse_unit_float_env("SPACEAGORA_RHS_FLAT_PACKET_OVERHEAD_DISABLE_RATIO", 0.10)`, constrained to (0, 1]; the default disables packets once overhead reaches 10 % of work. The comparison is only made after `_rhs_flat_packet_overhead_min_samples` observations. Stored in `RhsPlanEnvConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_rhs_flat_packet_overhead_disable_ratio`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:881-881`

**Downstream**

- `callees` → [[simulation.setup__parse_unit_float_env|_parse_unit_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:450-450`
<!-- vulcan:connections:end -->

## Limitations
The disable is sticky for the run; there is no re-enable if conditions improve. A ratio of 1.0 effectively never disables. Overhead measurement itself relies on wall-clock timers that include unrelated system noise.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 449.
