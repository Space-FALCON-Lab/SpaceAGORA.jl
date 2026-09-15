---
id: simulation.setup__rhs_flat_packet_scheduler_mode
label: _rhs_flat_packet_scheduler_mode
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_packet_scheduler_mode
  lines:
  - 433
  - 433
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
  type: Symbol
  units: n/a
  description: Return value of `_rhs_flat_packet_scheduler_mode`.
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

# _rhs_flat_packet_scheduler_mode

## Purpose
Chooses whether the flat-constellation effector queue uses cost-aware packet scheduling (on), plain static striding (off), or decides per step from measured heterogeneity and overhead (auto).

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_RHS_FLAT_PACKET_SCHEDULER"; default="auto")` as `:off`, `:on`, or `:auto`. In `:auto` the packet path is used only when total work exceeds `_rhs_flat_packet_work_ns_threshold`, cost heterogeneity exceeds `_rhs_flat_packet_heterogeneity_threshold`, and the observed overhead ratio has not tripped `_rhs_flat_packet_overhead_disable_ratio`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_rhs_flat_packet_scheduler_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:877-877`

**Downstream**

- `callees` → [[parallel.env_config_parse_parallel_mode_env|parse_parallel_mode_env]] · `callers` · call · `src/simulation/engine/setup.jl:434-434`
<!-- vulcan:connections:end -->

## Limitations
Reuses the generic parallel-mode parser, so the accepted synonyms (`threads`, `serial`) read oddly for a scheduler choice. Only consulted when the execution plan already selected the flat queue; it has no effect in `satellite_batch` or serial modes.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 433.
