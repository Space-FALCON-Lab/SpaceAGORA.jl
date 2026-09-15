---
id: simulation.setup__rhs_flat_packet_overhead_min_samples
label: _rhs_flat_packet_overhead_min_samples
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_packet_overhead_min_samples
  lines:
  - 453
  - 453
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
  type: Int
  units: n/a
  description: Return value of `_rhs_flat_packet_overhead_min_samples`.
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

# _rhs_flat_packet_overhead_min_samples

## Purpose
Number of packet-scheduled RHS evaluations that must be timed before the overhead ratio is trusted to disable packet scheduling.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_RHS_FLAT_PACKET_OVERHEAD_MIN_SAMPLES", 4)`, at least 1. Paired with `_rhs_flat_packet_overhead_disable_ratio` in `RhsPlanEnvConfig`; the flat-queue telemetry counter must reach this value before the ratio check fires.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_rhs_flat_packet_overhead_min_samples`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:882-882`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:454-454`
<!-- vulcan:connections:end -->

## Limitations
Four samples is small enough that a single cold-start outlier can trip the disable permanently. No upper bound; a very large value makes the overhead guard effectively inert.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 453.
