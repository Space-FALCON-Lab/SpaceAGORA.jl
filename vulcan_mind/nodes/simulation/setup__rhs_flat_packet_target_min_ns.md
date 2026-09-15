---
id: simulation.setup__rhs_flat_packet_target_min_ns
label: _rhs_flat_packet_target_min_ns
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_packet_target_min_ns
  lines:
  - 429
  - 429
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
  description: Return value of `_rhs_flat_packet_target_min_ns`.
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

# _rhs_flat_packet_target_min_ns

## Purpose
Minimum amount of work, in nanoseconds, that the flat-constellation effector queue tries to place in each scheduling packet so that per-packet dispatch overhead is amortised.

## Design & Implementation
Returns `_parse_positive_float_env("SPACEAGORA_RHS_FLAT_PACKET_TARGET_MIN_NS", 2.5e4)`, defaulting to 25 μs. Captured into `RhsPlanEnvConfig` and used by the flat-queue packetiser to group `(satellite, effector)` items whose summed estimated cost reaches the target.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_rhs_flat_packet_target_min_ns`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:876-876`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:430-430`
- `callees` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:430-430`
<!-- vulcan:connections:end -->

## Limitations
Interacts with `_rhs_flat_packet_min_items` and the per-effector cost estimates; if those estimates are stale the packet sizes are wrong in the same direction. There is no maximum, so a huge target collapses the queue into one packet per worker.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 429.
