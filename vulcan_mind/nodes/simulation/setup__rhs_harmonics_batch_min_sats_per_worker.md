---
id: simulation.setup__rhs_harmonics_batch_min_sats_per_worker
label: _rhs_harmonics_batch_min_sats_per_worker
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_harmonics_batch_min_sats_per_worker
  lines:
  - 818
  - 818
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
  description: Return value of `_rhs_harmonics_batch_min_sats_per_worker`.
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

# _rhs_harmonics_batch_min_sats_per_worker

## Purpose
Reads how many satellites each harmonics batch worker should own before an additional worker is spawned, keeping the per-worker coefficient-sharing benefit above the fork cost.

## Design & Implementation
Parses `SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER` through `parse_thread_threshold_env` with a default of 4. Declared `@inline`. The planner divides active satellites by this floor to compute viable worker count.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_rhs_harmonics_batch_min_sats_per_worker`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.rhs_calibration__rhs_plan_candidates|_rhs_plan_candidates]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:224-224`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:872-872`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:819-819`
<!-- vulcan:connections:end -->

## Limitations
Ignored entirely when `harmonics_batch_spin_barrier` is on, where the floor becomes one satellite per worker.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 818.
