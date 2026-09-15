---
id: simulation.rhs_calibration__rhs_plan_candidates
label: _rhs_plan_candidates
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _rhs_plan_candidates
  lines:
  - 220
  - 220
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: dynamic_effectors
  type: Any
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
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
  type: Any
  units: n/a
  description: Return value of `_rhs_plan_candidates`. Returns `candidates`.
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

# _rhs_plan_candidates

## Purpose
Enumerates the small set of execution plans worth timing for the current constellation size, thread budget, and effector configuration.

## Design & Implementation
Reads `budget = effective_inner_thread_budget()` and `active_sats = count(identity, p.is_active)`. The floor of satellites per worker is 1 when `harmonics_batch_spin_barrier_enabled()` is true, otherwise `_rhs_harmonics_batch_min_sats_per_worker()`; `viable_workers = fld(active_sats, max(1, floor))`. The list always starts with `_make_calib_satellite_batch_plan()`. Flat plans are added only if `viable_workers >= 2` and `_rhs_flat_supported(dynamic_effectors)`; allotments probed are 2, then `max(2, viable_workers ÷ 2)` and `viable_workers` when more than two workers are viable, plus `budget` when it exceeds `viable_workers`. After `sort!(unique!(...))` each allotment is clamped to `budget` and dropped if below 2. Returns a `Vector{Any}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `dynamic_effectors` | Any | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_rhs_plan_candidates`. Returns `candidates`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- [[simulation.rhs_calibration__run_rhs_sweep_bang|_run_rhs_sweep!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:255-255`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:231-231`
- `callees` → [[parallel.env_config_effective_inner_thread_budget|effective_inner_thread_budget]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:221-221`
- `callees` → [[simulation.rhs_calibration__make_calib_flat_plan|_make_calib_flat_plan]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:243-243`
- `callees` → [[simulation.rhs_calibration__make_calib_satellite_batch_plan|_make_calib_satellite_batch_plan]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:227-227`
- `callees` → [[simulation.setup__rhs_flat_supported|_rhs_flat_supported]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:229-229`
- `callees` → [[simulation.setup__rhs_harmonics_batch_min_sats_per_worker|_rhs_harmonics_batch_min_sats_per_worker]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:224-224`
<!-- vulcan:connections:end -->

## Limitations
At most five candidates are ever generated, so the true optimum may lie between the probed allotments. The result is untyped (`Any[]`), which is harmless here but prevents inference downstream. `_rhs_flat_supported` and `_rhs_harmonics_batch_min_sats_per_worker` are defined in sibling engine files and must be loaded first.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 220.
