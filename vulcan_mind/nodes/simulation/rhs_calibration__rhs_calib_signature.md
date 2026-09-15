---
id: simulation.rhs_calibration__rhs_calib_signature
label: _rhs_calib_signature
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _rhs_calib_signature
  lines:
  - 76
  - 76
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
  type: String
  units: n/a
  description: Return value of `_rhs_calib_signature`.
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

# _rhs_calib_signature

## Purpose
Builds the string key under which a calibration result is stored, capturing every factor that materially changes which execution plan is fastest.

## Design & Implementation
Takes the ODE parameter object `p` and the `dynamic_effectors` vector. It reads `SimulationModel.ParallelPolicy.effective_inner_thread_budget()`, counts `p.is_active` with `count(identity, ...)`, takes `n_eff = length(dynamic_effectors)`, and sets `has_harmonics` only when there is exactly one effector and it `isa SimulationModel.GravitationalHarmonicsModel`. The fields are joined with `|` in the fixed order `v1`, `machine=`, `budget=`, `sats=<bucket>`, `effs=`, `harm=`; the satellite count is coarsened through `_calib_sat_bucket`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `dynamic_effectors` | Any | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_rhs_calib_signature`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_rhs_calibration_calibrate_rhs_plan_if_needed__calibrate_rhs_plan_if_needed_bang|_calibrate_rhs_plan_if_needed!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:323-323`

**Downstream**

- `callees` → [[parallel.env_config_effective_inner_thread_budget|effective_inner_thread_budget]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:77-77`
- `callees` → [[simulation.rhs_calibration__calib_machine_label|_calib_machine_label]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:84-84`
- `callees` → [[simulation.rhs_calibration__calib_sat_bucket|_calib_sat_bucket]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:86-86`
<!-- vulcan:connections:end -->

## Limitations
The signature ignores effector types beyond the single-harmonics case, harmonics degree and order, the spacecraft state dimension, and integrator settings, so two workloads with the same counts but very different per-call cost share one plan. `count(identity, p.is_active)` assumes `p.is_active` is a boolean collection. The leading `v1` token is the only versioning hook; changing bucket boundaries without bumping it leaves stale entries live.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 76.
