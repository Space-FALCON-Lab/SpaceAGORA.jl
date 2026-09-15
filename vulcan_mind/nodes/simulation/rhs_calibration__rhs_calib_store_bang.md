---
id: simulation.rhs_calibration__rhs_calib_store_bang
label: _rhs_calib_store!
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _rhs_calib_store!
  lines:
  - 173
  - 173
inputs:
- id: sig
  type: String
  units: n/a
  required: true
  description: Positional argument `sig`.
- id: plan
  type: Any
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: elapsed_mean_ns
  type: Float64
  units: n/a
  required: true
  description: Positional argument `elapsed_mean_ns`.
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
  type: Nothing
  units: n/a
  description: Return value of `_rhs_calib_store!`; mutates `sig` in place.
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

# _rhs_calib_store!

## Purpose
Records a freshly measured winning plan in the in-process calibration cache under its signature, ready for `_rhs_calib_save!` to persist.

## Design & Implementation
Signature `_rhs_calib_store!(sig::String, plan, elapsed_mean_ns::Float64)::Nothing`. It deliberately calls `_rhs_calib_load!()` first so the one-time disk read cannot later clobber this newer in-memory entry, then under `_rhs_calib_lock` writes `_rhs_calib_cache[sig] = Dict("mode" => String(plan.mode), "allotment" => Int(plan.allotment), "elapsed_mean_ns" => elapsed_mean_ns)`. Mutates only the global cache.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sig` | String | n/a | yes | Positional argument `sig`. |
| in | `plan` | Any | n/a | yes | Positional argument `plan`. |
| in | `elapsed_mean_ns` | Float64 | n/a | yes | Positional argument `elapsed_mean_ns`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_rhs_calib_store!`; mutates `sig` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_rhs_calibration_calibrate_rhs_plan_if_needed__calibrate_rhs_plan_if_needed_bang|_calibrate_rhs_plan_if_needed!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:343-343`

**Downstream**

- `callees` → [[simulation.rhs_calibration__rhs_calib_load_bang|_rhs_calib_load!]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:176-176`
<!-- vulcan:connections:end -->

## Limitations
Nothing is persisted here; a caller that forgets `_rhs_calib_save!` loses the result at process exit. `String(plan.mode)` requires `plan.mode` to be a `Symbol` or string; other plan shapes throw `MethodError`. Storing overwrites any existing entry for the signature without comparing elapsed times, so a noisy re-run under `force` mode can replace a better earlier measurement.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 173.
