---
id: simulation.rhs_calibration__rhs_calibrate_n_warmup
label: _rhs_calibrate_n_warmup
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _rhs_calibrate_n_warmup
  lines:
  - 35
  - 35
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
  description: Return value of `_rhs_calibrate_n_warmup`.
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

# _rhs_calibrate_n_warmup

## Purpose
Reads how many untimed warm-up `spacecraft_dynamics!` calls precede the timed block for each candidate plan, absorbing JIT compilation and cache-population effects.

## Design & Implementation
`@inline` accessor that parses `SPACEAGORA_RHS_CALIBRATE_N_WARMUP` through `_engine_env_get` with default `"5"`, guarding `parse(Int, strip(...))` in a `try`/`catch` that yields 5 on failure, and finally applying `max(1, n)` so at least one warm-up call always runs. The returned count drives the first inner loop of `_run_rhs_sweep!`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_rhs_calibrate_n_warmup`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- [[simulation.rhs_calibration__run_rhs_sweep_bang|_run_rhs_sweep!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:253-253`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:36-36`
<!-- vulcan:connections:end -->

## Limitations
Warm-up cannot be disabled (minimum is 1), and parse errors degrade silently to the default. Five calls may be insufficient to fully warm thread pools or SPICE memo caches for very large constellations, biasing the first candidate measured (always satellite_batch) toward slower times.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 35.
