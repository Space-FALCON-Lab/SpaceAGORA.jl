---
id: simulation.rhs_calibration__rhs_calibration_mode
label: _rhs_calibration_mode
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _rhs_calibration_mode
  lines:
  - 25
  - 25
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
  description: Return value of `_rhs_calibration_mode`.
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

# _rhs_calibration_mode

## Purpose
Translates the `SPACEAGORA_RHS_CALIBRATE` environment variable into one of the three calibration policies `:auto`, `:force`, or `:off`.

## Design & Implementation
`@inline` function returning `Symbol`. It lowercases and strips the value from `_engine_env_get("SPACEAGORA_RHS_CALIBRATE", "auto")` and matches literal sets: `"auto"` or empty gives `:auto`; `"force"` or `"always"` gives `:force`; `"off"`, `"false"`, `"0"`, or `"none"` gives `:off`. Any other string throws `ArgumentError` quoting the offending value. `_calibrate_rhs_plan_if_needed!` calls it twice, once to short-circuit on `:off` and once to decide whether to consult the cache.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_rhs_calibration_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_rhs_calibration_calibrate_rhs_plan_if_needed__calibrate_rhs_plan_if_needed_bang|_calibrate_rhs_plan_if_needed!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:316-316`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:26-26`
<!-- vulcan:connections:end -->

## Limitations
The throw happens at solve setup rather than at configuration time, so a typo aborts the run late. Because the variable is read on each call, changing it between the two calls inside `_calibrate_rhs_plan_if_needed!` could yield inconsistent behaviour, though no code path does that in practice.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 25.
