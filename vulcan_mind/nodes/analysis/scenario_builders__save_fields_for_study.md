---
id: analysis.scenario_builders__save_fields_for_study
label: _save_fields_for_study
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _save_fields_for_study
  lines:
  - 698
  - 698
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
  type: AbstractArray
  units: n/a
  description: Return value of `_save_fields_for_study`. Returns `out` or `[`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# _save_fields_for_study

## Purpose
Declares the extra per-satellite state columns the verification study saves, so error tables can read position and velocity by name.

## Design & Implementation
Defines two getter closures over the integrator state that extract each satellite's `pos` and `vel` as `SVector{3,Float64}` into a vector, and returns two `SaveField`s named `position` and `velocity` with `per_satellite=true` and column prefixes `pos` and `vel`, yielding columns such as `sc1_pos_1`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractArray | n/a | — | Return value of `_save_fields_for_study`. Returns `out` or `[`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_simulation_dataframe|_run_simulation_dataframe]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:29-29`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[simulation.save_fields_savefield|SaveField]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:715-715`
<!-- vulcan:connections:end -->

## Limitations
The getters allocate a fresh vector on every save, and only position and velocity are exported; mass and attitude are not available to the error tables through this path.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 698.
