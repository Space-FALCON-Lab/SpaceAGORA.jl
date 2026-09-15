---
id: analysis.scenario_builders__harmonics_order
label: _harmonics_order
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _harmonics_order
  lines:
  - 29
  - 29
inputs:
- id: degree
  type: Int
  units: n/a
  required: true
  description: Positional argument `degree`.
- id: order
  type: Int
  units: n/a
  required: true
  description: Positional argument `order`.
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
  description: Return value of `_harmonics_order`.
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

# _harmonics_order

## Purpose
Resolves the tesseral order for a spherical-harmonics gravity model from the manifest's degree and order fields, with sentinel handling.

## Design & Implementation
Returns zero when degree is non-positive or order is exactly zero (zonal-only), returns the degree itself when order is negative (the full-field sentinel), and otherwise the smaller of order and degree so the order can never exceed the degree.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `degree` | Int | n/a | yes | Positional argument `degree`. |
| in | `order` | Int | n/a | yes | Positional argument `order`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_harmonics_order`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:135-135`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A negative order means full field by convention, which is not obvious from the manifest schema; a user writing `-1` to mean disabled gets the opposite.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 29.
