---
id: core.simulation_configuration_integrationtolerances
label: IntegrationTolerances
kind: struct
source:
  file: src/core/state/simulation_configuration.jl
  symbol: IntegrationTolerances
  lines:
  - 96
  - 96
inputs:
- id: reltol
  type: Float64
  units: n/a
  required: false
  description: Field `reltol` (default `1e-9`).
- id: abstol
  type: Float64
  units: n/a
  required: false
  description: Field `abstol` (default `1e-11`).
- id: reltol_orbit
  type: Float64
  units: n/a
  required: false
  description: Field `reltol_orbit` (default `1e-6`).
- id: abstol_orbit
  type: Float64
  units: n/a
  required: false
  description: Field `abstol_orbit` (default `1e-8`).
- id: reltol_atmosphere
  type: Float64
  units: n/a
  required: false
  description: Field `reltol_atmosphere` (default `1e-7`).
- id: abstol_atmosphere
  type: Float64
  units: n/a
  required: false
  description: Field `abstol_atmosphere` (default `1e-9`).
- id: reltol_quaternion
  type: Float64
  units: n/a
  required: false
  description: Field `reltol_quaternion` (default `1e-9`).
- id: abstol_quaternion
  type: Float64
  units: n/a
  required: false
  description: Field `abstol_quaternion` (default `1e-11`).
- id: reltol_mass
  type: Float64
  units: n/a
  required: false
  description: Field `reltol_mass` (default `1e-8`).
- id: abstol_mass
  type: Float64
  units: n/a
  required: false
  description: Field `abstol_mass` (default `1e-10`).
- id: reltol_heat_load
  type: Float64
  units: n/a
  required: false
  description: Field `reltol_heat_load` (default `1e-7`).
- id: abstol_heat_load
  type: Float64
  units: n/a
  required: false
  description: Field `abstol_heat_load` (default `1e-9`).
- id: reltol_angular_rate
  type: Float64
  units: n/a
  required: false
  description: Field `reltol_angular_rate` (default `1e-8`).
- id: abstol_angular_rate
  type: Float64
  units: n/a
  required: false
  description: Field `abstol_angular_rate` (default `1e-10`).
- id: dt_max
  type: Float64
  units: n/a
  required: false
  description: Field `dt_max` (default `1.0`).
- id: dt_max_orbit
  type: Float64
  units: n/a
  required: false
  description: Field `dt_max_orbit` (default `30.0`).
- id: dt_max_atmosphere
  type: Float64
  units: n/a
  required: false
  description: Field `dt_max_atmosphere` (default `1.0`).
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
  type: IntegrationTolerances
  units: n/a
  description: Constructed `IntegrationTolerances` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# IntegrationTolerances

## Purpose
`IntegrationTolerances` collects the relative and absolute error tolerances per state block and the maximum step sizes the ODE solver may take in each mission phase. It is an optional field of `SimulationConfiguration` and is read by the engine when building solver options.

## Design & Implementation
A `@kwdef struct` of seventeen `Float64` fields. Global defaults are `reltol = 1e-9`, `abstol = 1e-11`; phase-specific pairs are `reltol_orbit/abstol_orbit = 1e-6/1e-8` for the exo-atmospheric arc and `reltol_atmosphere/abstol_atmosphere = 1e-7/1e-9` for drag passages. Per-state-block pairs cover quaternion (`1e-9/1e-11`), mass (`1e-8/1e-10`), heat load (`1e-7/1e-9`) and angular rate (`1e-8/1e-10`). Step caps are `dt_max = 1.0` s, `dt_max_orbit = 30.0` s and `dt_max_atmosphere = 1.0` s. All values are plain fields with no coupling logic.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `reltol` | Float64 | n/a | no | Field `reltol` (default `1e-9`). |
| in | `abstol` | Float64 | n/a | no | Field `abstol` (default `1e-11`). |
| in | `reltol_orbit` | Float64 | n/a | no | Field `reltol_orbit` (default `1e-6`). |
| in | `abstol_orbit` | Float64 | n/a | no | Field `abstol_orbit` (default `1e-8`). |
| in | `reltol_atmosphere` | Float64 | n/a | no | Field `reltol_atmosphere` (default `1e-7`). |
| in | `abstol_atmosphere` | Float64 | n/a | no | Field `abstol_atmosphere` (default `1e-9`). |
| in | `reltol_quaternion` | Float64 | n/a | no | Field `reltol_quaternion` (default `1e-9`). |
| in | `abstol_quaternion` | Float64 | n/a | no | Field `abstol_quaternion` (default `1e-11`). |
| in | `reltol_mass` | Float64 | n/a | no | Field `reltol_mass` (default `1e-8`). |
| in | `abstol_mass` | Float64 | n/a | no | Field `abstol_mass` (default `1e-10`). |
| in | `reltol_heat_load` | Float64 | n/a | no | Field `reltol_heat_load` (default `1e-7`). |
| in | `abstol_heat_load` | Float64 | n/a | no | Field `abstol_heat_load` (default `1e-9`). |
| in | `reltol_angular_rate` | Float64 | n/a | no | Field `reltol_angular_rate` (default `1e-8`). |
| in | `abstol_angular_rate` | Float64 | n/a | no | Field `abstol_angular_rate` (default `1e-10`). |
| in | `dt_max` | Float64 | n/a | no | Field `dt_max` (default `1.0`). |
| in | `dt_max_orbit` | Float64 | n/a | no | Field `dt_max_orbit` (default `30.0`). |
| in | `dt_max_atmosphere` | Float64 | n/a | no | Field `dt_max_atmosphere` (default `1.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | IntegrationTolerances | n/a | — | Constructed `IntegrationTolerances` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__with_study_settings|_with_study_settings]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:687-687`
- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:165-165`
- [[parcore.simulation_configuration_simulationconfiguration|SimulationConfiguration]] · `callees` → `callers` · call · `src/core/state/simulation_configuration.jl:245-245`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Values are not validated for positivity, ordering (`abstol` smaller than `reltol`), or finiteness; a zero or negative tolerance is only rejected by the solver. Which fields are honoured depends on the engine's solver mode; fixed-step symplectic modes ignore the tolerance fields entirely. The per-block tolerances assume the engine's state layout and become meaningless if the layout changes.

## Provenance
Mapped from `src/core/state/simulation_configuration.jl` line 96.
