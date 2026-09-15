---
id: dynamics.aerodynamic_wrench_models_aerodynamiccoefficientfm
label: AerodynamicCoefficientfM
kind: struct
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: AerodynamicCoefficientfM
  lines:
  - 97
  - 97
inputs:
- id: per_link_atmosphere
  type: Bool
  units: n/a
  required: false
  description: Field `per_link_atmosphere` (default `false`).
- id: fixed_attitude_incidence
  type: Symbol
  units: n/a
  required: false
  description: Field `fixed_attitude_incidence` (default `:max_drag`).
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
  type: AerodynamicCoefficientfM
  units: n/a
  description: Constructed `AerodynamicCoefficientfM` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# AerodynamicCoefficientfM

## Purpose
Effector type selecting the Hart et al. free-molecular rectangular-prism aerodynamics, evaluated per link and summed, with configurable fixed-attitude incidence handling.

## Design & Implementation
`@kwdef struct` with `per_link_atmosphere::Bool=false` (instance-scoped opt-in for per-link density sampling) and `fixed_attitude_incidence::Symbol=:max_drag`. The incidence modes, used only when `orientation_sim=false`, are `:max_drag` (every link flow-normal at full `ref_area`), `:attitude` (incidence from the stored quaternions relative to a flow-aligned frame, zero sideslip), and `:tumbling_average` (normal incidence on Cauchy mean projected area, surface area / 4). The `wrench` methods call `_aero_pure_wrench(:fm, ...)`; `calcForceTorque` has a threaded per-link implementation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `per_link_atmosphere` | Bool | n/a | no | Field `per_link_atmosphere` (default `false`). |
| in | `fixed_attitude_incidence` | Symbol | n/a | no | Field `fixed_attitude_incidence` (default `:max_drag`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AerodynamicCoefficientfM | n/a | — | Constructed `AerodynamicCoefficientfM` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:165-165`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The docstring itself records two known convention issues: `bus_ram_face=:legacy` reference area mismatches the Hart normalisation face, and `planet.R`/`planet.γ` are sea-level air values that overstate the molecular speed ratio at exospheric altitude. `:attitude` mode cannot represent yaw-only attitudes (degenerates to `:max_drag`). The symbol is validated at evaluation time, not construction.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 97.
