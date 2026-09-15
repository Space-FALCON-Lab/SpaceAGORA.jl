---
id: envana.env_aerodynamic_effectors_aerodynamiceffectors
label: AerodynamicEffectors
kind: struct
source:
  file: src/environment/aerodynamics/aerodynamic_effectors.jl
  symbol: AerodynamicEffectors
  lines:
  - 1
  - 7
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: EnvironmentModels namespace under which the aerodynamic effector facade
    is loaded.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: aero_symbols
  type: Function
  units: n/a
  description: Re-exported aerodynamic_coefficient_fM force/torque model symbol.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- envana
origin: agent
---
# AerodynamicEffectors

## Purpose
`AerodynamicEffectors` is the environment-side facade module that republishes the aerodynamic force/torque effector defined in the dynamics tree. It exists so that environment-facing code can reach `aerodynamic_coefficient_fM` without importing the whole `DynamicEffectors` subtree, keeping the dependency direction one-way from environment code into dynamics.

## Theory & Math
The re-exported effector evaluates the aerodynamic wrench from the free-stream dynamic pressure `q = 0.5 * rho * V_rel^2`, in pascals, where `rho` is atmospheric density in kg/m^3 and `V_rel` is the planet-relative speed in m/s. The resulting force magnitude follows `F = q * A_ref * C_F`, with reference area `A_ref` in m^2 and the dimensionless force coefficient `C_F`; the torque follows `T = q * A_ref * L_ref * C_M` with reference length `L_ref` in metres and dimensionless moment coefficient `C_M`. All symbols are evaluated per vehicle at each integrator stage.

## Model & Assumptions
The facade assumes the dynamics package has already been loaded, since the `using ..DynamicEffectors.AerodynamicEffectors` clause resolves at module definition time. It performs no numerical work of its own and adds no methods, so the numerical assumptions of the underlying coefficient model (continuum flow, quasi-steady coefficients, rigid body) carry over unchanged.

## Design & Implementation
The whole module is seven lines: a `module` header, one `using` clause importing `aerodynamic_coefficient_fM`, one `export` statement, and the closing `end`. Because it only forwards a binding, method dispatch and specialisation happen entirely inside the dynamics implementation, so there is zero runtime cost.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | EnvironmentModels namespace under which the aerodynamic effector facade is loaded. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `aero_symbols` | Function | n/a | — | Re-exported aerodynamic_coefficient_fM force/torque model symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/aerodynamics/aerodynamic_effectors.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only `aerodynamic_coefficient_fM` is forwarded; other aerodynamic effectors in the dynamics tree remain invisible through this path. The module cannot be loaded standalone because its `using` clause references a relative parent path. Adding a new effector requires editing both the import and the export line.

## Provenance
Read directly from `src/environment/aerodynamics/aerodynamic_effectors.jl:1-7`; the module body was inspected in full.
