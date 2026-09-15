---
id: envana.env_gravity_effectors_gravityeffectors
label: GravityEffectors
kind: struct
source:
  file: src/environment/gravity/gravity_effectors.jl
  symbol: GravityEffectors
  lines:
  - 1
  - 9
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: EnvironmentModels namespace under which the gravity effector facade
    is loaded.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: gravity_symbols
  type: Module
  units: n/a
  description: Re-exported gravity model types and the aerobraking and J2 secular
    helpers.
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
# GravityEffectors

## Purpose
`GravityEffectors` is the environment-side facade that republishes the three gravity force/torque models and two gravity analysis functions defined in the dynamics tree, giving environment and analysis code a single import point for gravity symbols.

## Theory & Math
The re-exported models span three fidelity levels. `ConstantGravityModel` applies a fixed acceleration vector in m/s^2. `InverseSquaredGravityModel` applies `a = -mu * r / |r|^3`, with gravitational parameter `mu` in m^3/s^2 and inertial position `r` in metres. `InverseSquaredJ2GravityModel` adds the oblateness term whose radial component is `-1.5 * J2 * mu * R_e^2 / |r|^4 * (1 - 5 (z/|r|)^2)`, where `J2` is the dimensionless second zonal coefficient, `R_e` is the equatorial radius in metres, and `z` is the component along the spin axis. `j2_secular_rates` returns the mean rates of node and perigee drift implied by that same term, in rad/s.

## Model & Assumptions
The facade assumes the dynamics package is loaded first, because both `using` clauses resolve relative parent paths at definition time. It introduces no methods and no state, so all numerical assumptions belong to the underlying implementations: point-mass or oblate-spheroid gravity, inertial frame position input, and mass taken from the vehicle state vector.

## Design & Implementation
Nine lines total: a `module` header, two `using` clauses splitting model types from analysis helpers, two matching `export` lines, and the closing `end`. Keeping the type imports and the function imports on separate lines makes it obvious which symbols are dispatch targets and which are plain computations.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | EnvironmentModels namespace under which the gravity effector facade is loaded. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `gravity_symbols` | Module | n/a | — | Re-exported gravity model types and the aerobraking and J2 secular helpers. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/gravity/gravity_effectors.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the five listed symbols are forwarded, so spherical-harmonic gravity and third-body terms are not reachable through this facade. Because it re-exports rather than wraps, a change in the dynamics-side signature propagates immediately with no compatibility shim.

## Provenance
Read directly from `src/environment/gravity/gravity_effectors.jl:1-9`; both `using` clauses and both `export` lines were inspected.
