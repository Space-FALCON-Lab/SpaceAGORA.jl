---
id: dynamics.perturbations_nbodygravitymodel
label: NBodyGravityModel
kind: struct
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: NBodyGravityModel
  lines:
  - 86
  - 86
inputs:
- id: body_names
  type: NP
  units: n/a
  required: true
  description: Field `body_names`.
- id: body_mus
  type: NM
  units: n/a
  required: true
  description: Field `body_mus`.
- id: primary_body_name
  type: String
  units: n/a
  required: true
  description: Field `primary_body_name`.
- id: planet
  type: P
  units: n/a
  required: true
  description: Field `planet`.
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
  type: NBodyGravityModel
  units: n/a
  description: Constructed `NBodyGravityModel`.
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

# NBodyGravityModel

## Purpose
The third-body gravity effector, listing perturbing bodies by name with their gravitational parameters resolved at construction.

## Design & Implementation
Immutable, parameterised on planet type and on the tuple types of `body_names` and `body_mus`, with `primary_body_name` and the planet. The keyword constructor resolves each GM through `_resolve_third_body_mu`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `body_names` | NP | n/a | yes | Field `body_names`. |
| in | `body_mus` | NM | n/a | yes | Field `body_mus`. |
| in | `primary_body_name` | String | n/a | yes | Field `primary_body_name`. |
| in | `planet` | P | n/a | yes | Field `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NBodyGravityModel | n/a | — | Constructed `NBodyGravityModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:149-149`
- [[dynamics.perturbations_solarradiationpressuremodel|SolarRadiationPressureModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:629-629`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Tuple-typed body lists make the type depend on the body count, so each configuration is a fresh specialisation.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 86.
