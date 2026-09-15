---
id: dynamics.aerodynamic_wrench_models__validate_fm_incidence
label: _validate_fm_incidence
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _validate_fm_incidence
  lines:
  - 266
  - 266
inputs:
- id: incidence
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `incidence`.
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
  description: Return value of `_validate_fm_incidence`. Returns `nothing`.
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

# _validate_fm_incidence

## Purpose
Rejects any `fixed_attitude_incidence` symbol other than `:max_drag`, `:attitude`, or `:tumbling_average` with a descriptive `ArgumentError`.

## Design & Implementation
A guard expression using `===` comparisons; on failure throws `ArgumentError("AerodynamicCoefficientfM fixed_attitude_incidence must be ... got :$(incidence).")`. Returns `nothing`. Called at the top of `_aero_pure_wrench` and the fM `calcForceTorque`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `incidence` | Symbol | n/a | yes | Positional argument `incidence`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_validate_fm_incidence`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:717-717`
- [[dynx.coupled_aerodynamic_wrench_models_aero_pure_wrench|_aero_pure_wrench]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:363-363`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Validation happens per evaluation rather than at struct construction, so a bad configuration surfaces only when the first aerodynamic step runs. The check is also applied to the `:constant` path where the value is ignored.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 266.
