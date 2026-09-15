---
id: vehx.thermal_models_getheatrate
label: getHeatRate
kind: function
source:
  file: src/vehicle/thermal/thermal_models.jl
  symbol: getHeatRate
  lines:
  - 93
  - 115
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Vehicle package namespace that loads this file and exports the symbol.
- id: flow_state
  type: Float64
  units: n/a
  required: true
  description: Speed ratio S, free-stream temperature T, density rho, speed v and
    incidence angle alpha.
- id: model
  type: MaxwellianHeat
  units: n/a
  required: true
  description: Thermal model carrying the accommodation factor, planet gas constants
    and contact flag.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: heat_rate
  type: Float64
  units: W/cm^2
  description: Free-molecular heat flux delivered to the surface at the given incidence
    angle.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
- thermal
charts:
- vehx
origin: agent
---

# getHeatRate

## Purpose
`getHeatRate` evaluates the free-molecular aerothermal heating of a spacecraft surface. In very low Earth orbit and during aerobraking the residual atmosphere delivers a measurable flux to ram-facing panels, and this function converts the local flow state into the heat rate that thermal analysis and material limits are checked against.

## Theory & Math
With speed ratio $S = v/\sqrt{2RT}$ and incidence $\alpha$, the recovery temperature is $T_r = T + \frac{\gamma}{\gamma+1} r' (T_0 - T)$ where $T_0 = T\left(1 + \frac{\gamma-1}{\gamma}S^2\right)$. The flux evaluated is $q = \alpha_T \rho R T_p \sqrt{\frac{R T_p}{2\pi}}\left[\left(S^2 + \frac{\gamma}{\gamma-1} - \frac{\gamma+1}{2(\gamma-1)}\frac{T_w}{T_p}\right)\left(e^{-(S\sin\alpha)^2} + \sqrt{\pi} S\sin\alpha\,(1 + \mathrm{erf}(S\sin\alpha))\right) - \tfrac{1}{2}e^{-(S\sin\alpha)^2}\right]$.

## Model & Assumptions
The gas is treated as a Maxwellian free stream in the collisionless regime, so molecules strike the surface without an intervening shock layer and the flux follows from the drifting Maxwellian distribution rather than a continuum boundary layer. The molecular speed ratio sets the balance between thermal and bulk motion, and the incidence angle projects the drift onto the surface normal. A thermal accommodation factor scales the fraction of incident energy actually retained on reflection. The `thermal_contact` flag selects between two closely related recovery formulations, altering how the error function term enters the recovery factor and the Stanton number.

## Design & Implementation
The method dispatches on `MaxwellianHeat` and returns a bare `Float64`, with the return type annotated so the caller sees a concrete type. It first computes the recovery factor and Stanton-like coefficient in the branch on `model.thermal_contact`, then reads the ratio of specific heats and gas constant from the planet object carried by the model, forms the stagnation and recovery temperatures, and evaluates the flux expression. `erf` comes from `SpecialFunctions`, imported at the top of the file. The trailing factor of 1e-4 converts from watts per square metre to watts per square centimetre, which is the unit the returned value carries.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Vehicle package namespace that loads this file and exports the symbol. |
| in | `flow_state` | Float64 | n/a | yes | Speed ratio S, free-stream temperature T, density rho, speed v and incidence angle alpha. |
| in | `model` | MaxwellianHeat | n/a | yes | Thermal model carrying the accommodation factor, planet gas constants and contact flag. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `heat_rate` | Float64 | W/cm^2 | — | Free-molecular heat flux delivered to the surface at the given incidence angle. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/thermal/thermal_models.jl`
- [[simulation.setup__validate_thermal_model_support_bang|_validate_thermal_model_support!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:60-60`
- [[simulation.thermal_callbacks__compute_stage_heat_rates_bang|_compute_stage_heat_rates!]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:54-54`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Wall temperature is set equal to the free-stream-derived plate temperature rather than being solved from a surface energy balance, so radiative equilibrium and conduction into the structure are not represented. The formulation is singular as the speed ratio approaches zero because several terms divide by it, and it is only valid in the free-molecular limit, breaking down as the Knudsen number falls. Surface curvature, shadowing and re-radiation between facets are absent, and the returned unit differs from the watts per square metre used by the convective and radiative helpers earlier in the same file.

## Provenance
Mapped from `src/vehicle/thermal/thermal_models.jl:93-115`, with `MaxwellianHeat` at line 5 and the `heatrate_convective`, `heatrate_radiative` and `heatrate_convective_radiative` helpers preceding it.
