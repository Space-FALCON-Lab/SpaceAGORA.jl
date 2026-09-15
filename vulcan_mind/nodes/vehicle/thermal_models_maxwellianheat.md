---
id: vehicle.thermal_models_maxwellianheat
label: MaxwellianHeat
kind: struct
source:
  file: src/vehicle/thermal/thermal_models.jl
  symbol: MaxwellianHeat
  lines:
  - 5
  - 5
inputs:
- id: thermal_accomodation_factor
  type: Float64
  units: n/a
  required: true
  description: Field `thermal_accomodation_factor`.
- id: planet
  type: P
  units: n/a
  required: true
  description: Field `planet`.
- id: thermal_contact
  type: Bool
  units: n/a
  required: false
  description: Field `thermal_contact` (default `false`).
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
  type: MaxwellianHeat
  units: n/a
  description: Constructed `MaxwellianHeat` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# MaxwellianHeat

## Purpose

`MaxwellianHeat{P <: AbstractPlanet}` is the free-molecular aerothermal model of the vehicle thermal stack. It is an `@kwdef` struct subtyping `AbstractThermalModel`, and it carries the three quantities `getHeatRate` needs to evaluate stagnation heating from a drifting Maxwellian gas: the wall accommodation factor, the planet whose gas constant and specific-heat ratio set the thermodynamics, and a flag selecting the wall boundary condition.

## Design & Implementation

Three fields are stored. `thermal_accomodation_factor::Float64` is the dimensionless energy accommodation coefficient multiplying the whole heat-flux expression. `planet::P` is the parameterised `AbstractPlanet` supplying `planet.R` (specific gas constant, J/(kg K)) and `planet.γ` (ratio of specific heats). `thermal_contact::Bool` defaults to `false` and switches `getHeatRate` between two algebraic forms of the recovery factor `r_prime` and Stanton-like term `St_prime`: the non-contact branch uses `erf(S sin α)` alone, the contact branch uses `1 + erf(S sin α)`. Being immutable and concretely typed, the struct dispatches without allocation inside the heating loop.

## Theory & Math

For a drifting Maxwellian at speed ratio $S = v/\sqrt{2 R T}$ (dimensionless), incidence $\alpha$ (rad), density $\rho$ (kg/m$^3$), gas constant $R$ (J/(kg K)) and specific-heat ratio $\gamma$, the wall flux implemented in `getHeatRate` is

$$q = \alpha_T\,\rho R T_p \sqrt{\frac{R T_p}{2\pi}}\left[\left(S^2 + \frac{\gamma}{\gamma-1} - \frac{\gamma+1}{2(\gamma-1)}\frac{T_w}{T_p}\right)\left(e^{-(S\sin\alpha)^2} + \sqrt{\pi}\,S\sin\alpha\,(1+\operatorname{erf}(S\sin\alpha))\right) - \tfrac12 e^{-(S\sin\alpha)^2}\right]$$

where $\alpha_T$ is `thermal_accomodation_factor` (dimensionless), $T_p$ the free-stream/particle temperature (K) and $T_w$ the wall temperature (K). The stagnation temperature used for the recovery calculation is $T_0 = T\left(1 + \frac{\gamma-1}{\gamma}S^2\right)$.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `thermal_accomodation_factor` | Float64 | n/a | yes | Field `thermal_accomodation_factor`. |
| in | `planet` | P | n/a | yes | Field `planet`. |
| in | `thermal_contact` | Bool | n/a | no | Field `thermal_contact` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | MaxwellianHeat | n/a | — | Constructed `MaxwellianHeat` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:156-156`
- [[grp.src_core_state|core/state/]] · `members_out` → `callers` · call · `src/core/state/no_gram_presets.jl:85-85`
- [[parcore.no_gram_presets_make_no_gram_environment|make_no_gram_environment]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:85-85`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The model is free-molecular only, so it is valid in rarefied flow and not in the continuum regime the Sutton-Graves and Tauber-Sutton correlations in the same file target. Wall temperature is not an independent field: `getHeatRate` sets `T_p = T` and `T_w = T_p`, so a radiative-equilibrium or cold-wall condition cannot be expressed. The computed `T_r` recovery temperature is never used. The returned flux is scaled by `1e-4` to W/cm², which differs from the W/m² units the correlation functions in this file document.

## Provenance
Mapped from `src/vehicle/thermal/thermal_models.jl` line 5.
