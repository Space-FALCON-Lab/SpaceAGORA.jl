---
id: gnc.heat_rate_control__edg_maxwellian_heat_rate
label: _edg_maxwellian_heat_rate
kind: function
source:
  file: src/gnc/control/heat_rate_control.jl
  symbol: _edg_maxwellian_heat_rate
  lines:
  - 108
  - 108
inputs:
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: env
  type: Any
  units: n/a
  required: true
  description: Positional argument `env`.
- id: alpha
  type: Float64
  units: n/a
  required: true
  description: Positional argument `alpha`.
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
  type: Float64
  units: n/a
  description: Return value of `_edg_maxwellian_heat_rate`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _edg_maxwellian_heat_rate

## Purpose
Computes the instantaneous free-molecular heat rate in W/cm^2 for the current environment sample and a commanded angle of attack, reading the accommodation factor and planet constants straight out of the integrator parameter object. It is the form of the heat-rate model used for monitoring and logging during energy-depletion aerobraking.

## Theory & Math
$\dot q = K\left[\left(S^2 + \frac{\gamma}{\gamma-1} - \frac{\gamma+1}{2(\gamma-1)}\frac{T_w}{T_p}\right)\left(e^{-s^2} + \sqrt{\pi}s(1+\mathrm{erf}\,s)\right) - \tfrac12 e^{-s^2}\right]$ with $s = S\sin\alpha$ and $K = 10^{-4}\,\rho\,\tau\,R\,T_p\sqrt{R T_p/(2\pi)}$, giving W/cm^2 for $\rho$ in kg/m^3, $R$ in J/kg/K and $T_p$ in K.

## Design & Implementation
Signature `_edg_maxwellian_heat_rate(p::ODEParams, env, alpha::Float64)::Float64`. It pulls `p.args.environment_model.thermal_model` and reads `thermal_accomodation_factor` through `hasproperty`, defaulting to `1.0` when the thermal model does not expose it. Density `rho`, temperature `T_p` and molecular speed ratio `S` come from the `env` sample; `gamma` and `R` from `p.args.environment_model.planet`. If any of `S`, `T_p`, `rho` is non-finite or non-positive it returns `0.0` immediately, treating an absent atmosphere as zero heating. Otherwise it evaluates the same Maxwellian expression as `_energy_depletion_heat_rate_calc` with `T_w = T_p`, and returns the value only when it is finite and strictly positive, otherwise `0.0`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `env` | Any | n/a | yes | Positional argument `env`. |
| in | `alpha` | Float64 | n/a | yes | Positional argument `alpha`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_edg_maxwellian_heat_rate`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_command_alpha_bang|_edg_command_alpha!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:223-223`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_rate_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/heat_rate_control.jl:110-110`
<!-- vulcan:connections:end -->

## Limitations
The physics is duplicated from `_energy_depletion_heat_rate_calc` rather than delegated, so the two can drift apart under edit. Clamping a negative or non-finite result to `0.0` hides genuine model breakdown — a very cold wall relative to free stream can drive the bracket negative and that case reports as no heating at all. The `env` argument is untyped, so a sample missing `molecular_speed_ratio` fails only at the property access. The misspelled field name `thermal_accomodation_factor` must be reproduced exactly by any thermal model, or the silent `1.0` default applies.

## Provenance
Mapped from `src/gnc/control/heat_rate_control.jl` line 108.
