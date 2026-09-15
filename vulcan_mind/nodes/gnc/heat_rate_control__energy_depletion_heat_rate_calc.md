---
id: gnc.heat_rate_control__energy_depletion_heat_rate_calc
label: _energy_depletion_heat_rate_calc
kind: function
source:
  file: src/gnc/control/heat_rate_control.jl
  symbol: _energy_depletion_heat_rate_calc
  lines:
  - 4
  - 4
inputs:
- id: taf
  type: Float64
  units: n/a
  required: true
  description: Positional argument `taf`.
- id: rho
  type: Float64
  units: n/a
  required: true
  description: Positional argument `rho`.
- id: T_w
  type: Float64
  units: n/a
  required: true
  description: Positional argument `T_w`.
- id: T_p
  type: Float64
  units: n/a
  required: true
  description: Positional argument `T_p`.
- id: R
  type: Float64
  units: n/a
  required: true
  description: Positional argument `R`.
- id: gamma
  type: Float64
  units: n/a
  required: true
  description: Positional argument `gamma`.
- id: S
  type: Float64
  units: n/a
  required: true
  description: Positional argument `S`.
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
  description: Return value of `_energy_depletion_heat_rate_calc`.
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

# _energy_depletion_heat_rate_calc

## Purpose
Evaluates the free-molecular convective heat rate on the vehicle for a given angle of attack, using the Maxwellian accommodation model. It is the scalar physics kernel that the energy-depletion guidance law both limits against and inverts when it needs the attitude that just meets a heat-rate constraint.

## Theory & Math
With $s = S\sin\alpha$, $A = e^{-s^2} + \sqrt{\pi}\,s\,(1+\mathrm{erf}\,s)$ and $K = \rho\,\tau\,R\,T_p\sqrt{R T_p/(2\pi)}\times 10^{-4}$, the heat rate is $\dot q = K\left[\left(S^2 + \frac{\gamma}{\gamma-1} - \frac{\gamma+1}{2(\gamma-1)}\frac{T_w}{T_p}\right)A - \tfrac12 e^{-s^2}\right]$, where $S$ is the molecular speed ratio, $\alpha$ the angle of attack in radians, $\rho$ the free-stream density, $\tau$ the thermal accommodation factor, $T_w$ the wall temperature and $T_p$ the free-stream temperature.

## Design & Implementation
Takes eight positional `Float64` arguments — `taf` (thermal accommodation factor), `rho` (kg/m^3), wall temperature `T_w` and free-stream temperature `T_p` (K), specific gas constant `R` (J/kg/K), ratio of specific heats `gamma`, molecular speed ratio `S`, and `alpha` (rad) — and returns a `Float64` in W/cm^2. The scale factor is `first_term = rho * 1e-4 * taf * R * T_p * sqrt(R*T_p/(2pi))`, where the `1e-4` converts W/m^2 to W/cm^2. With `s_sin = S*sin(alpha)` it forms `term_a = exp(-s_sin^2) + sqrt(pi)*s_sin*(1+erf(s_sin))` using `SpecialFunctions.erf`, scales it by the energy bracket `S^2 + gamma/(gamma-1) - (gamma+1)/(2*(gamma-1))*(T_w/T_p)`, subtracts `0.5*exp(-s_sin^2)` for the re-emitted flux, and multiplies by `first_term`. A non-finite result is replaced by `Inf` so that a constraint check against it fails conservatively rather than propagating `NaN`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `taf` | Float64 | n/a | yes | Positional argument `taf`. |
| in | `rho` | Float64 | n/a | yes | Positional argument `rho`. |
| in | `T_w` | Float64 | n/a | yes | Positional argument `T_w`. |
| in | `T_p` | Float64 | n/a | yes | Positional argument `T_p`. |
| in | `R` | Float64 | n/a | yes | Positional argument `R`. |
| in | `gamma` | Float64 | n/a | yes | Positional argument `gamma`. |
| in | `S` | Float64 | n/a | yes | Positional argument `S`. |
| in | `alpha` | Float64 | n/a | yes | Positional argument `alpha`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_energy_depletion_heat_rate_calc`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_profile_heat_rates|_edg_profile_heat_rates]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:395-395`
- [[gnc.heat_rate_control__energy_depletion_heatrate_root_alpha|_energy_depletion_heatrate_root_alpha]] · `callees` → `callers` · call · `src/gnc/control/heat_rate_control.jl:47-47`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_rate_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The expression is singular as `gamma` approaches 1 — both `gamma/(gamma-1)` and `(gamma+1)/(2*(gamma-1))` blow up — and nothing guards that case. Negative `rho`, `T_p` or `R` are accepted and `sqrt(R*T_p)` will throw a `DomainError` for negative arguments rather than returning a value. The continuum and transitional regimes are outside the free-molecular assumption, so the number is only meaningful at very high Knudsen number. Returning `Inf` for a non-finite input erases the distinction between an overflow and a `NaN` input.

## Provenance
Mapped from `src/gnc/control/heat_rate_control.jl` line 4.
