---
id: gncx.heat_rate_models_lambdas
label: lambdas
kind: function
source:
  file: src/gnc/guidance/aerobraking/common/heat_rate_models.jl
  symbol: lambdas
  lines:
  - 1
  - 49
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: GuidanceHooks namespace supplying the costate propagator with the mission,
    angle-of-attack profile, guidance gain, and closed-form trajectory histories.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: switching_function
  type: Vector{Float64}
  units: mixed
  description: Lambda switching function, velocity costate history, and initial costate
    triple for the energy-depletion guidance law.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncx
origin: agent
---
# lambdas

## Purpose
`lambdas` propagates the adjoint variables of the aerobraking energy-depletion optimal control problem and evaluates the switching function that decides when the solar panels change angle of attack. It is the mathematical core of the closed-form guidance law: everything the switch solvers do reduces to finding where this function changes sign.

## Theory & Math
With state $(h, \gamma, v)$ and control $\alpha$, the costates satisfy $\dot{\boldsymbol\lambda} = -(\partial f/\partial x)^\top \boldsymbol\lambda - \partial L/\partial x$, propagated backwards from $\lambda_v(t_f) = \nu_E v_f$, $\lambda_\gamma(t_f) = 0$, $\lambda_h(t_f) = \nu_E \mu/(R_p + h_f)^2$. The switching function

$$\lambda_{\text{sw}} = \frac{2 k\, m\, v}{S_{\text{ref}}\, C_{D,\text{slope}}\, \pi}$$

compared against $\lambda_v$ yields the bang-bang angle-of-attack schedule.

## Model & Assumptions
The dynamic state is altitude, flight-path angle, and velocity, so there are three costates: `lambdav`, `lambdag`, and `lambdah`. Vehicle properties are pulled from the structural model, with the reference area from `get_spacecraft_reference_area` and the mass from `get_spacecraft_mass` after traversing the body tree from its root. Aerodynamics enter through the triple of drag-coefficient slope, zero-angle lift coefficient, and zero-angle drag coefficient, matching the linear-in-angle coefficient model used by the control-side routines. Density comes from `density_polyfit` and gravity from the inverse-square law referenced to the planet surface value.

## Design & Implementation
Integration runs backwards. Terminal conditions are set from the terminal multiplier `nu_E`, giving a velocity costate proportional to terminal velocity, a zero flight-path-angle costate, and an altitude costate equal to the terminal gravitational acceleration scaled by the multiplier. The loop then steps from the final index down to the second, forming each costate derivative from the current trajectory sample and applying an explicit Euler step over the local time increment. The switching function itself is algebraic rather than integrated: it is the guidance gain times twice the mass times velocity, divided by the product of reference area, drag-coefficient slope, and pi. The routine returns that switching function, the full velocity costate history, and the initial costate triple taken at the second index.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | GuidanceHooks namespace supplying the costate propagator with the mission, angle-of-attack profile, guidance gain, and closed-form trajectory histories. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `switching_function` | Vector{Float64} | mixed | — | Lambda switching function, velocity costate history, and initial costate triple for the energy-depletion guidance law. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_rate_models_aoa|aoa]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl:56-56`
- [[gnc.switch_window_solver_switch_calculation|switch_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:100-100`
- [[gnc.targeting_solver_control_solarpanels_targeting_closed_form|control_solarpanels_targeting_closed_form]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:327-327`

**Downstream**

- `callees` → [[environment.density_models_density_polyfit|density_polyfit]] · `callers` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl:23-23`
<!-- vulcan:connections:end -->

## Limitations
Backward Euler stepping on a stiff adjoint system is only first-order accurate and inherits whatever grid the closed-form trajectory used. The costate equations are written for the linear-in-angle aerodynamic model and a single controlled surface, so they do not extend to a vehicle whose coefficients vary non-linearly with attitude. The initial costate is read at index two rather than one, so the very first sample of the trajectory is not represented in the returned triple.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/common/heat_rate_models.jl:1-49`.
