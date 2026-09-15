---
id: gnc.targeting_control_solve_apoapsis_switch
label: solve_apoapsis_switch
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: solve_apoapsis_switch
  lines:
  - 1048
  - 1048
inputs:
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
  type: Roots.find_zero
  units: n/a
  description: Return value of `solve_apoapsis_switch`. Returns `t_low + clamp(frac,
    0.0, 1.0) * (t_high - t_low)` or `solve_energy_switch()` or `Roots.find_zero(apoapsis_residual,
    (t_low, t_high), Roots.Brent(); rtol=1e-7)`.
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

# solve_apoapsis_switch

## Purpose
Finds the switch time that hits the target apoapsis, falling back to the energy solve when no apoapsis target exists or the bracket does not straddle it.

## Design & Implementation
Returns `solve_energy_switch()` immediately if the target apoapsis is not finite and positive. Otherwise it checks the apoapsis residuals at the bracket ends; a non-finite pair or a shared sign leads to linear interpolation on apoapsis if the denominator is usable, else the energy solve. A proper bracket is solved with Brent at 1e-7.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Roots.find_zero | n/a | — | Return value of `solve_apoapsis_switch`. Returns `t_low + clamp(frac, 0.0, 1.0) * (t_high - t_low)` or `solve_energy_switch()` or `Roots.find_zero(apoapsis_residual, (t_low, t_high), Roots.Brent(); rtol=1e-7)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.targeting_control_evaluate_candidate|evaluate_candidate]] · `callers` · call · `src/gnc/control/targeting_control.jl:1065-1065`
- `callees` → [[gnc.targeting_control_solve_energy_switch|solve_energy_switch]] · `callers` · call · `src/gnc/control/targeting_control.jl:1049-1049`
<!-- vulcan:connections:end -->

## Limitations
Because apoapsis is far more sensitive than energy to the final state, the bracket check fails more often here and the energy solve is used as the fallback more than the name suggests.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 1048.
