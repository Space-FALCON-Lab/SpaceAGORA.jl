---
id: gnc.targeting_control_solve_energy_switch
label: solve_energy_switch
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: solve_energy_switch
  lines:
  - 1037
  - 1037
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
  description: Return value of `solve_energy_switch`. Returns `t_low + clamp(frac,
    0.0, 1.0) * (t_high - t_low)` or `Roots.find_zero(energy_residual, (t_low, t_high),
    Roots.Brent(); rtol=1e-7)`.
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

# solve_energy_switch

## Purpose
Finds the switch time that hits the target energy, by Brent's method when the certified bracket straddles it and by linear interpolation otherwise.

## Design & Implementation
Forms the residuals at the bracket ends. If either is non-finite or they share a sign it interpolates linearly between `t_low` and `t_high` by the energy fraction, clamped to the bracket. Otherwise it calls `Roots.find_zero` with `Roots.Brent()` at relative tolerance 1e-7.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Roots.find_zero | n/a | — | Return value of `solve_energy_switch`. Returns `t_low + clamp(frac, 0.0, 1.0) * (t_high - t_low)` or `Roots.find_zero(energy_residual, (t_low, t_high), Roots.Brent(); rtol=1e-7)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control_solve_apoapsis_switch|solve_apoapsis_switch]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:1049-1049`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The interpolation fallback assumes energy is linear in switch time, which it is not; it is a reasonable guess rather than a solution when the bracket does not contain the target.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 1037.
