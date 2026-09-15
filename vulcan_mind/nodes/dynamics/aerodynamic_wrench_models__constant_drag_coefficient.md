---
id: dynamics.aerodynamic_wrench_models__constant_drag_coefficient
label: _constant_drag_coefficient
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _constant_drag_coefficient
  lines:
  - 231
  - 231
inputs:
- id: alpha_rad
  type: Float64
  units: n/a
  required: true
  description: Positional argument `alpha_rad`.
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
  description: Return value of `_constant_drag_coefficient`.
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

# _constant_drag_coefficient

## Purpose
Linear drag-coefficient law for the constant model: CD rises from 0.8 at zero incidence to 2.2 at 90 degrees.

## Theory & Math
$C_D(\alpha) = 0.8 + \frac{2(2.2 - 0.8)}{\pi}\,\alpha$ for $\alpha \in [0, \pi/2]$ in radians.

## Design & Implementation
Returns `2 * (2.2 - 0.8) / pi * alpha_rad + 0.8`, i.e. slope `2.8/π` per radian. Consumed by `_aero_pure_wrench` after the incidence has been folded by `_fold_constant_incidence`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `alpha_rad` | Float64 | n/a | yes | Positional argument `alpha_rad`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_constant_drag_coefficient`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynx.coupled_aerodynamic_wrench_models_aero_pure_wrench|_aero_pure_wrench]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:461-461`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only physical on `[0, π/2]`; unfolded negative incidence gives negative CD (thrust). The end points 0.8 and 2.2 are hard-coded literals with no configuration hook and no dependence on Knudsen number or geometry.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 231.
