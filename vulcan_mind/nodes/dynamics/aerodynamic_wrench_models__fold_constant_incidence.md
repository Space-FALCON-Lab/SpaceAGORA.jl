---
id: dynamics.aerodynamic_wrench_models__fold_constant_incidence
label: _fold_constant_incidence
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _fold_constant_incidence
  lines:
  - 241
  - 241
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
  description: Return value of `_fold_constant_incidence`.
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

# _fold_constant_incidence

## Purpose
Folds a signed incidence angle in `(-π, π]` into `[0, π/2]` using the box's x/z symmetry so the constant-CD law never produces negative drag.

## Design & Implementation
Computes `a = abs(alpha_rad)` and returns `min(a, π - a)`. The comment records the motivating bug (PR #86): raw `atan2` incidence of −π/2 gave CD = −0.6, pumping orbital energy into a tumbling spacecraft.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `alpha_rad` | Float64 | n/a | yes | Positional argument `alpha_rad`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_fold_constant_incidence`. |
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
Assumes a body symmetric under flow reversal; asymmetric shapes lose the distinction between front and back faces. Only applied in the `:constant` coefficient branch; the fM model keeps signed angles.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 241.
