---
id: dynamics.aerodynamic_wrench_models__aero_link_area
label: _aero_link_area
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _aero_link_area
  lines:
  - 277
  - 277
inputs:
- id: body
  type: Any
  units: n/a
  required: true
  description: Positional argument `body`.
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
  type: Float64
  units: n/a
  description: Return value of `_aero_link_area`.
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

# _aero_link_area

## Purpose
Returns the effective per-link reference area (m^2) for fixed-attitude fM forces: the configured `ref_area`, or the Cauchy mean projected area in `:tumbling_average` mode.

## Theory & Math
Cauchy's projection theorem for a convex body: $\bar A = \tfrac{1}{4} A_{surf}$, and for a box with sides $d_1,d_2,d_3$, $\bar A = \tfrac{1}{2}(d_1 d_2 + d_1 d_3 + d_2 d_3)$.

## Design & Implementation
If `incidence === :tumbling_average`, reads `d = body.dims` and returns `0.5 * (d1*d2 + d1*d3 + d2*d3)`, which equals total box surface area `2(d1d2 + d1d3 + d2d3)` divided by 4. Otherwise returns `body.ref_area`. Called only when `orientation_sim=false`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `body` | Any | n/a | yes | Positional argument `body`. |
| in | `incidence` | Symbol | n/a | yes | Positional argument `incidence`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_aero_link_area`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:808-808`
- [[dynx.coupled_aerodynamic_wrench_models_aero_pure_wrench|_aero_pure_wrench]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:420-420`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The Cauchy result is exact only for convex bodies; a multi-link spacecraft with shadowing panels is over-counted. Requires `body.dims` to be a 3-element box even for links defined with a `ref_area` alone.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 277.
