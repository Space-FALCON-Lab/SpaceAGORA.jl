---
id: vehicle.thermal_models_heatrate_convective_radiative
label: heatrate_convective_radiative
kind: function
source:
  file: src/vehicle/thermal/thermal_models.jl
  symbol: heatrate_convective_radiative
  lines:
  - 71
  - 71
inputs:
- id: S
  type: Any
  units: n/a
  required: true
  description: Positional argument `S`.
- id: T
  type: Any
  units: n/a
  required: true
  description: Positional argument `T`.
- id: m
  type: Any
  units: n/a
  required: true
  description: Positional argument `m`.
- id: rho
  type: Any
  units: n/a
  required: true
  description: Positional argument `ρ`.
- id: v
  type: Any
  units: n/a
  required: true
  description: Positional argument `v`.
- id: alpha
  type: Any
  units: n/a
  required: true
  description: Positional argument `α`.
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
  type: Any
  units: n/a
  description: Return value of `heatrate_convective_radiative`. Returns `q_conv`.
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

# heatrate_convective_radiative

## Purpose

`heatrate_convective_radiative(S, T, m, ρ, v, α)` is intended to return the total blunt-body stagnation heat rate as the sum of the convective and radiative correlations defined alongside it, giving callers a single entry point for continuum entry heating without choosing a mechanism.

## Design & Implementation

The body forwards all six arguments unchanged to `heatrate_convective(S, T, m, ρ, v, α)`, binds the result to `q_conv`, and returns `q_conv`. Because the argument list is identical to the two component functions, the three can be swapped at a call site by name alone. The in-body comment block documents the intended return as the sum of convective and radiative flux.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `S` | Any | n/a | yes | Positional argument `S`. |
| in | `T` | Any | n/a | yes | Positional argument `T`. |
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `rho` | Any | n/a | yes | Positional argument `ρ`. |
| in | `v` | Any | n/a | yes | Positional argument `v`. |
| in | `alpha` | Any | n/a | yes | Positional argument `α`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `heatrate_convective_radiative`. Returns `q_conv`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:228-228`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:228-228`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:228-228`

**Downstream**

- `callees` → [[vehicle.thermal_models_heatrate_convective|heatrate_convective]] · `callers` · call · `src/vehicle/thermal/thermal_models.jl:89-89`
<!-- vulcan:connections:end -->

## Limitations

As written the function is a defect: `heatrate_radiative` is never called and the radiative contribution is never added, so the returned value is exactly `heatrate_convective`. For high-velocity entry, where radiative flux is a significant fraction of the total, callers will therefore see an under-prediction. Every limitation of `heatrate_convective` applies transitively, including the hard-coded 0.25 m nose-radius offset, the unused `S`, `T` and `α` arguments, and the W/cm² return units.

## Provenance
Mapped from `src/vehicle/thermal/thermal_models.jl` line 71.
