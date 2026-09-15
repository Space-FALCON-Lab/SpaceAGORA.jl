---
id: gnc.targeting_control__edg_targeting_aero_acceleration
label: _edg_targeting_aero_acceleration
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_targeting_aero_acceleration
  lines:
  - 415
  - 415
inputs:
- id: config
  type: AerobrakingEnergyDepletionConfig
  units: n/a
  required: true
  description: Positional argument `config`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: r
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `r`.
- id: v
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `v`.
- id: mass
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass`.
- id: t_abs
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t_abs`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_edg_targeting_aero_acceleration`.
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

# _edg_targeting_aero_acceleration

## Purpose
Computes the inertial aerodynamic acceleration on the vehicle at a predicted state and angle, summing free-molecular lift and drag over every link.

## Theory & Math
$$
\vec{a}_{aero} = \frac{1}{m}\, L_{PI}^\top \sum_{k} q\, A_k \left( C_{D,k}\, \hat{d} + C_{L,k}\, \hat{l} \right),\qquad \hat{d} = -\hat{v}_{rw},\quad \hat{l} = \hat{h} \times \hat{v}_{rw}
$$

## Design & Implementation
Samples the prediction environment and returns zero if density or speed vanish. It forms the drag direction opposite airspeed and the lift direction from the orbit normal crossed with the airspeed direction. For each link with positive area it uses `π/2` for the root, the commanded `alpha` for controlled links and the link's own clamped `α` otherwise, evaluates `aerodynamic_coefficient_fM` at the temperature, speed ratio and link angles, and accumulates `q A (C_D drag + C_L lift)` in the planet frame. The sum is rotated to inertial through `l_pi'` and divided by mass.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | AerobrakingEnergyDepletionConfig | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `r` | SVector{3, Float64} | n/a | yes | Positional argument `r`. |
| in | `v` | SVector{3, Float64} | n/a | yes | Positional argument `v`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `t_abs` | Float64 | n/a | yes | Positional argument `t_abs`. |
| in | `alpha` | Float64 | n/a | yes | Positional argument `alpha`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_edg_targeting_aero_acceleration`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_integrated_max_energy_depletion_trajectory|_edg_integrated_max_energy_depletion_trajectory]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:588-588`
- [[gnc.targeting_control_acceleration|acceleration]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:482-482`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/targeting_control.jl:444-444`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/control/targeting_control.jl:448-448`
- `callees` → [[gnc.targeting_control__edg_targeting_prediction_environment|_edg_targeting_prediction_environment]] · `callers` · call · `src/gnc/control/targeting_control.jl:426-426`
<!-- vulcan:connections:end -->

## Limitations
Lift is confined to the orbit plane by construction, so out-of-plane forces from asymmetric panel angles are not represented; the controlled-link set is rebuilt as a `Set` on every call.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 415.
