---
id: dynamics.perturbations_lvlhcascadeattitudecontrolmodel
label: LVLHCascadeAttitudeControlModel
kind: struct
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: LVLHCascadeAttitudeControlModel
  lines:
  - 2130
  - 2130
inputs:
- id: q_cmd_lb
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Field `q_cmd_lb`.
- id: k_out
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `k_out`.
- id: w_max
  type: Float64
  units: n/a
  required: true
  description: Field `w_max`.
- id: k_rate
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `k_rate`.
- id: tau_max
  type: Float64
  units: n/a
  required: true
  description: Field `tau_max`.
- id: tau_ff
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `tau_ff`.
- id: q_cmd_lb_2
  type: AbstractVector{<:Real}
  units: n/a
  required: false
  description: Field `q_cmd_lb` (default `SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),`).
- id: k_out_2
  type: AbstractVector{<:Real},
  units: n/a
  required: true
  description: Field `k_out`.
- id: w_max_2
  type: Real,
  units: n/a
  required: true
  description: Field `w_max`.
- id: k_rate_2
  type: AbstractVector{<:Real},
  units: n/a
  required: true
  description: Field `k_rate`.
- id: tau_max_2
  type: Real,
  units: n/a
  required: true
  description: Field `tau_max`.
- id: tau_ff_2
  type: AbstractVector{<:Real}
  units: n/a
  required: false
  description: Field `tau_ff` (default `SVector{3, Float64}(0.0, 0.0, 0.0),`).
- id: q
  type: Any
  units: n/a
  required: false
  description: Field `q` (default `SVector{4, Float64}(q_cmd_lb...)`).
- id: qn
  type: Any
  units: n/a
  required: false
  description: Field `qn` (default `norm(q)`).
- id: ko
  type: Any
  units: n/a
  required: false
  description: Field `ko` (default `SVector{3, Float64}(k_out...)`).
- id: kr
  type: Any
  units: n/a
  required: false
  description: Field `kr` (default `SVector{3, Float64}(k_rate...)`).
- id: tf
  type: Any
  units: n/a
  required: false
  description: Field `tf` (default `SVector{3, Float64}(tau_ff...)`).
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
  type: LVLHCascadeAttitudeControlModel
  units: n/a
  description: Constructed `LVLHCascadeAttitudeControlModel`.
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

# LVLHCascadeAttitudeControlModel

## Purpose
Effector implementing a cascaded LVLH attitude controller as a torque source in the dynamics.

## Design & Implementation
Immutable with the commanded LVLH-to-body quaternion, outer-loop gains `k_out`, rate limit `w_max`, inner-loop gains `k_rate`, torque limit `tau_max` and feedforward `tau_ff`. The constructor requires a finite nonzero quaternion and finite non-negative gains and limits.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q_cmd_lb` | SVector{4, Float64} | n/a | yes | Field `q_cmd_lb`. |
| in | `k_out` | SVector{3, Float64} | n/a | yes | Field `k_out`. |
| in | `w_max` | Float64 | n/a | yes | Field `w_max`. |
| in | `k_rate` | SVector{3, Float64} | n/a | yes | Field `k_rate`. |
| in | `tau_max` | Float64 | n/a | yes | Field `tau_max`. |
| in | `tau_ff` | SVector{3, Float64} | n/a | yes | Field `tau_ff`. |
| in | `q_cmd_lb_2` | AbstractVector{<:Real} | n/a | no | Field `q_cmd_lb` (default `SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),`). |
| in | `k_out_2` | AbstractVector{<:Real}, | n/a | yes | Field `k_out`. |
| in | `w_max_2` | Real, | n/a | yes | Field `w_max`. |
| in | `k_rate_2` | AbstractVector{<:Real}, | n/a | yes | Field `k_rate`. |
| in | `tau_max_2` | Real, | n/a | yes | Field `tau_max`. |
| in | `tau_ff_2` | AbstractVector{<:Real} | n/a | no | Field `tau_ff` (default `SVector{3, Float64}(0.0, 0.0, 0.0),`). |
| in | `q` | Any | n/a | no | Field `q` (default `SVector{4, Float64}(q_cmd_lb...)`). |
| in | `qn` | Any | n/a | no | Field `qn` (default `norm(q)`). |
| in | `ko` | Any | n/a | no | Field `ko` (default `SVector{3, Float64}(k_out...)`). |
| in | `kr` | Any | n/a | no | Field `kr` (default `SVector{3, Float64}(k_rate...)`). |
| in | `tf` | Any | n/a | no | Field `tf` (default `SVector{3, Float64}(tau_ff...)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | LVLHCascadeAttitudeControlModel | n/a | — | Constructed `LVLHCascadeAttitudeControlModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2163-2163`
- `callees` → [[dynamics.perturbations__lvlh_cascade_torque|_lvlh_cascade_torque]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:2168-2168`
<!-- vulcan:connections:end -->

## Limitations
Applies torque directly to the body without modelling an actuator, so it is an idealised controller rather than a reaction-wheel or thruster model.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 2130.
