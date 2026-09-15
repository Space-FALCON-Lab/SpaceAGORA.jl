---
id: gnc.targeting_control__edg_integrated_max_energy_depletion_trajectory
label: _edg_integrated_max_energy_depletion_trajectory
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_integrated_max_energy_depletion_trajectory
  lines:
  - 563
  - 563
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
- id: pos0
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos0`.
- id: vel0
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel0`.
- id: mass
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: times
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `times`.
- id: heat_load_switches
  type: NTuple{2, Float64}
  units: n/a
  required: true
  description: Positional argument `heat_load_switches`.
- id: heat_rate_control
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `heat_rate_control`.
- id: structural_control
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `structural_control`.
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
  description: Return value of `_edg_integrated_max_energy_depletion_trajectory`.
    Returns `gravity + aero` or `alpha` or `(`.
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

# _edg_integrated_max_energy_depletion_trajectory

## Purpose
Integrates a drag passage under the maximum-energy-depletion profile, which holds maximum drag except inside a heat-load-limited low-drag window.

## Design & Implementation
Same structure as the targeting integrator, but the per-step angle comes from a local `max_energy_alpha` closure that returns the minimum angle inside the `heat_load_switches` window when the heat-load sub-mode is active, and otherwise constrains the maximum angle. The same RK4 and post-pass sampling follow.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | AerobrakingEnergyDepletionConfig | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `pos0` | SVector{3, Float64} | n/a | yes | Positional argument `pos0`. |
| in | `vel0` | SVector{3, Float64} | n/a | yes | Positional argument `vel0`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `times` | Vector{Float64} | n/a | yes | Positional argument `times`. |
| in | `heat_load_switches` | NTuple{2, Float64} | n/a | yes | Positional argument `heat_load_switches`. |
| in | `heat_rate_control` | Bool | n/a | yes | Keyword argument `heat_rate_control`. |
| in | `structural_control` | Bool | n/a | yes | Keyword argument `structural_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_integrated_max_energy_depletion_trajectory`. Returns `gravity + aero` or `alpha` or `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_predict_max_energy_depletion_outcome|_edg_predict_max_energy_depletion_outcome]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:742-742`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.heat_load_control_acceleration|acceleration]] · `callers` · call · `src/gnc/control/targeting_control.jl:585-585`
- `callees` → [[gnc.targeting_control__edg_targeting_aero_acceleration|_edg_targeting_aero_acceleration]] · `callers` · call · `src/gnc/control/targeting_control.jl:588-588`
- `callees` → [[gnc.targeting_control_acceleration|acceleration]] · `callers` · call · `src/gnc/control/targeting_control.jl:585-585`
<!-- vulcan:connections:end -->

## Limitations
Shares the fixed-angle-per-step and double-sampling costs of the targeting integrator; the two functions duplicate roughly sixty lines and must be edited together.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 563.
