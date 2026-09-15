---
id: gnc.heat_load_control__edg_weighted_aero_coefficients
label: _edg_weighted_aero_coefficients
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_weighted_aero_coefficients
  lines:
  - 50
  - 50
inputs:
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: temperature
  type: Float64
  units: n/a
  required: true
  description: Positional argument `temperature`.
- id: speed_ratio
  type: Float64
  units: n/a
  required: true
  description: Positional argument `speed_ratio`.
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
  type: Any
  units: n/a
  description: Return value of `_edg_weighted_aero_coefficients`. Returns `lift_area
    / area, max(drag_area / area, eps(Float64))`.
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

# _edg_weighted_aero_coefficients

## Purpose
Computes area-weighted lift and drag coefficients for the whole spacecraft at a given angle of attack using the free-molecular panel model, giving the heat-load predictor a single effective CL and CD.

## Design & Implementation
Signature `(spacecraft, temperature::Float64, speed_ratio::Float64, alpha::Float64)`. For each link with positive `ref_area` it calls `aerodynamic_coefficient_fM(link, temperature, speed_ratio, alpha, Float64(link.β), Float64(link.θ))`, accumulating `coeffs[1] * link_area` into `lift_area` and `max(0, coeffs[2]) * link_area` into `drag_area`. Both sums are divided by `_edg_total_ref_area(spacecraft)`; the drag coefficient is floored at `eps(Float64)`. Returns the tuple `(CL, CD)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `temperature` | Float64 | n/a | yes | Positional argument `temperature`. |
| in | `speed_ratio` | Float64 | n/a | yes | Positional argument `speed_ratio`. |
| in | `alpha` | Float64 | n/a | yes | Positional argument `alpha`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_weighted_aero_coefficients`. Returns `lift_area / area, max(drag_area / area, eps(Float64))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_heat_load_coefficients|_edg_heat_load_coefficients]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:67-67`
- [[gnc.heat_load_control_acceleration|acceleration]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:204-204`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/heat_load_control.jl:55-55`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/control/heat_load_control.jl:57-57`
- `callees` → [[gnc.heat_load_control__edg_total_ref_area|_edg_total_ref_area]] · `callers` · call · `src/gnc/control/heat_load_control.jl:51-51`
<!-- vulcan:connections:end -->

## Limitations
All links are evaluated at the same `alpha`, so per-panel orientation offsets are only captured through `link.β` and `link.θ`. Negative drag contributions are clipped to zero while negative lift is kept. Each call loops over every link and calls the panel model, which is expensive inside the RK4 integrator used by `_edg_integrated_heat_load_trajectory`.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 50.
