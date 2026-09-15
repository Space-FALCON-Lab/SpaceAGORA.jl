---
id: gnc.heat_load_control__edg_sample_prediction_atmosphere
label: _edg_sample_prediction_atmosphere
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_sample_prediction_atmosphere
  lines:
  - 111
  - 111
inputs:
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: altitude
  type: Float64
  units: n/a
  required: true
  description: Positional argument `altitude`.
- id: t_abs
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t_abs`.
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
  description: Return value of `_edg_sample_prediction_atmosphere`. Returns `max(0.0,
    rho), max(temperature, eps(Float64))`.
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

# _edg_sample_prediction_atmosphere

## Purpose
Queries the environment density model for density and temperature at a predicted altitude and absolute time, providing the atmospheric inputs for the heat-load trajectory predictors.

## Design & Implementation
Signature `(p::ODEParams, altitude::Float64, t_abs::Float64)`. Calls `getDensity(density_model, max(0, altitude), 0.0, 0.0, t_abs, wind, p)` with latitude and longitude fixed at zero, discards the third return value, and returns `(max(0, rho), max(temperature, eps(Float64)))` in kg/m^3 and K.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `altitude` | Float64 | n/a | yes | Positional argument `altitude`. |
| in | `t_abs` | Float64 | n/a | yes | Positional argument `t_abs`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_sample_prediction_atmosphere`. Returns `max(0.0, rho), max(temperature, eps(Float64))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_closed_form_heat_load_trajectory|_edg_closed_form_heat_load_trajectory]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:160-160`
- [[gnc.heat_load_control_acceleration|acceleration]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:197-197`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `src/gnc/control/heat_load_control.jl:112-112`
<!-- vulcan:connections:end -->

## Limitations
Latitude and longitude are hard-coded to 0.0, so latitude-dependent GRAM density variations along the predicted ground track are ignored. Negative altitude is clamped to zero rather than flagged. Each call goes through the full density model, which for GRAM-backed models can be the dominant cost of a prediction.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 111.
