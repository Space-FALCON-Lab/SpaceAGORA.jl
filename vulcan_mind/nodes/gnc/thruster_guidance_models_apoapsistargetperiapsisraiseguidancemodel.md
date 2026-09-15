---
id: gnc.thruster_guidance_models_apoapsistargetperiapsisraiseguidancemodel
label: ApoapsisTargetPeriapsisRaiseGuidanceModel
kind: struct
source:
  file: src/gnc/guidance/thruster_guidance/thruster_guidance_models.jl
  symbol: ApoapsisTargetPeriapsisRaiseGuidanceModel
  lines:
  - 12
  - 12
inputs:
- id: target_apoapsis_radius_m
  type: Float64
  units: n/a
  required: true
  description: Field `target_apoapsis_radius_m`.
- id: target_periapsis_altitude_m
  type: Float64
  units: n/a
  required: false
  description: Field `target_periapsis_altitude_m` (default `200.0e3`).
- id: apoapsis_tolerance_m
  type: Float64
  units: n/a
  required: false
  description: Field `apoapsis_tolerance_m` (default `0.0`).
- id: apoapsis_window_rad
  type: Float64
  units: n/a
  required: false
  description: Field `apoapsis_window_rad` (default `deg2rad(30.0)`).
- id: command_state
  type: Vector{Int64}
  units: n/a
  required: false
  description: Field `command_state` (default `Int64[]`).
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
  type: ApoapsisTargetPeriapsisRaiseGuidanceModel
  units: n/a
  description: Constructed `ApoapsisTargetPeriapsisRaiseGuidanceModel` (keyword constructor
    via @kwdef).
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

# ApoapsisTargetPeriapsisRaiseGuidanceModel

## Purpose
Configures a guidance model that raises periapsis by burning near apoapsis until a target apoapsis radius is reached, ending an aerobraking campaign.

## Design & Implementation
A `@kwdef mutable struct` subtyping `AbstractGuidanceModel`. `target_apoapsis_radius_m` is required; `target_periapsis_altitude_m` defaults to 200 km and `apoapsis_tolerance_m` to zero. `apoapsis_window_rad` defaults to 30 degrees converted at construction by `deg2rad`, bounding the true-anomaly arc around apoapsis in which a burn is permitted. `command_state` is a mutable integer vector carrying per-satellite command state across calls.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `target_apoapsis_radius_m` | Float64 | n/a | yes | Field `target_apoapsis_radius_m`. |
| in | `target_periapsis_altitude_m` | Float64 | n/a | no | Field `target_periapsis_altitude_m` (default `200.0e3`). |
| in | `apoapsis_tolerance_m` | Float64 | n/a | no | Field `apoapsis_tolerance_m` (default `0.0`). |
| in | `apoapsis_window_rad` | Float64 | n/a | no | Field `apoapsis_window_rad` (default `deg2rad(30.0)`). |
| in | `command_state` | Vector{Int64} | n/a | no | Field `command_state` (default `Int64[]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ApoapsisTargetPeriapsisRaiseGuidanceModel | n/a | — | Constructed `ApoapsisTargetPeriapsisRaiseGuidanceModel` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A zero default `apoapsis_tolerance_m` means the apoapsis target is only satisfied exactly, so termination depends entirely on the window check rather than on a tolerance band.

## Provenance
Mapped from `src/gnc/guidance/thruster_guidance/thruster_guidance_models.jl` line 12.
