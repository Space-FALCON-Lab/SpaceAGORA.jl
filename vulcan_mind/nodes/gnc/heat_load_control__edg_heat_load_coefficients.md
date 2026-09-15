---
id: gnc.heat_load_control__edg_heat_load_coefficients
label: _edg_heat_load_coefficients
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_heat_load_coefficients
  lines:
  - 64
  - 64
inputs:
- id: config
  type: Any
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
- id: env
  type: Any
  units: n/a
  required: true
  description: Positional argument `env`.
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
  type: Tuple
  units: n/a
  description: Return value of `_edg_heat_load_coefficients`. Returns `(cd_slope=cd_slope,
    cl_low=cl_low, cd_low=cd_low)`.
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

# _edg_heat_load_coefficients

## Purpose
Linearises the drag coefficient in angle of attack between the minimum and maximum alpha bounds so the heat-load costate equations can use `cd = cd_low + alpha * cd_slope`.

## Design & Implementation
Takes `config`, `p::ODEParams`, `spacecraft`, and an `env` record with `temperature` and `molecular_speed_ratio`. It evaluates `_edg_weighted_aero_coefficients` at `config.min_alpha_rad` (keeping `cl_low` and `cd_low`) and at `config.max_alpha_rad` (keeping `cd_high`), then computes `cd_slope = (cd_high - cd_low) / max(max_alpha - min_alpha, eps)`. If the slope is non-finite or its magnitude is below `eps(Float64)`, it substitutes the generic values `cd_slope = (2.2 - 0.8) / (pi/2)` and `cd_low = 0.8`. Returns a NamedTuple `(cd_slope, cl_low, cd_low)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `env` | Any | n/a | yes | Positional argument `env`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_edg_heat_load_coefficients`. Returns `(cd_slope=cd_slope, cl_low=cl_low, cd_low=cd_low)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:633-633`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_weighted_aero_coefficients|_edg_weighted_aero_coefficients]] · `callers` · call · `src/gnc/control/heat_load_control.jl:67-67`
<!-- vulcan:connections:end -->

## Limitations
The linear fit uses only two endpoints and ignores the actual (roughly sinusoidal) CD variation between them. The fallback constants 0.8 and 2.2 are hard-coded plate-like values that may not match the vehicle. `cl_low` at minimum alpha is treated as constant across the profile.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 64.
