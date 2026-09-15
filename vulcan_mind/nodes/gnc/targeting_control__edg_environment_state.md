---
id: gnc.targeting_control__edg_environment_state
label: _edg_environment_state
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_environment_state
  lines:
  - 64
  - 64
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: i
  type: Int
  units: n/a
  required: true
  description: Positional argument `i`.
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
  description: Return value of `_edg_environment_state`. Returns `(`.
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

# _edg_environment_state

## Purpose
Evaluates the atmospheric conditions the controller needs at the current state: altitude, density, temperature, airspeed, molecular speed ratio and dynamic pressure.

## Theory & Math
$$
a = \sqrt{\gamma R T},\qquad S = \sqrt{\tfrac{\gamma}{2}}\,\frac{V}{a},\qquad q = \tfrac{1}{2}\rho V^2
$$

## Design & Implementation
Rotates the inertial state into the planet frame at the ephemeris time, converts to geodetic coordinates, samples the density model with wind, and rotates the east-north-up wind into the planet frame via the NED basis. Airspeed is the planet-relative velocity plus wind. Sound speed is `sqrt(γ R T)` and the molecular speed ratio `sqrt(γ/2) V / a`. Returns a named tuple with density and temperature floored at zero and machine epsilon.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_environment_state`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:260-260`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/targeting_control.jl:82-82`
- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/gnc/control/targeting_control.jl:73-73`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/gnc/control/targeting_control.jl:70-70`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/gnc/control/targeting_control.jl:71-71`
- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `src/gnc/control/targeting_control.jl:72-72`
- `callees` → [[gnc.targeting_control__edg_control_pos_vel_mass|_edg_control_pos_vel_mass]] · `callers` · call · `src/gnc/control/targeting_control.jl:66-66`
- `callees` → [[gnc.targeting_control__edg_control_sat_state|_edg_control_sat_state]] · `callers` · call · `src/gnc/control/targeting_control.jl:65-65`
- `callees` → [[gnc.targeting_control__edg_ephemeris_time|_edg_ephemeris_time]] · `callers` · call · `src/gnc/control/targeting_control.jl:69-69`
<!-- vulcan:connections:end -->

## Limitations
It calls `getDensity` directly rather than reading the shared buffers, so under native GRAM every control tick pays a full density sample and takes the GRAM lock, in addition to the RHS's own sampling.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 64.
