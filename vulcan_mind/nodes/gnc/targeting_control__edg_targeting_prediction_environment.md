---
id: gnc.targeting_control__edg_targeting_prediction_environment
label: _edg_targeting_prediction_environment
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_targeting_prediction_environment
  lines:
  - 328
  - 328
inputs:
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
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
  description: Return value of `_edg_targeting_prediction_environment`. Returns `(`.
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

# _edg_targeting_prediction_environment

## Purpose
The prediction-time counterpart of `_edg_environment_state`, additionally returning the planet-fixed state and rotation so the aerodynamic acceleration can be assembled.

## Design & Implementation
Performs the same rotation, geodetic conversion, density sample and wind composition as the control-time version, but at an arbitrary future `t_abs`, and returns `l_pi`, `pos_pp`, `vel_pp` and `vel_pp_rw` alongside altitude, density, temperature, speed, speed ratio and dynamic pressure.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `r` | SVector{3, Float64} | n/a | yes | Positional argument `r`. |
| in | `v` | SVector{3, Float64} | n/a | yes | Positional argument `v`. |
| in | `t_abs` | Float64 | n/a | yes | Positional argument `t_abs`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_targeting_prediction_environment`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_targeting_aero_acceleration|_edg_targeting_aero_acceleration]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:426-426`
- [[gnc.targeting_control__edg_targeting_bracket_outcomes|_edg_targeting_bracket_outcomes]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:831-831`
- [[gnc.targeting_control_acceleration|acceleration]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:491-491`
- [[gnc.targeting_control_max_energy_alpha|max_energy_alpha]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:601-601`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/targeting_control.jl:356-356`
- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/gnc/control/targeting_control.jl:344-344`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/gnc/control/targeting_control.jl:333-333`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/gnc/control/targeting_control.jl:334-334`
- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `src/gnc/control/targeting_control.jl:335-335`
- `callees` → [[gnc.targeting_control__edg_ephemeris_time|_edg_ephemeris_time]] · `callers` · call · `src/gnc/control/targeting_control.jl:331-331`
- `callees` → [[gnc.targeting_control__edg_planet_frame_lpi|_edg_planet_frame_lpi]] · `callers` · call · `src/gnc/control/targeting_control.jl:332-332`
<!-- vulcan:connections:end -->

## Limitations
Each call is a full density sample; the targeting solve evaluates it four times per RK4 step and once more per output sample, so a 20,000-point grid with certification over many candidates can issue millions of density queries.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 328.
