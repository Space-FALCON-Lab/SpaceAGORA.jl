---
id: gnc.heat_load_control__edg_closed_form_heat_load_trajectory
label: _edg_closed_form_heat_load_trajectory
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_closed_form_heat_load_trajectory
  lines:
  - 124
  - 124
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
- id: pos
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos`.
- id: vel
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel`.
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
  type: Any
  units: n/a
  description: Return value of `_edg_closed_form_heat_load_trajectory`. Returns `(time=times,
    h=h, gamma=gamma, speed=speed, rho=rho, temperature=temperature, sp`.
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

# _edg_closed_form_heat_load_trajectory

## Purpose
Generates a fast analytic approximation of the altitude, flight-path angle, speed and atmosphere history over the drag pass, used as the baseline track for the heat-load switch solver.

## Theory & Math
$$h(\tau) = h_0 + v_0\gamma_0\left(\tau - \frac{\tau^2}{2 t_p}\right),\qquad v(\tau) = \sqrt{\mu\left(\frac{2}{R_e + h} - \frac{1}{a}\right)},\qquad S = \sqrt{\tfrac{\gamma}{2}}\,\frac{v}{\sqrt{\gamma R T}}$$ where $h_0$ is the current altitude (m), $\gamma_0$ the initial flight-path angle (rad), $t_p = \Delta t / 2$ the assumed time to periapsis (s), $R_e$ the equatorial radius (m), and $S$ the molecular speed ratio.

## Design & Implementation
Computes `duration` via `_edg_drag_passage_duration`, the grid via `_edg_prediction_time_grid`, and initial `r0`, `v0`, and `gamma0 = asin(dot(pos, vel) / (r0 v0))`. With `a` the osculating semi-major axis (or `r0` when invalid) and `t_peri = duration / 2`, altitude follows the parabola `altitude0 + v0 gamma0 (tau - tau^2 / (2 t_peri))`, speed comes from the vis-viva equation `sqrt(μ (2/r - 1/a))`, and the flight-path angle from `asin(hdot / v)` with `hdot = v0 gamma0 (1 - tau / t_peri)`. Atmosphere is sampled per point with `_edg_sample_prediction_atmosphere`; the molecular speed ratio is `sqrt(γ/2) v / c` with `c = sqrt(γ R T)`. Returns a NamedTuple of seven `Vector{Float64}` fields `time, h, gamma, speed, rho, temperature, speed_ratio`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `env` | Any | n/a | yes | Positional argument `env`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_closed_form_heat_load_trajectory`. Returns `(time=times, h=h, gamma=gamma, speed=speed, rho=rho, temperature=temperature, sp`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_heat_load_profile_for_k|_edg_heat_load_profile_for_k]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:560-560`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_drag_passage_duration|_edg_drag_passage_duration]] · `callers` · call · `src/gnc/control/heat_load_control.jl:134-134`
- `callees` → [[gnc.heat_load_control__edg_prediction_time_grid|_edg_prediction_time_grid]] · `callers` · call · `src/gnc/control/heat_load_control.jl:135-135`
- `callees` → [[gnc.heat_load_control__edg_sample_prediction_atmosphere|_edg_sample_prediction_atmosphere]] · `callers` · call · `src/gnc/control/heat_load_control.jl:160-160`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/control/heat_load_control.jl:148-148`
<!-- vulcan:connections:end -->

## Limitations
Drag is entirely neglected in the trajectory, and periapsis is assumed to occur at exactly half the pass duration, which is wrong whenever the vehicle is already past periapsis. Altitude is measured from `planet.Rp_e` (spherical), and the vis-viva speed uses that radius, introducing errors for oblate planets.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 124.
