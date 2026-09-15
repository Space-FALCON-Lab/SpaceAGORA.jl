---
id: gnc.heat_load_control__edg_drag_passage_duration
label: _edg_drag_passage_duration
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_drag_passage_duration
  lines:
  - 87
  - 87
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
  description: Return value of `_edg_drag_passage_duration`. Returns `config.planning_horizon_s`
    or `min(config.planning_horizon_s, 1_000.0)` or `min(duration, config.planning_horizon_s)`.
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

# _edg_drag_passage_duration

## Purpose
Estimates how many seconds remain until the spacecraft exits the current atmospheric pass by symmetry about periapsis, bounding the prediction horizon for the heat-load solver.

## Theory & Math
$$\Delta t = \frac{(M_{\mathrm{exit}} - M_0) \bmod 2\pi}{n},\qquad n = \sqrt{\mu / a^3}$$ where $M_0$ and $M_{\mathrm{exit}}$ are mean anomalies at the current and mirrored true anomaly, $\mu$ the planet gravitational parameter (m^3/s^2), and $a$ the semi-major axis (m).

## Design & Implementation
Computes osculating elements with `rvtoorbitalelement(pos, vel, mass, planet)` and reads `a = oe[1]`, `e = oe[2]`, `nu = oe[6]`. If the orbit is not a valid ellipse it returns `config.planning_horizon_s`. The exit true anomaly is mirrored as `nu_exit = nu > pi ? 2pi - nu : nu`, mean motion is `n = sqrt(planet.μ / a^3)`, and the duration is `mod(M_exit - M0, 2pi) / n` using `_edg_mean_anomaly_from_true`. Durations that are non-finite or below 1 s return `min(planning_horizon_s, 1000)`; otherwise the result is capped at `planning_horizon_s`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_drag_passage_duration`. Returns `config.planning_horizon_s` or `min(config.planning_horizon_s, 1_000.0)` or `min(duration, config.planning_horizon_s)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_closed_form_heat_load_trajectory|_edg_closed_form_heat_load_trajectory]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:134-134`
- [[gnc.heat_load_control_residual|residual]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:750-750`
- [[gnc.targeting_control__edg_predict_max_energy_depletion_outcome|_edg_predict_max_energy_depletion_outcome]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:725-725`
- [[gnc.targeting_control__edg_predict_targeting_outcome|_edg_predict_targeting_outcome]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:683-683`
- [[gnc.targeting_control__edg_solve_targeting_switch|_edg_solve_targeting_switch]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:967-967`
- [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:750-750`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_mean_anomaly_from_true|_edg_mean_anomaly_from_true]] · `callers` · call · `src/gnc/control/heat_load_control.jl:96-96`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/control/heat_load_control.jl:89-89`
<!-- vulcan:connections:end -->

## Limitations
Assumes the pass is symmetric about periapsis with no drag, so the true exit time is shorter than predicted for a decelerating vehicle. The mirror rule treats any `nu <= pi` as pre-periapsis and returns nearly a full orbit for post-periapsis points close to `pi`. Requires `oe[6]` to be true anomaly in radians in `[0, 2pi)`.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 87.
