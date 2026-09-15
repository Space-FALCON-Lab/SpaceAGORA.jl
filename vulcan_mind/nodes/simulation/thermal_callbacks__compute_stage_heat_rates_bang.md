---
id: simulation.thermal_callbacks__compute_stage_heat_rates_bang
label: _compute_stage_heat_rates!
kind: function
source:
  file: src/simulation/callbacks/thermal_callbacks.jl
  symbol: _compute_stage_heat_rates!
  lines:
  - 12
  - 12
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: use_buffered_density
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `use_buffered_density` (default `false`).
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
  description: Return value of `_compute_stage_heat_rates!`; mutates `p` in place.
    Returns `heat_rates`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _compute_stage_heat_rates!

## Purpose
Computes the convective heat rate on every link of spacecraft `sat_idx` at time `t` for the current integrator stage, writing the results into the shared per-satellite heat-rate buffer and returning it.

## Theory & Math
Speed of sound $a=\sqrt{\gamma R T}$, Mach number $M=v/a$ with $v=\lVert \mathbf{v}_{pp}-\mathbf{w}_{pp}\rVert$, and molecular speed ratio $S=M\sqrt{\gamma/2}=v/\sqrt{2RT}$. Here $\gamma$ is the ratio of specific heats, $R$ the specific gas constant in J/(kg K), $T$ the atmospheric temperature in kelvin, and $\mathbf{w}_{pp}$ the wind expressed in the planet-fixed frame. Each link's heat rate $\dot q(S,T,\rho,v,\alpha)$ comes from the configured thermal model.

## Design & Implementation
After an early return for spacecraft with no links, it samples the planet-relative frame via `engine.sample_planet_frame` and the atmosphere either from the buffered values (`use_buffered_density=true`) or by a fresh `sample_atmosphere` call with `write_buffers=false`. It converts the GRAM east/north/up wind into the planet frame using the NED basis from `latlongtoNED((alt_m, lat_rad, lon_rad))` as `wN*uN + wE*uE - wU*uD`, subtracts it from `planet_frame.vel_pp` to get the relative velocity, and forms the speed ratio `S = sqrt(planet.γ * 0.5) * mach`. Each link's `getHeatRate(thermal_model, S, T, rho, v, alpha)` result is clamped: non-finite or non-positive values are stored as `0.0`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `use_buffered_density` | Bool | n/a | no | Keyword argument `use_buffered_density` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_compute_stage_heat_rates!`; mutates `p` in place. Returns `heat_rates`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1371-1371`
- [[simulation.dynamics_rhs_spacecraft_dynamics_explicit_remainder_bang|spacecraft_dynamics_explicit_remainder!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2118-2118`
- [[simulation.dynamics_rhs_spacecraft_dynamics_slow_bang|spacecraft_dynamics_slow!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1861-1861`
- [[simulation.save_fields__save_heat_rate|_save_heat_rate]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:132-132`
- [[simulation.thermal_callbacks_update_thermal_sat_bang|update_thermal_sat!]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:62-62`
- [[simulation_a.thermal_callbacks_get_thermal_callback|get_thermal_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:62-62`
- [[simx.engine_dynamics_rhs_spacecraft_dynamics_bang|spacecraft_dynamics!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1749-1749`

**Downstream**

- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:37-37`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:22-22`
- `callees` → [[simulation.thermal_callbacks__heat_rate_buffer_for_sat_bang|_heat_rate_buffer_for_sat!]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:20-20`
- `callees` → [[vehx.thermal_models_getheatrate|getHeatRate]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:54-54`
<!-- vulcan:connections:end -->

## Limitations
When density or temperature is non-finite or non-positive, or when relative speed or sound speed is non-positive, the routine returns a fully zeroed buffer, so a bad atmosphere sample is indistinguishable from a genuinely zero heat load. Links whose angle of attack `α` is non-finite are skipped and left at zero. The returned vector aliases shared buffer storage and is invalidated by the next call for the same spacecraft.

## Provenance
Mapped from `src/simulation/callbacks/thermal_callbacks.jl` line 12.
