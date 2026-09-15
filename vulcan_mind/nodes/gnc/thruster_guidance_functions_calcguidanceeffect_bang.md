---
id: gnc.thruster_guidance_functions_calcguidanceeffect_bang
label: calcGuidanceEffect!
kind: function
source:
  file: src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl
  symbol: calcGuidanceEffect!
  lines:
  - 111
  - 111
inputs:
- id: guidanceAlg
  type: AerobrakingCampaignPropulsiveManeuverGuidanceModel
  units: n/a
  required: true
  description: Positional argument `guidanceAlg`.
- id: u
  type: ComponentVector
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
  type: Int64
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
  type: Nothing
  units: n/a
  description: Return value of `calcGuidanceEffect!`; mutates `guidanceAlg` in place.
    Returns `nothing`.
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

# calcGuidanceEffect!

## Purpose
Guidance step for propulsive aerobraking maneuvers. Two methods write a `PropulsiveManeuverCommand` into `p.shared_buffers.maneuver_commands[i]`: one replays a precomputed campaign table of per-orbit delta-v (optionally rescaled by flight apoapsis), the other computes a single periapsis-raise burn at apoapsis once the apoapsis has decayed to a target radius.

## Theory & Math
Periapsis-raise burn at apoapsis radius $r_a$ (m): target semi-major axis $a_t=\tfrac{1}{2}(r_a+r_{p,t})$ where $r_{p,t}$ is the target periapsis radius (m), and $\Delta v = \sqrt{\mu\left(\frac{2}{r_a}-\frac{1}{a_t}\right)}-\sqrt{\mu\left(\frac{2}{r_a}-\frac{1}{a}\right)}$ with $\mu$ the gravitational parameter (m^3/s^2) and $a$ the current semi-major axis (m).

## Design & Implementation
Both methods take `(guidanceAlg, u::ComponentVector, p::ODEParams, t::Float64, i::Int64)` and return `nothing` without writing when `i` is outside `maneuver_commands`. The campaign method looks up `p.orbit_counter[i]` in `guidanceAlg.maneuver_orbit_number`; a hit yields `delta_v_cmd = maneuver_Δv[idx] * _flight_apoapsis_ratio_scale(...)`, stored as `abs(delta_v_cmd)` with `direction_rad` 0.0 for positive and π for negative, otherwise a valid command with zero delta-v is stored. The apoapsis-target method is a three-state machine on `guidanceAlg.command_state[i]` (grown via `_ensure_apo_target_state!`): `DISCARDED` clears the command; `COMMAND_ISSUED` clears and transitions to `DISCARDED` once `maneuver_burn_plans[i].valid` shows the controller has locked the plan; `IDLE` computes osculating elements, requires `a(1+e) <= target_apoapsis_radius_m + apoapsis_tolerance_m`, requires the wrapped true anomaly to be pre-apoapsis within `apoapsis_window_rad`, converts `target_periapsis_altitude_m` to a radius along the planet-fixed periapsis direction (via `planet_frame_lpi` at `et_start[] + t` and `_radius_for_oblate_altitude`), and issues a vis-viva delta-v `sqrt(μ(2/r_a - 1/a_target)) - sqrt(μ(2/r_a - 1/a))` with `a_target = (r_a + r_p_target)/2`. A non-positive delta-v discards the request permanently.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `guidanceAlg` | AerobrakingCampaignPropulsiveManeuverGuidanceModel | n/a | yes | Positional argument `guidanceAlg`. |
| in | `u` | ComponentVector | n/a | yes | Positional argument `u`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `calcGuidanceEffect!`; mutates `guidanceAlg` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`
- [[simulation.control_callbacks__run_guidance_for_thruster_schedule_bang|_run_guidance_for_thruster_schedule!]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:30-30`
- [[simulation.navigation_guidance_callbacks_get_guidance_callbacks|get_guidance_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/navigation_guidance_callbacks.jl:27-27`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:124-124`
- `callees` → [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:180-180`
- `callees` → [[gnc.propulsive_maneuver_command|PropulsiveManeuverCommand]] · `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:122-122`
- `callees` → [[gnc.thruster_guidance_functions__ensure_apo_target_state_bang|_ensure_apo_target_state!]] · `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:142-142`
- `callees` → [[gnc.thruster_guidance_functions__flight_apoapsis_ratio_scale|_flight_apoapsis_ratio_scale]] · `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:125-125`
- `callees` → [[gnc.thruster_guidance_functions__radius_for_oblate_altitude|_radius_for_oblate_altitude]] · `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:183-183`
- `callees` → [[gnc.thruster_guidance_functions__wrap_2pi_guidance|_wrap_2pi_guidance]] · `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:173-173`
- `callees` → [[gncz.thruster_guidance_functions__osculating_elements_and_periapsis_direction|_osculating_elements_and_periapsis_direction]] · `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:161-161`
<!-- vulcan:connections:end -->

## Limitations
The apoapsis-target machine fires at most once per spacecraft and never re-arms, so a second periapsis raise requires a fresh model instance. Direction is encoded only as 0 or π (along/anti velocity), and the delta-v is computed impulsively at exact apoapsis even though the command is issued up to `apoapsis_window_rad` early. Both methods mutate `p.shared_buffers.maneuver_commands` and the apoapsis method mutates `guidanceAlg.command_state`; neither is safe if the same `guidanceAlg` is shared across concurrently integrated spacecraft. Non-finite elements or an unbound orbit silently return without writing.

## Provenance
Mapped from `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl` line 111.
