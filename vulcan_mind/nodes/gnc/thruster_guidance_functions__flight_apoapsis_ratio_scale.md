---
id: gnc.thruster_guidance_functions__flight_apoapsis_ratio_scale
label: _flight_apoapsis_ratio_scale
kind: function
source:
  file: src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl
  symbol: _flight_apoapsis_ratio_scale
  lines:
  - 92
  - 92
inputs:
- id: guidanceAlg
  type: AerobrakingCampaignPropulsiveManeuverGuidanceModel
  units: n/a
  required: true
  description: Positional argument `guidanceAlg`.
- id: maneuver_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `maneuver_idx`.
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
  type: Float64
  units: n/a
  description: Return value of `_flight_apoapsis_ratio_scale`.
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

# _flight_apoapsis_ratio_scale

## Purpose
Computes a scalar multiplier applied to the pre-tabulated campaign burn `maneuver_Δv` so that the simulated maneuver is rescaled by the ratio of the flight-observed apoapsis radius to the simulation's current osculating apoapsis radius.

## Design & Implementation
Takes the `AerobrakingCampaignPropulsiveManeuverGuidanceModel`, the table index `maneuver_idx`, the state `u::ComponentVector`, `p::ODEParams` and spacecraft index `i`. It returns 1.0 immediately if `guidanceAlg.maneuver_flight_apoapsis_radius_m` is empty or `maneuver_idx` is out of range. Otherwise it extracts `u.sc[i].pos` and `u.sc[i].vel` as `SVector{3,Float64}`, computes osculating elements with `_osculating_elements_and_periapsis_direction` (returns 1.0 if that yields `nothing`, e.g. hyperbolic state), forms `r_a_sim_m = a*(1+e)`, and delegates to `_flight_ratio_scale`, which returns `clamp(r_a_flight/r_a_sim, 0.1, 10.0)` when both radii are finite and positive and 1.0 otherwise. The comment notes the control layer locks the plan at the first post-atmosphere tick, where the osculating apoapsis already reflects the drag pass.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `guidanceAlg` | AerobrakingCampaignPropulsiveManeuverGuidanceModel | n/a | yes | Positional argument `guidanceAlg`. |
| in | `maneuver_idx` | Int | n/a | yes | Positional argument `maneuver_idx`. |
| in | `u` | ComponentVector | n/a | yes | Positional argument `u`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_flight_apoapsis_ratio_scale`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.thruster_guidance_functions_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:125-125`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:108-108`
- `callees` → [[gncz.thruster_guidance_functions__osculating_elements_and_periapsis_direction|_osculating_elements_and_periapsis_direction]] · `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:105-105`
<!-- vulcan:connections:end -->

## Limitations
The ratio is clamped to the hard-coded band [0.1, 10.0] with no diagnostic when clamping engages. The scale is linear in apoapsis radius ratio, which is only a first-order proxy for the actual delta-v sensitivity. Any failure mode (unbound orbit, non-finite state, out-of-range index) silently returns 1.0 rather than signalling that the flight correction was skipped.

## Provenance
Mapped from `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl` line 92.
