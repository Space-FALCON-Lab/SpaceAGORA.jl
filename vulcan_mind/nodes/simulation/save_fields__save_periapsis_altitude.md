---
id: simulation.save_fields__save_periapsis_altitude
label: _save_periapsis_altitude
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_periapsis_altitude
  lines:
  - 73
  - 73
inputs:
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
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
  description: Return value of `_save_periapsis_altitude`. Returns `periapsis_altitudes`.
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

# _save_periapsis_altitude

## Purpose
Save-time getter for each spacecraft's osculating periapsis altitude in metres, the standard figure of merit for tracking orbit decay during an aerobraking campaign.

## Theory & Math
Periapsis radius $r_p = a(1-e)$, where $a$ is the osculating semi-major axis in metres and $e$ the eccentricity from `rvtoorbitalelement`. The saved altitude is $h_p = a(1-e) - R_{p,e}$ with $R_{p,e}$ the planet's equatorial radius. For a hyperbolic state $a<0$ and $e>1$, so the product remains positive but the two-body interpretation no longer describes a closed orbit.

## Design & Implementation
Marked `@inline`. For each spacecraft it reads the inertial position and velocity out of the state, converts them with `rvtoorbitalelement(pos, vel, planet)`, and forms `oe[1] * (1.0 - oe[2]) - planet.Rp_e`, that is semi-major axis times one minus eccentricity, minus the planet's equatorial radius. The result is a plain `Vector{Float64}` of length `num_sats`. Because the elements are osculating, the value changes within a single orbit as perturbations act, which is exactly what makes it useful for watching a drag pass reduce apoapsis energy.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_save_periapsis_altitude`. Returns `periapsis_altitudes`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:183-183`

**Downstream**

- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:79-79`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:77-77`
<!-- vulcan:connections:end -->

## Limitations
Subtracting the equatorial radius ignores planetary oblateness, so for a polar or high-latitude periapsis the reported altitude is biased by the difference between the equatorial and local radius. Nothing rejects a hyperbolic or near-parabolic state, so escape trajectories still produce a number. Deep in a drag pass the osculating elements oscillate strongly and the saved value is not a per-orbit periapsis.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 73.
