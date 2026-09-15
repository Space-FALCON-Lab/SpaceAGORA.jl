---
id: simulation.planet_frame__planet_relative_state
label: _planet_relative_state
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/planet_frame.jl
  symbol: _planet_relative_state
  lines:
  - 54
  - 54
inputs:
- id: pos_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_ii`.
- id: vel_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel_ii`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: l_pi
  type: SMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Positional argument `l_pi`.
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
  type: Tuple{SVector{3,
  units: n/a
  description: Return value of `_planet_relative_state`.
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

# _planet_relative_state

## Purpose
`_planet_relative_state(pos_ii, vel_ii, planet, l_pi)` converts an inertial (J2000) position and velocity into the rotating planet-fixed frame, producing the planet-relative state that atmospheric density, wind and aerodynamic angle calculations require. The velocity conversion is the part that matters: the atmosphere co-rotates with the body, so drag depends on velocity relative to the rotating air mass, not on inertial velocity.

## Theory & Math
For a frame rotating with angular velocity $\boldsymbol{\omega}$ expressed in body-fixed axes, and $L_{PI}$ the direction-cosine matrix from inertial to planet-fixed,

$$\mathbf{r}_{PP} = L_{PI}\,\mathbf{r}_{II}, \qquad \mathbf{v}_{PP} = L_{PI}\,\mathbf{v}_{II} - \boldsymbol{\omega} \times \mathbf{r}_{PP}.$$

The second term is the transport velocity of the rotating frame. Because $\boldsymbol{\omega}$ is given in body-fixed components, $\mathbf{r}$ must be in the same basis before the cross product; the equivalent inertial-frame form would be $L_{PI}(\mathbf{v}_{II} - \boldsymbol{\omega}_{II} \times \mathbf{r}_{II})$ with $\boldsymbol{\omega}_{II} = L_{PI}^{\mathsf{T}} \boldsymbol{\omega}$. Units are metres and metres per second throughout, with $\boldsymbol{\omega}$ in rad/s.

## Design & Implementation
Position is a pure rotation, `pos_pp = l_pi * pos_ii`. Velocity applies the transport theorem, `vel_pp = l_pi * vel_ii - cross(planet.ω, pos_pp)`. A prominent source comment records why the ordering is load-bearing: `planet.ω` is expressed in the planet-fixed frame, where the pole is `+z`, so the cross product must be taken *after* rotating the position into that frame. Forming the cross product in J2000 axes instead would point the co-rotation velocity along the J2000 pole and misdirect roughly 250 m/s at Mars periapsis. Inputs and outputs are `SVector{3, Float64}` and the rotation is an `SMatrix{3, 3, Float64}`, so the whole computation is stack-allocated; the function is `@inline` and returns a `Tuple` of the two vectors.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos_ii` | SVector{3, Float64} | n/a | yes | Positional argument `pos_ii`. |
| in | `vel_ii` | SVector{3, Float64} | n/a | yes | Positional argument `vel_ii`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `l_pi` | SMatrix{3, 3, Float64} | n/a | yes | Positional argument `l_pi`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{SVector{3, | n/a | — | Return value of `_planet_relative_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl`
- [[simulation.effector_sampling_sample_planet_frame|sample_planet_frame]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:46-46`
- [[simulation.effector_sampling_sample_planet_frame_with_lpi|sample_planet_frame_with_lpi]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:55-55`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`planet.ω` is treated as a constant vector, so precession, nutation and length-of-day variation are ignored and the rotation rate is assumed rigid. Only the first-order transport term is applied — no centripetal or Coriolis correction appears, which is correct for a velocity transform but means callers must not differentiate `vel_pp` to get a planet-fixed acceleration. The function assumes without checking that `l_pi` is orthonormal; a stale or non-orthonormal matrix silently distorts both outputs. The frame convention is implicit in the argument, so passing the transpose produces a plausible-looking but wrong state with no diagnostic.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/planet_frame.jl` line 54.
