---
id: simulation.effector_sampling_sample_planet_frame
label: sample_planet_frame
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: sample_planet_frame
  lines:
  - 42
  - 42
inputs:
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  type: PlanetFrameSample
  units: n/a
  description: Return value of `sample_planet_frame`.
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

# sample_planet_frame

## Purpose
Computes the planet-fixed state and geodetic coordinates of one satellite at time `t`, the foundation for atmosphere and gravity-harmonics sampling.

## Design & Implementation
Extracts inertial position and velocity, fetches the rotation `l_pi` at `t`, rotates into the planet-fixed frame with `_planet_relative_state`, and converts to altitude, latitude and longitude with `rtolatlong`. Returns a `PlanetFrameSample` bundling all six results. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PlanetFrameSample | n/a | — | Return value of `sample_planet_frame`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__spacecraft_outside_atmosphere_for_current_state|_spacecraft_outside_atmosphere_for_current_state]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1957-1957`
- [[simulation.effector_sampling__sample_reusable_planet_frame|_sample_reusable_planet_frame]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:262-262`
- [[simulation.effector_sampling_sample_atmosphere|sample_atmosphere]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:116-116`
- [[simx.engine_effector_sampling_sample_environment|sample_environment]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:245-245`

**Downstream**

- `callees` → [[core.effector_sampling_planetframesample|PlanetFrameSample]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:48-48`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:47-47`
- `callees` → [[simulation.effector_sampling__extract_sample_pos_vel|_extract_sample_pos_vel]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:43-43`
- `callees` → [[simulation.effector_sampling__planet_lpi_at_engine|_planet_lpi_at_engine]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:45-45`
- `callees` → [[simulation.planet_frame__planet_relative_state|_planet_relative_state]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:46-46`
<!-- vulcan:connections:end -->

## Limitations
Every call recomputes `l_pi`, which for a SPICE-backed planet is a kernel lookup; callers evaluating many satellites at the same `t` should compute `l_pi` once and use the `_with_lpi` variant.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 42.
