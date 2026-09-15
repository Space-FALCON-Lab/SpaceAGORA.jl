---
id: simx.engine_effector_sampling_sample_environment
label: sample_environment
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: sample_environment
  lines:
  - 234
  - 256
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: engine_api
  type: Module
  units: n/a
  required: true
  description: RuntimeServices namespace and the SimulationEngine/SimulationCampaigns
    scope it aggregates, supplying the SPICE lock and the SimulationModel types used
    here.
- id: requirements
  type: EffectorEnvironmentRequirements
  units: n/a
  required: true
  description: 'Per-effector declaration of which environment products are needed:
    planet frame, atmosphere, solar ephemeris and third-body names.'
- id: state_x
  type: AbstractVector
  units: m,m/s
  required: true
  description: Spacecraft state view for the satellite being sampled, supplying inertial
    position and velocity.
- id: sat_index
  type: Int
  units: index
  required: true
  description: One-based spacecraft index used to select per-satellite cache and buffer
    slots.
- id: time_t
  type: Float64
  units: s
  required: true
  description: Mission-elapsed time converted downstream into SPICE ephemeris time.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: environment_sample
  type: EnvironmentSample
  units: kg/m^3,m,K
  description: Bundle carrying the planet plus the optional planet-frame, atmosphere,
    solar and third-body samples actually requested.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simx
origin: agent
---
# sample_environment

## Purpose
`sample_environment` is the demand-driven environment gather for one spacecraft at one instant. Each effector declares what it needs through an `EffectorEnvironmentRequirements` value, and this function evaluates exactly those products and nothing else, so a gravity-only run never pays for an atmosphere query or a SPICE solar lookup.

## Model & Assumptions
There is one implicit dependency in the requirement set: the atmosphere sample is expressed in planet-fixed coordinates, so requesting `atmosphere` forces the planet frame to be computed even when `planet_frame` itself was not requested. The function captures this with `need_planet_frame = req.planet_frame || req.atmosphere`, computes the frame, and then returns it in the sample only when the effector actually asked for it. Third-body sampling is keyed off `req.third_body_names` being non-empty rather than a boolean flag.

## Design & Implementation
Each optional product is a ternary guarded by its requirement flag and defaults to `nothing`, so `EnvironmentSample` fields are `Union{Nothing,T}` and a consumer that reads a product it did not request fails immediately instead of silently using stale data. The planet comes from `p.args.environment_model.planet`. Atmosphere sampling routes through `_sample_atmosphere_from_planet_frame` with the already-computed frame, avoiding a second frame transform, and carries the `write_buffers` keyword, which defaults to false here so speculative sampling does not pollute the per-satellite buffered-atmosphere slots used by the callback path. The function is `@inline` and its sibling `sample_environment_with_reusable_buffers` shares products across effectors within one RHS call.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `requirements` | EffectorEnvironmentRequirements | n/a | yes | Per-effector declaration of which environment products are needed: planet frame, atmosphere, solar ephemeris and third-body names. |
| in | `state_x` | AbstractVector | m,m/s | yes | Spacecraft state view for the satellite being sampled, supplying inertial position and velocity. |
| in | `sat_index` | Int | index | yes | One-based spacecraft index used to select per-satellite cache and buffer slots. |
| in | `time_t` | Float64 | s | yes | Mission-elapsed time converted downstream into SPICE ephemeris time. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `environment_sample` | EnvironmentSample | kg/m^3,m,K | — | Bundle carrying the planet plus the optional planet-frame, atmosphere, solar and third-body samples actually requested. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__gravity_backbone_core_acceleration|_gravity_backbone_core_acceleration]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1501-1501`
- [[simulation.dynamics_rhs__gravity_backbone_kick_acceleration|_gravity_backbone_kick_acceleration]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1518-1518`

**Downstream**

- `callees` → [[parcore.effector_sampling_environmentsample|EnvironmentSample]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:249-249`
- `callees` → [[simulation.effector_sampling__sample_atmosphere_from_planet_frame|_sample_atmosphere_from_planet_frame]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:246-246`
- `callees` → [[simulation.effector_sampling_sample_planet_frame|sample_planet_frame]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:245-245`
- `callees` → [[simulation.effector_sampling_sample_solar_ephemeris|sample_solar_ephemeris]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:247-247`
- `callees` → [[simulation.effector_sampling_sample_third_body_ephemerides|sample_third_body_ephemerides]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:248-248`
<!-- vulcan:connections:end -->

## Limitations
Solar and third-body sampling hit SPICE and must be serialised behind the shared lock, so an effector set that requests them on every spacecraft turns the RHS into a lock-bound serial section unless the ephemeris caches are enabled. The function has no notion of time tolerance: two calls at times differing by a nanosecond both perform full lookups, which is why the buffered variants exist.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl:234-256`, with the reusable-buffer variant at line 285 and the atmosphere kernel at line 60 of the same file.
