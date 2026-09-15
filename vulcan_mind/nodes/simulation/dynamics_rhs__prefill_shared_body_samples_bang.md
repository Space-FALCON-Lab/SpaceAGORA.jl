---
id: simulation.dynamics_rhs__prefill_shared_body_samples_bang
label: _prefill_shared_body_samples!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _prefill_shared_body_samples!
  lines:
  - 1212
  - 1212
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: sc_state
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_state`.
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
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
  description: Return value of `_prefill_shared_body_samples!`; mutates `p` in place.
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

# _prefill_shared_body_samples!

## Purpose
Samples the Sun and third-body positions once per RHS call using the first active satellite, so all satellites share them.

## Design & Implementation
Finds the first active satellite, samples solar ephemeris if any effector needs it and marks the solar prefill, samples third bodies for the first N-body effector, and stores.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `sc_state` | Any | n/a | yes | Positional argument `sc_state`. |
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_prefill_shared_body_samples!`; mutates `p` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1312-1312`
- [[simulation.dynamics_rhs_spacecraft_dynamics_explicit_remainder_bang|spacecraft_dynamics_explicit_remainder!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2100-2100`
- [[simulation.dynamics_rhs_spacecraft_dynamics_implicit_atmosphere_bang|spacecraft_dynamics_implicit_atmosphere!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2001-2001`
- [[simulation.dynamics_rhs_spacecraft_dynamics_slow_bang|spacecraft_dynamics_slow!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1845-1845`
- [[simx.engine_dynamics_rhs_spacecraft_dynamics_bang|spacecraft_dynamics!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1731-1731`

**Downstream**

- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1220-1220`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1220-1220`
- `callees` → [[dynamics.perturbations__harmonics_lpi_at_bang|_harmonics_lpi_at!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1240-1240`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1220-1220`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1220-1220`
- `callees` → [[simulation.effector_sampling__wrench_method_available|_wrench_method_available]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1220-1220`
- `callees` → [[simulation.effector_sampling_sample_solar_ephemeris|sample_solar_ephemeris]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1222-1222`
- `callees` → [[simulation.effector_sampling_sample_third_body_ephemerides|sample_third_body_ephemerides]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1232-1232`
<!-- vulcan:connections:end -->

## Limitations
Third-body samples are computed but the flat path re-fetches them per effector kernel.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1212.
