---
id: environment.simple_ephemerides_spice_position_j2000_m
label: spice_position_j2000_m
kind: function
source:
  file: src/environment/ephemerides/simple_ephemerides.jl
  symbol: spice_position_j2000_m
  lines:
  - 27
  - 27
inputs:
- id: target
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `target`.
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
- id: observer
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `observer`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `spice_position_j2000_m`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# spice_position_j2000_m

## Purpose
Thread-safe public entry point for querying a body's J2000 position in metres from SPICE, used by n-body perturbation, SRP sun-vector, and ephemeris prewarm code paths.

## Design & Implementation
Acquires the global `SPICE_LOCK` via `lock(SPICE_LOCK) do ... end` and inside calls `_spice_position_j2000_m_unlocked(target, et, observer)`, which performs `spkpos(target, et, "J2000", "none", observer)` and scales by `_SPICE_POSITION_KM_TO_M`. The return is `SVector{3,Float64}`. Being `@inline`, the closure and lock overhead are the only costs beyond the SPICE call itself; the lock is a `ReentrantLock`, so nested calls from a locked context do not deadlock.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `target` | AbstractString | n/a | yes | Positional argument `target`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `observer` | AbstractString | n/a | yes | Positional argument `observer`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `spice_position_j2000_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__nbody_body_position_from_spice_direct_j2000_m|_nbody_body_position_from_spice_direct_j2000_m]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1049-1049`
- [[dynamics.perturbations__nbody_body_position_from_spice_memoized_j2000_m|_nbody_body_position_from_spice_memoized_j2000_m]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1065-1065`
- [[dynamics.perturbations__srp_sun_position_from_spice_direct_j2000_m|_srp_sun_position_from_spice_direct_j2000_m]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1421-1421`
- [[dynamics.perturbations__srp_sun_position_from_spice_memoized_j2000_m|_srp_sun_position_from_spice_memoized_j2000_m]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1437-1437`
- [[dynx.coupled_perturbations_srp|srp]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1269-1269`
- [[grp.src_dynamics_coupled|dynamics/coupled/]] · `members_out` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1049-1049`

**Downstream**

- `callees` → [[environment.simple_ephemerides__spice_position_j2000_m_unlocked|_spice_position_j2000_m_unlocked]] · `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:33-33`
<!-- vulcan:connections:end -->

## Limitations
Every call serialises on the global lock, so multi-threaded RHS evaluations contend heavily; the engine's n-body and SRP ephemeris caches exist specifically to avoid calling this in hot loops. SPICE errors propagate as exceptions. Geometric (uncorrected) positions are returned, so light-time and stellar aberration are neglected, an error of order 1e-4 rad in direction for the Sun.

## Provenance
Mapped from `src/environment/ephemerides/simple_ephemerides.jl` line 27.
