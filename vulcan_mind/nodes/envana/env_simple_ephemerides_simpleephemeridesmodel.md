---
id: envana.env_simple_ephemerides_simpleephemeridesmodel
label: SimpleEphemeridesModel
kind: struct
source:
  file: src/environment/ephemerides/simple_ephemerides.jl
  symbol: SimpleEphemeridesModel
  lines:
  - 4
  - 11
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: EphemeridesModels namespace that includes this file and exports the
    type.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: model
  type: SimpleEphemeridesModel
  units: n/a
  description: Immutable analytic ephemerides model carrying reference epoch and prime
    meridian angle.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- envana
origin: agent
---
# SimpleEphemeridesModel

## Purpose
`SimpleEphemeridesModel` is the kernel-free alternative to `SpiceEphemeridesModel`. It describes planetary rotation analytically so that unit tests, regression scenarios, and portable runs can compute body-fixed frames without furnishing SPICE kernels.

## Theory & Math
The model rotates the local planet-inertial frame about the spin axis through `theta = pm + omega_3 * (et - t_ref)`, where `pm` is `prime_meridian_at_reference_rad` in radians, `omega_3` is the planet spin rate in rad/s taken from the planet struct, `et` is ephemeris time in seconds past J2000 TDB, and `t_ref` is `reference_epoch_seconds` in the same units. The resulting rotation matrix has the standard third-axis form with entries `cos(theta)` and `sin(theta)`, so its determinant is exactly one and its inverse is its transpose.

## Model & Assumptions
`prime_meridian_at_reference_rad` defaults to `NaN`, which the file documents as a sentinel selecting the planet-true prime meridian rather than an angle. Under that sentinel, Earth uses Greenwich mean sidereal time computed in `planet_frame_lpi`, while other planets take a zero angle at the reference epoch. Supplying any finite value instead selects the legacy linear model with that exact offset. `reference_epoch_seconds` defaults to 0.0, meaning J2000 unless overridden.

## Design & Implementation
The type is declared with `@kwdef` so both fields are keyword-constructible with defaults, and it subtypes `AbstractEphemeridesModel` so the exported query functions dispatch on it. Being an immutable struct of two `Float64` fields, it is isbits and can be stored inside integrator parameter structs without heap allocation. The sentinel comparison is a `NaN` check, not an equality test, since `NaN == NaN` is false in IEEE arithmetic.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | EphemeridesModels namespace that includes this file and exports the type. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `model` | SimpleEphemeridesModel | n/a | — | Immutable analytic ephemerides model carrying reference epoch and prime meridian angle. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.no_gram_presets_make_no_gram_density_model|make_no_gram_density_model]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:68-68`
- [[grp.src_core_state|core/state/]] · `members_out` → `callers` · call · `src/core/state/no_gram_presets.jl:68-68`
- [[parcore.no_gram_presets_make_no_gram_environment|make_no_gram_environment]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:92-92`
- [[spaceagora.precompile_workload__spaceagora_precompile_args|_spaceagora_precompile_args]] · `callees` → `callers` · call · `src/precompile_workload.jl:29-29`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:243-243`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The model captures uniform spin only: it ignores precession, nutation, polar motion, and libration, so body-fixed positions drift relative to a SPICE-derived frame over long propagations. It provides no planetary translational ephemeris, so third-body positions still require another source. The `NaN` sentinel means an accidentally propagated `NaN` in a configuration file is interpreted as a valid mode rather than reported as an error.

## Provenance
Read directly from `src/environment/ephemerides/simple_ephemerides.jl:4-11`, including the inline comment block documenting the `NaN` sentinel and the linear rotation law.
