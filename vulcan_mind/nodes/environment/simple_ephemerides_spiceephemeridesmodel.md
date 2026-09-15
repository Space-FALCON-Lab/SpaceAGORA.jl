---
id: environment.simple_ephemerides_spiceephemeridesmodel
label: SpiceEphemeridesModel
kind: struct
source:
  file: src/environment/ephemerides/simple_ephemerides.jl
  symbol: SpiceEphemeridesModel
  lines:
  - 1
  - 1
inputs:
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
  type: SpiceEphemeridesModel
  units: n/a
  description: Constructed `SpiceEphemeridesModel`.
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

# SpiceEphemeridesModel

## Purpose
Marker type selecting the SPICE-backed ephemeris implementation, in which planet frame orientations and body positions are obtained from furnished NAIF kernels rather than analytic approximations.

## Design & Implementation
Declared as an empty `struct SpiceEphemeridesModel <: AbstractEphemeridesModel end`; all behaviour comes from method dispatch on the type. Associated methods in this file are `ephemerides_requires_spice(::SpiceEphemeridesModel) = true`, `ephemerides_time_seconds(initial_time, ::SpiceEphemeridesModel)` (UTC to ET via `utc2et`), `planet_frame_lpi(planet, et, ::SpiceEphemeridesModel)` (via `pxform`), and `ephemerides_cache_key(::SpiceEphemeridesModel) = (:spice,)`. Because it is a singleton it is `isbits` and can be embedded in configuration structs without allocation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SpiceEphemeridesModel | n/a | — | Constructed `SpiceEphemeridesModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__initial_condition_in_j2000|_initial_condition_in_j2000]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:509-509`
- [[core.simulation_configuration_environmentmodel|EnvironmentModel]] · `callees` → `callers` · call · `src/core/state/simulation_configuration.jl:200-200`
- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:126-126`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- [[simulation.setup__validate_ephemerides_support_bang|_validate_ephemerides_support!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:114-114`
- [[vehicle.model__initial_condition_lpi|_initial_condition_lpi]] · `callees` → `callers` · call · `src/vehicle/spacecraft/model.jl:98-98`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The type carries no reference to which kernels are loaded, so correctness depends entirely on the global SPICE state at call time; the cache key cannot distinguish kernel sets. All SPICE access is serialised by `SPICE_LOCK`. Selecting this model without furnished LSK, SPK, and PCK kernels causes `utc2et`, `spkpos`, or `pxform` to throw at simulation setup or in the first RHS evaluation.

## Provenance
Mapped from `src/environment/ephemerides/simple_ephemerides.jl` line 1.
