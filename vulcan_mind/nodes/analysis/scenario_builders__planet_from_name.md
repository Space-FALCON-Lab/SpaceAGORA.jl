---
id: analysis.scenario_builders__planet_from_name
label: _planet_from_name
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _planet_from_name
  lines:
  - 6
  - 6
inputs:
- id: planet_name
  type: String
  units: n/a
  required: true
  description: Positional argument `planet_name`.
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
  type: Union{Earth, Mars, Moon, Venus}
  units: n/a
  description: Return value of `_planet_from_name`. Returns `Mars("", SPICE_PATH)`
    or `Venus("", SPICE_PATH)` or `Earth("", SPICE_PATH)` or `Moon("", SPICE_PATH)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# _planet_from_name

## Purpose
Instantiates the planet object a telemetry scenario names, with SPICE kernels furnished from the repository's kernel path.

## Design & Implementation
Lowercases and strips the name, then returns `Mars`, `Venus`, `Earth` or `Moon` constructed with an empty label and `SPICE_PATH`. Any other name raises `ArgumentError`. Constructing through the kernel-furnishing constructor is what makes downstream ephemeris and frame calls valid.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Positional argument `planet_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Earth, Mars, Moon, Venus} | n/a | — | Return value of `_planet_from_name`. Returns `Mars("", SPICE_PATH)` or `Venus("", SPICE_PATH)` or `Earth("", SPICE_PATH)` or `Moon("", SPICE_PATH)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_orbit_args|_make_orbit_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:555-555`
- [[analysis.scenario_builders__make_time_aligned_args|_make_time_aligned_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:596-596`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[environment.planets_earth|Earth]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:13-13`
- `callees` → [[environment.planets_mars|Mars]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:9-9`
- `callees` → [[environment.planets_moon|Moon]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:15-15`
- `callees` → [[environment.planets_venus|Venus]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:11-11`
<!-- vulcan:connections:end -->

## Limitations
Four bodies only; `titan`, which `_nbody_primary_name` accepts as an N-body primary, is not constructible here, so a Titan scenario fails at planet resolution.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 6.
