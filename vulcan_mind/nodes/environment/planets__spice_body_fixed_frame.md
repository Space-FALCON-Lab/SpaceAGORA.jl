---
id: environment.planets__spice_body_fixed_frame
label: _spice_body_fixed_frame
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _spice_body_fixed_frame
  lines:
  - 338
  - 338
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
  type: String
  units: n/a
  description: Return value of `_spice_body_fixed_frame`.
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

# _spice_body_fixed_frame

## Purpose
Returns the name of the body-fixed reference frame SPICE associates with a planet, used when transforming between J2000 and planet-fixed coordinates.

## Design & Implementation
`@inline` function returning `String`. For `"Moon"` it hard-codes `"MOON_PA_DE421"` (the principal-axes lunar frame from the loaded DE421 frame kernel); for any other name it takes `SPICE_LOCK` and calls `cnmfrm(planet_name)`, discarding the frame ID and returning the resolved frame name such as `IAU_MARS` or `ITRF93` when the association kernel is loaded.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Positional argument `planet_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_spice_body_fixed_frame`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system__body_fixed_to_j2000_state|_body_fixed_to_j2000_state]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:84-84`
- [[core.reference_system__j2000_to_body_fixed_state|_j2000_to_body_fixed_state]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:67-67`
- [[environment.simple_ephemerides__spice_planet_frame_lpi|_spice_planet_frame_lpi]] · `callees` → `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:88-88`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets__furnsh_first_existing|_furnsh_first_existing]] · `callers` · call · `src/environment/ephemerides/planets.jl:440-440`
- `callees` → [[environment.planets__furnsh_first_existing_if_available|_furnsh_first_existing_if_available]] · `callers` · call · `src/environment/ephemerides/planets.jl:382-382`
- `callees` → [[environment.planets__furnsh_mars_pck|_furnsh_mars_pck]] · `callers` · call · `src/environment/ephemerides/planets.jl:405-405`
- `callees` → [[environment.planets__furnsh_mars_system_kernel|_furnsh_mars_system_kernel]] · `callers` · call · `src/environment/ephemerides/planets.jl:408-408`
- `callees` → [[environment.planets__furnsh_planetary_kernel|_furnsh_planetary_kernel]] · `callers` · call · `src/environment/ephemerides/planets.jl:377-377`
- `callees` → [[environment.planets__furnsh_required|_furnsh_required]] · `callers` · call · `src/environment/ephemerides/planets.jl:375-375`
- `callees` → [[environment.planets__gravity_constants_kernel_if_available|_gravity_constants_kernel_if_available]] · `callers` · call · `src/environment/ephemerides/planets.jl:378-378`
- `callees` → [[environment.planets__spice_backed_planet_kwargs|_spice_backed_planet_kwargs]] · `callers` · call · `src/environment/ephemerides/planets.jl:410-410`
- `callees` → [[environment.planets_earth|Earth]] · `callers` · call · `src/environment/ephemerides/planets.jl:371-371`
- `callees` → [[environment.planets_mars|Mars]] · `callers` · call · `src/environment/ephemerides/planets.jl:401-401`
- `callees` → [[environment.planets_moon|Moon]] · `callers` · call · `src/environment/ephemerides/planets.jl:448-448`
- `callees` → [[environment.planets_titan|Titan]] · `callers` · call · `src/environment/ephemerides/planets.jl:432-432`
- `callees` → [[environment.planets_venus|Venus]] · `callers` · call · `src/environment/ephemerides/planets.jl:417-417`
<!-- vulcan:connections:end -->

## Limitations
The lunar frame name assumes `SPICELunaFrameKernel.tf` defines `MOON_PA_DE421`; a different lunar frame kernel breaks the hard-coded constant. `cnmfrm` returns an empty or invalid frame if no PCK for the body is loaded, and this is not checked.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 338.
