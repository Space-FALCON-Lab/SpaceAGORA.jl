---
id: ext.spaceagoragramsuiteext__gram_spice_ephemeris_state
label: _gram_spice_ephemeris_state
kind: function
source:
  file: ext/SpaceAGORAGRAMSuiteExt.jl
  symbol: _gram_spice_ephemeris_state
  lines:
  - 73
  - 73
inputs:
- id: planet_name
  type: String
  units: n/a
  required: true
  description: Positional argument `planet_name`.
- id: initial_time
  type: Any
  units: n/a
  required: true
  description: Positional argument `initial_time`.
- id: el_time
  type: Float64
  units: n/a
  required: true
  description: Positional argument `el_time`.
- id: lat_deg
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lat_deg`.
- id: lon_deg
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lon_deg`.
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
  description: Return value of `_gram_spice_ephemeris_state`. Returns `(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- ext
charts:
- ext
origin: agent
---

# _gram_spice_ephemeris_state

## Purpose
Computes the solar geometry GRAM needs at a query point using SpaceAGORA's working SPICE bindings, bypassing the vendored library's private CSPICE instance whose default kernels fail for non-Earth bodies.

## Theory & Math
With subsolar latitude $\phi_s$, subsolar longitude $\lambda_s$ and query point $(\phi, \lambda)$, the solar zenith angle is

$$
\theta_z = \arccos\left( \sin\phi \sin\phi_s + \cos\phi \cos\phi_s \cos(\lambda - \lambda_s) \right)
$$

and local solar time in hours is $t_{\odot} = \left(12 + H/15\right) \bmod 24$ where $H$ is the hour angle in degrees, wrapped to the interval from $-180$ to $180$.

## Design & Implementation
Rejects any body absent from the `_GRAM_SECONDS_PER_SOL` table by returning `nothing`. Otherwise it converts the epoch through `SPICE.utc2et`, adds `el_time`, and assembles eight quantities matching `EphemerisStateC` field conventions: orbital radius from `spkpos` relative to the Sun scaled by the astronomical unit in kilometres, planetocentric solar longitude from `lspcn`, one-way light time to NAIF body 399 from `ltime` divided by sixty, and the subsolar point from `subslr` decomposed by `reclat`. Local solar time comes from the hour angle at fifteen degrees per hour, and solar zenith angle from the spherical cosine law with the argument clamped to the closed interval from minus one to one before `acos`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Positional argument `planet_name`. |
| in | `initial_time` | Any | n/a | yes | Positional argument `initial_time`. |
| in | `el_time` | Float64 | n/a | yes | Positional argument `el_time`. |
| in | `lat_deg` | Float64 | n/a | yes | Positional argument `lat_deg`. |
| in | `lon_deg` | Float64 | n/a | yes | Positional argument `lon_deg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gram_spice_ephemeris_state`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.ext|SpaceAGORAGRAMSuiteExt]] · `api` → `module_api` · call · `ext/SpaceAGORAGRAMSuiteExt.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:138-138`
- `callees` → [[environment.density_models__gram_core_density_state|_gram_core_density_state]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:234-234`
- `callees` → [[environment.density_models__gram_point_density|_gram_point_density]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:255-255`
- `callees` → [[environment.density_models_density_polyfit|density_polyfit]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:288-288`
- `callees` → [[environment.density_models_gramatmospheremodel|GRAMAtmosphereModel]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:122-122`
- `callees` → [[environment.density_models_gramatmospheremodelsurrogate|GRAMAtmosphereModelSurrogate]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:126-126`
- `callees` → [[environment.density_models_precompute_gram_static_grids_bang|precompute_gram_static_grids!]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:149-149`
- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:271-271`
- `callees` → [[ext.gram_core_density_state|_gram_core_density_state]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:234-234`
- `callees` → [[ext.spaceagoragramsuiteext__gram_call_lock|_gram_call_lock]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:264-264`
- `callees` → [[ext.spaceagoragramsuiteext__gram_utc_string|_gram_utc_string]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:83-83`
- `callees` → [[grp.src_analysis_verification|analysis/verification/]] · `members_in` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:138-138`
- `callees` → [[grp.src_environment_atmosphere|environment/atmosphere/]] · `members_in` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:271-271`
<!-- vulcan:connections:end -->

## Limitations
Seconds per sol are fixed literals copied from the vendored C++ rather than derived from the ephemeris, so they do not track any body not in the table; the whole computation depends on the relevant SPICE kernels being furnished, and `utc2et` throws rather than returning `nothing` when they are not.

## Provenance
Mapped from `ext/SpaceAGORAGRAMSuiteExt.jl` line 73.
