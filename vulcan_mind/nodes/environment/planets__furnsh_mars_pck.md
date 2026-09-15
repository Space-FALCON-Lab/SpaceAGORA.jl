---
id: environment.planets__furnsh_mars_pck
label: _furnsh_mars_pck
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _furnsh_mars_pck
  lines:
  - 321
  - 321
inputs:
- id: spice_path
  type: String
  units: n/a
  required: true
  description: Positional argument `spice_path`.
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
  description: Return value of `_furnsh_mars_pck`. Returns `kernel_path` or `_furnsh_required(spice_path,
    "pck/mars_iau2000_m2_quadratic_patch.tpc")`.
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

# _furnsh_mars_pck

## Purpose
Ensures the Mars body-fixed frame and radii constants are available in the SPICE kernel pool, preferring modern generic PCKs and falling back to the legacy pck00008 plus a local IAU2000 quadratic patch.

## Design & Implementation
First loops over `("pck/pck00011.tpc", "pck/pck00010.tpc")`, and on the first `isfile` hit calls `_furnsh_once` and returns the path. If neither exists it calls `_furnsh_required(spice_path, "pck/pck00008.tpc")` followed by `_furnsh_required(spice_path, "pck/mars_iau2000_m2_quadratic_patch.tpc")`, returning the patch path. Either `_furnsh_required` call throws `ArgumentError` when its file is missing.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spice_path` | String | n/a | yes | Positional argument `spice_path`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_furnsh_mars_pck`. Returns `kernel_path` or `_furnsh_required(spice_path, "pck/mars_iau2000_m2_quadratic_patch.tpc")`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:405-405`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets__furnsh_once|_furnsh_once]] · `callers` · call · `src/environment/ephemerides/planets.jl:327-327`
- `callees` → [[environment.planets__furnsh_required|_furnsh_required]] · `callers` · call · `src/environment/ephemerides/planets.jl:334-334`
<!-- vulcan:connections:end -->

## Limitations
The two branches load different rotational models for Mars (IAU 2015 versus IAU 2000 with a patch), so results differ slightly depending on which files the bundle ships, with no warning emitted. The function loads at most one modern PCK; if both pck00011 and pck00010 are present only pck00011 is furnished.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 321.
