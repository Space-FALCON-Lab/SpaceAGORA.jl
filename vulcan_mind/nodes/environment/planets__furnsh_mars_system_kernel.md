---
id: environment.planets__furnsh_mars_system_kernel
label: _furnsh_mars_system_kernel
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _furnsh_mars_system_kernel
  lines:
  - 311
  - 311
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
  description: Return value of `_furnsh_mars_system_kernel`. Returns `_furnsh_first_existing(`.
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

# _furnsh_mars_system_kernel

## Purpose
Loads the Mars satellite-system SPK (`mar097`) so that Mars barycentre-to-body offsets and Phobos/Deimos positions resolve for Mars-centred simulations.

## Design & Implementation
`@inline` wrapper that calls `_furnsh_first_existing(spice_path, ("spk/satellites/mar097_GRAM.bsp", "spk/satellites/mar097.bsp"))`, preferring the GRAM-trimmed variant shipped with the GRAM Suite bundle. Throws `ArgumentError` via the callee when neither file exists. Only the `Mars` constructor calls it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spice_path` | String | n/a | yes | Positional argument `spice_path`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_furnsh_mars_system_kernel`. Returns `_furnsh_first_existing(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:408-408`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets__furnsh_first_existing|_furnsh_first_existing]] · `callers` · call · `src/environment/ephemerides/planets.jl:312-312`
<!-- vulcan:connections:end -->

## Limitations
The kernel is mandatory even for missions that never query Phobos or Deimos, so a bundle missing `mar097` prevents constructing a `Mars` at all. Time coverage of the trimmed `_GRAM` file is narrower than the full `mar097.bsp`, and the choice is made purely by file presence.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 311.
