---
id: environment.planets__furnsh_planetary_kernel
label: _furnsh_planetary_kernel
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _furnsh_planetary_kernel
  lines:
  - 293
  - 293
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
  description: Return value of `_furnsh_planetary_kernel`. Returns `_furnsh_required(spice_path,
    override_relpath)` or `_furnsh_first_existing(`.
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

# _furnsh_planetary_kernel

## Purpose
Loads the planetary ephemeris SPK (a JPL DE series file) needed for heliocentric and barycentric body positions, honouring an environment override for the exact file.

## Design & Implementation
Reads `_planetary_kernel_override_relpath()`; if non-empty it delegates to `_furnsh_required(spice_path, override)` which throws when the file is absent. Otherwise it calls `_furnsh_first_existing` with the fixed preference tuple `de430.bsp`, `de421.bsp`, `de442s.bsp`, `de442.bsp`, `de440s.bsp`, `de440_GRAM.bsp` under `spk/planets/`. Returns the loaded kernel path. Called by every planet constructor.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spice_path` | String | n/a | yes | Positional argument `spice_path`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_furnsh_planetary_kernel`. Returns `_furnsh_required(spice_path, override_relpath)` or `_furnsh_first_existing(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:377-377`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets__furnsh_first_existing|_furnsh_first_existing]] · `callers` · call · `src/environment/ephemerides/planets.jl:298-298`
- `callees` → [[environment.planets__furnsh_required|_furnsh_required]] · `callers` · call · `src/environment/ephemerides/planets.jl:296-296`
- `callees` → [[environment.planets__planetary_kernel_override_relpath|_planetary_kernel_override_relpath]] · `callers` · call · `src/environment/ephemerides/planets.jl:294-294`
<!-- vulcan:connections:end -->

## Limitations
The preference order places de430 ahead of the newer de440 family, which is a deliberate compatibility choice but means a bundle containing both silently uses the older ephemeris. Only one planetary SPK is ever loaded; missions spanning beyond the chosen file's time coverage fail at `spkezr` time rather than at construction.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 293.
