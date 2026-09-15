---
id: environment.planets__planetary_kernel_override_relpath
label: _planetary_kernel_override_relpath
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _planetary_kernel_override_relpath
  lines:
  - 289
  - 289
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
  type: String
  units: n/a
  description: Return value of `_planetary_kernel_override_relpath`.
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

# _planetary_kernel_override_relpath

## Purpose
Reads the optional `SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH` environment variable that lets a user force a specific planetary ephemeris file relative to the SPICE bundle root.

## Design & Implementation
`@inline` accessor returning `strip(get(ENV, "SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH", ""))` as a `String`. An empty result means no override, and `_furnsh_planetary_kernel` then falls back to its built-in preference list. The variable is read on every call rather than cached.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_planetary_kernel_override_relpath`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets__furnsh_planetary_kernel|_furnsh_planetary_kernel]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:294-294`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The value is used verbatim as a relative path under `spice_path`; absolute paths are not supported and would be joined incorrectly by `joinpath` only on POSIX semantics (Julia's `joinpath` does honour a second absolute argument, but that behaviour is undocumented here). Because planet instances are cached by `(topo_harmonics_file, spice_path)`, changing this variable mid-process has no effect on already-constructed planets.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 289.
