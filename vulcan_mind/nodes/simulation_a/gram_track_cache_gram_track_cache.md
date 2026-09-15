---
id: simulation_a.gram_track_cache_gram_track_cache
label: gram_track_cache
kind: struct
source:
  file: src/simulation/callbacks/gram_track_cache.jl
  symbol: gram_track_cache
  lines:
  - 1
  - 4
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: family_includes
  type: Module
  units: n/a
  required: true
  description: 'The four GRAM track-cache files: configuration parsing, interpolation
    kernels, endpoint targeting and segment refresh.'
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: track_cache_family
  type: Module
  units: n/a
  description: Track-cache definitions injected into the enclosing `SimulationCallbacks`
    namespace and consumed by the density runtime.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation-a
origin: agent
---
# gram_track_cache

## Purpose
`gram_track_cache` is the include manifest for the GRAM track-cache subsystem, the optional acceleration layer that replaces per-step GRAM atmosphere queries with interpolation along a precomputed ground track. Keeping the manifest separate from `density_callbacks.jl` lets the whole feature be located, reviewed and disabled as one unit.

## Model & Assumptions
The ordering reflects a strict dependency chain. `config.jl` must load first because it defines the `GramTrackCache` constructor and the environment parsers that everything else calls. `interpolation.jl` follows with the segment lookup and evaluation kernels. `targeting.jl` supplies the endpoint predictors — Keplerian, periapsis, orbital-period and Allen-Eggers entry — that decide how far ahead a segment should reach. `refresh.jl` comes last because it orchestrates all three.

## Design & Implementation
Each line is an `include(joinpath(@__DIR__, "gram_track_cache", ...))` call, resolved relative to the containing file so the include tree is independent of the process working directory. The feature is off by default: `_gram_track_cache_mode` returns `:off` unless `SPACEAGORA_GRAM_TRACK_CACHE` or the legacy `SPACEAGORA_GRAM_SEGMENT_CACHE` variable says otherwise, because benchmark measurements recorded in `config.jl` show refresh cost can dominate on entry cases.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `family_includes` | Module | n/a | yes | The four GRAM track-cache files: configuration parsing, interpolation kernels, endpoint targeting and segment refresh. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `track_cache_family` | Module | n/a | — | Track-cache definitions injected into the enclosing `SimulationCallbacks` namespace and consumed by the density runtime. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/gram_track_cache.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
As an include manifest the file carries no logic and no tests; a new file placed in the `gram_track_cache/` directory is absent from the build until it is listed here. Load-time errors in any listed file are attributed to this include line rather than to the definition that failed.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache.jl:1-4`.
