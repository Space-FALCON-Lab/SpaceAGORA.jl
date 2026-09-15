---
id: simulation_a.density_callbacks_density_callbacks
label: density_callbacks
kind: struct
source:
  file: src/simulation/callbacks/density_callbacks.jl
  symbol: density_callbacks
  lines:
  - 1
  - 6
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
  description: 'The six-file density-callback family: configuration snapshot, model
    selection, callback assembly, planet-frame update, vacuum-predicted GRAM cache
    and runtime evaluation.'
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: density_family
  type: Module
  units: n/a
  description: Density callback definitions injected into the enclosing `SimulationCallbacks`
    namespace.
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
# density_callbacks

## Purpose
`density_callbacks` is the include manifest for the atmosphere-density callback family. It exists so that `callbacks.jl` can pull in six related files with one line while the family itself keeps a clear internal ordering, and so that the density subsystem can be reasoned about, reviewed and grown independently of thermal, event and control callbacks.

## Model & Assumptions
The file assumes it is included from inside `module SimulationCallbacks`, after `registry.jl` has performed the shared imports of `EnvironmentModels`, `ParallelPolicy` and the `ConfigTypes` cache structs. Its ordering places `config.jl` first so that `CallbackEnvConfig` construction helpers exist, then `model_selection.jl`, then `assembly.jl`, `planet_frame.jl`, `vacuum_predicted_gram.jl` and finally `runtime.jl`, which consumes definitions from all five predecessors.

## Design & Implementation
Each entry uses `include(joinpath(@__DIR__, "density_callbacks", ...))` so paths resolve relative to the file rather than the process working directory. The subdirectory split mirrors the responsibilities of the family: environment-variable snapshotting, choosing which atmosphere model instance serves a given spacecraft, assembling the full `CallbackSet`, refreshing the planet-fixed rotation matrix, the predictive vacuum GRAM spline cache, and the per-step density evaluation itself.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `family_includes` | Module | n/a | yes | The six-file density-callback family: configuration snapshot, model selection, callback assembly, planet-frame update, vacuum-predicted GRAM cache and runtime evaluation. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `density_family` | Module | n/a | — | Density callback definitions injected into the enclosing `SimulationCallbacks` namespace. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The manifest carries no behaviour and no tests of its own, so a file added to the `density_callbacks/` directory but not listed here is silently absent from the build. Load errors raised by any listed file are reported against this include line.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks.jl:1-6`.
