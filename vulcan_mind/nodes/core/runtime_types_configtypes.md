---
id: core.runtime_types_configtypes
label: ConfigTypes
kind: module
source:
  file: src/core/types/runtime_types.jl
  symbol: ConfigTypes
  lines:
  - 2
  - 2
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
  type: Any
  units: n/a
  description: Value produced by this symbol.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# ConfigTypes

## Purpose
Module that declares the runtime data types shared between simulation setup, the ODE right-hand side, callbacks, and result persistence, including `ODEParams`, `SharedBuffers`, ephemeris and atmosphere caches, and the run-scoped environment configuration snapshots.

## Design & Implementation
`module ConfigTypes` imports `SpacecraftModel`, `SimulationConfiguration`, the GRAM atmosphere model types, legacy enum codes with `_compat_enum_parse`, propulsive command types, `StateSample`, `StaticArrays`, `AstroTime`, and `OrdinaryDiffEq`. It exports the legacy Python-port structs (`Mission`, `InitialParameters`, `Model`, `Cnf`, `Solution`), the per-step `IntermediateSolution`, the ephemeris caches, `SpiceRuntimeCounters`, `SpiceRhsMemo`, the GRAM caches, the scratch workspaces, the three env-config structs, and the `RhsEffectorDecision`/`RhsExecutionPlan` NamedTuple aliases. `SharedBuffers` is `@kwdef` with `n_sats` as a runtime field; `ODEParams{A}` deliberately is not `@kwdef` and exposes a keyword outer constructor that validates `args isa SimulationConfiguration`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The module carries a large amount of legacy structure (commented-out `Planet`, `Cnf`, `Orientation`) that duplicates newer typed models, so consumers must know which path is live. `SaveData` is `Dict{Symbol, Any}` by design, giving no type safety at the persistence boundary. Many `Ref{Any}` fields (`harmonics_lpi_key`, `rhs_harmonics_batch_pool`) defeat inference where they are read.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 2.
