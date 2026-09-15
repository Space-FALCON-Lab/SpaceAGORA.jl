---
id: core.simulation_configuration_simconfig
label: SimConfig
kind: module
source:
  file: src/core/state/simulation_configuration.jl
  symbol: SimConfig
  lines:
  - 1
  - 1
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

# SimConfig

## Purpose
`SimConfig` is the module that defines the typed configuration tree for a SpaceAGORA run: the `MissionType` enum, the leaf settings structs (`SolverConfig`, `InitialTime`, `IntegrationTolerances`, `FilePaths`, `SimulationSettings`, `MissionConfiguration`, `EnvironmentModel`) and the root `SimulationConfiguration` that the engine consumes. It also owns the deprecation shims that accept string or symbol mission types.

## Design & Implementation
The module imports abstract planet, density, thermal and ephemerides types from `AbstractTypes`, the `DynamicsModel`/`GuidanceModel`/`ControlModel`/`NavigationModel` containers from `SpacecraftModels`, `SpiceEphemeridesModel` as a default, and `Earth`. Structs are declared with `Base.@kwdef` so every field has a keyword default, except `MissionConfiguration` and `EnvironmentModel`, which use inner constructors for range validation and a separate keyword outer constructor. `MissionType` is a `UInt8` enum (`MissionTime = 0x01`, `MissionOrbits = 0x02`); `Base.:(==)` is overloaded between `MissionType` and `AbstractString`/`Symbol` so legacy `mission_type == "Time"` checks still work. A process-global `Ref(false)` throttles the deprecation warning to one emission. Everything listed in the two `export` lines is re-exported by the package.

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

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/state/simulation_configuration.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Overloading `==` between an enum and strings is type-piracy-adjacent and can surprise generic code (for example `"Time" == MissionTime` is `true` but `hash` differs, so `Dict` lookups by either key disagree). The `_deprecated_mission_type_input_warned` flag is a plain `Ref`, not thread-safe, so concurrent first calls may warn twice. A source comment above `EnvironmentModel` records an unfinished migration from string dispatch to abstract types. `Earth` is imported but unused within this module.

## Provenance
Mapped from `src/core/state/simulation_configuration.jl` line 1.
