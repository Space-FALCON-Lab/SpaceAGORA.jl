---
id: core.simulation_configuration__parse_mission_type
label: _parse_mission_type
kind: function
source:
  file: src/core/state/simulation_configuration.jl
  symbol: _parse_mission_type
  lines:
  - 25
  - 25
inputs:
- id: mission_type
  type: MissionType
  units: n/a
  required: true
  description: Positional argument `mission_type`.
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
  type: MissionType
  units: n/a
  description: Return value of `_parse_mission_type`.
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

# _parse_mission_type

## Purpose
`_parse_mission_type` normalises the `mission_type` keyword accepted by `MissionConfiguration` into the `MissionType` enum, accepting the enum itself, a `Symbol`, or a string for backward compatibility. It is also the basis of the `==` overloads that let `MissionType` compare equal to `"Time"` or `:Orbits`.

## Design & Implementation
Three `@inline` methods, all returning `MissionType`. The `MissionType` method is the identity. The `Symbol` method converts to `String` and recurses. The `AbstractString` method computes `key = lowercase(strip(mission_type))`; `"time"` maps to `MissionTime` and either `"orbits"` or `"orbit"` maps to `MissionOrbits`, each after calling `_warn_deprecated_mission_type_input!`. Any other key throws `ArgumentError` naming the value and listing `"Time"` and `"Orbits"` as valid. The `==` overloads wrap the string method in `try/catch` so an invalid string compares `false` rather than throwing.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `mission_type` | MissionType | n/a | yes | Positional argument `mission_type`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | MissionType | n/a | — | Return value of `_parse_mission_type`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.simulation_configuration_missionconfiguration|MissionConfiguration]] · `callees` → `callers` · call · `src/core/state/simulation_configuration.jl:184-184`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/state/simulation_configuration.jl`

**Downstream**

- `callees` → [[core.simulation_configuration__warn_deprecated_mission_type_input_bang|_warn_deprecated_mission_type_input!]] · `callers` · call · `src/core/state/simulation_configuration.jl:36-36`
<!-- vulcan:connections:end -->

## Limitations
Accepted spellings are hard-coded; abbreviations like `"t"` or localised names are rejected. The `try/catch` in the equality overloads swallows every exception type, not only `ArgumentError`. Non-ASCII whitespace is not stripped by `strip`'s default, so a value padded with such characters throws. The deprecation warning fires only on success, so invalid legacy strings give an error without the migration hint.

## Provenance
Mapped from `src/core/state/simulation_configuration.jl` line 25.
