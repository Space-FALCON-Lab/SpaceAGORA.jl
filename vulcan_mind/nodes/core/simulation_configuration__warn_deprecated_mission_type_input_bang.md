---
id: core.simulation_configuration__warn_deprecated_mission_type_input_bang
label: _warn_deprecated_mission_type_input!
kind: function
source:
  file: src/core/state/simulation_configuration.jl
  symbol: _warn_deprecated_mission_type_input!
  lines:
  - 16
  - 16
inputs:
- id: mission_type
  type: Any
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
  type: Nothing
  units: n/a
  description: Return value of `_warn_deprecated_mission_type_input!`; mutates `mission_type`
    in place. Returns `nothing`.
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

# _warn_deprecated_mission_type_input!

## Purpose
`_warn_deprecated_mission_type_input!` emits, at most once per Julia process, a `@warn` telling the user that passing `mission_type` as a `String` or `Symbol` is deprecated in favour of the `MissionType` enum. `_parse_mission_type` calls it whenever it successfully parses a legacy string.

## Design & Implementation
Declared `@inline` with signature `(mission_type)`. It returns `nothing` immediately if `_warn_deprecated_config_enabled()` is false or the module-level `_deprecated_mission_type_input_warned[]` flag is already `true`. Otherwise it sets the flag to `true` (the mutation the `!` denotes) and logs `@warn` with `repr(mission_type)` embedded in the message, recommending `MissionTime`/`MissionOrbits`. The return value is always `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `mission_type` | Any | n/a | yes | Positional argument `mission_type`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_warn_deprecated_mission_type_input!`; mutates `mission_type` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.simulation_configuration__parse_mission_type|_parse_mission_type]] · `callees` → `callers` · call · `src/core/state/simulation_configuration.jl:36-36`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/state/simulation_configuration.jl`

**Downstream**

- `callees` → [[core.simulation_configuration__warn_deprecated_config_enabled|_warn_deprecated_config_enabled]] · `callers` · call · `src/core/state/simulation_configuration.jl:17-17`
<!-- vulcan:connections:end -->

## Limitations
The once-only flag is a non-atomic `Ref{Bool}`; two threads parsing legacy mission types simultaneously can both pass the check and warn twice. Because the flag is never reset, tests that assert on the warning must run first or reset the `Ref` manually. The message includes only the offending value, not the call site, so locating the deprecated caller requires the logger's backtrace.

## Provenance
Mapped from `src/core/state/simulation_configuration.jl` line 16.
