---
id: gnc.propulsive_maneuvers__trace_bool_enabled
label: _trace_bool_enabled
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _trace_bool_enabled
  lines:
  - 25
  - 25
inputs:
- id: raw
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `raw`.
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
  type: Bool
  units: n/a
  description: Return value of `_trace_bool_enabled`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _trace_bool_enabled

## Purpose
Normalises a raw environment string into a boolean for the maneuver-trace switches.

## Design & Implementation
Lowercases and strips `raw`, then returns whether the resulting token is a member of the tuple `("1", "true", "yes", "on")`. Unlike the density-callback parser, an unrecognised token is not an error: anything outside that set, including the empty string, yields `false`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | AbstractString | n/a | yes | Positional argument `raw`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_trace_bool_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers__maneuver_trace_enabled|_maneuver_trace_enabled]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:31-31`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Silent rejection means a typo such as `ture` or a value of `off` both disable tracing with no diagnostic, so a user who believes tracing is on gets an empty trace file and no explanation. The accepted vocabulary is narrower than `_parse_bool_env`'s, which additionally recognises the negative spellings.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 25.
