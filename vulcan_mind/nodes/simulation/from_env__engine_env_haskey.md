---
id: simulation.from_env__engine_env_haskey
label: _engine_env_haskey
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _engine_env_haskey
  lines:
  - 207
  - 207
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
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
  description: Return value of `_engine_env_haskey`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _engine_env_haskey

## Purpose
Presence test paired with `_engine_env_get`: reports whether a knob is defined in the active override scope, or in `ENV` when no scope is installed.

## Design & Implementation
`_engine_env_haskey(name::String)::Bool`, `@inline`. If `_engine_active_overrides_ref[]` is not `nothing` it returns `haskey(active_overrides, name)` and never touches `ENV`; otherwise it returns `haskey(ENV, name)`. Used by engine code that needs to distinguish an explicitly set knob from one merely defaulting.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_engine_env_haskey`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/adapters/from_env.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Within an override scope, environment variables that were not folded into the override dict report as absent even though `_with_engine_env_overrides` leaves them in `ENV`. The global `Ref` makes the answer depend on which task most recently installed or restored a scope.

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 207.
