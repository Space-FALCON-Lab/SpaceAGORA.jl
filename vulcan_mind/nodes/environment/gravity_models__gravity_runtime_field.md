---
id: environment.gravity_models__gravity_runtime_field
label: _gravity_runtime_field
kind: function
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: _gravity_runtime_field
  lines:
  - 25
  - 25
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: name
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: Any
  units: n/a
  required: true
  description: Positional argument `default`.
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
  description: Return value of `_gravity_runtime_field`. Returns `getproperty(args,
    name)` or `get(args, name, default)` or `getindex(args, name)` or `default`.
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

# _gravity_runtime_field

## Purpose
Tolerant accessor that reads a named runtime option (such as `:n_bodies`, `:gravity_harmonics`, `:L`, `:M`) from the loosely typed `args` object passed to `aerobraking_gravity_force_ii`, which may be a struct, a `Dict{Symbol,Any}`, a `NamedTuple` or `nothing`.

## Design & Implementation
Marked `@inline`. Returns `default` immediately when `args === nothing`. Otherwise it tries, in order, `hasproperty(args, name)` then `getproperty`; `applicable(get, args, name, default)` then `get`; `applicable(getindex, args, name)` then `getindex`. Each probe is guarded so the first strategy the container supports wins. The result type is whatever the container stores, so callers convert with `Int(...)` where needed.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `name` | Symbol | n/a | yes | Positional argument `name`. |
| in | `default` | Any | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gravity_runtime_field`. Returns `getproperty(args, name)` or `get(args, name, default)` or `getindex(args, name)` or `default`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.gravity_models_aerobraking_gravity_force_ii|aerobraking_gravity_force_ii]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:168-168`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/gravity/gravity_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
For a `Dict` that lacks the key, the `getindex` branch is never reached because `get` succeeds with the default, but for a container supporting `getindex` and not `get`, a missing key throws `KeyError` instead of returning `default`. `applicable` checks are dynamic and defeat inlining benefits on hot paths; the aerobraking predictor calls this several times per derivative evaluation. There is no type check on the returned value.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 25.
