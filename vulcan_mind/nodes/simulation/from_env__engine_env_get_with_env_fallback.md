---
id: simulation.from_env__engine_env_get_with_env_fallback
label: _engine_env_get_with_env_fallback
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _engine_env_get_with_env_fallback
  lines:
  - 220
  - 220
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: String
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
  type: String
  units: n/a
  description: Return value of `_engine_env_get_with_env_fallback`.
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

# _engine_env_get_with_env_fallback

## Purpose
Environment lookup variant for knobs outside the canonical override set (for example the solver `SAVE_*` switches): it prefers an explicit override when present but otherwise honours the process `ENV`, preserving historical behaviour inside a `SimulationEngineConfig` scope.

## Design & Implementation
`_engine_env_get_with_env_fallback(name::String, default::String)::String`, `@inline`. It reads `_engine_active_overrides_ref[]`; if a dict is active and `haskey(active_overrides, name)` it returns `String(active_overrides[name])`. In every other case it returns `String(get(ENV, name, default))`. The difference from `_engine_env_get` is that an active scope no longer shadows unrelated environment variables. Note that `default` has no default value here, unlike `_engine_env_get`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | String | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_engine_env_get_with_env_fallback`. |
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
Because `_with_engine_env_overrides` also writes overrides into `ENV`, the two branches normally agree; they differ only for the empty-override fast path or if `ENV` is mutated concurrently. The global `Ref` is shared across tasks. Callers must choose between this and `_engine_env_get` by convention; nothing enforces which knobs belong to which set.

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 220.
