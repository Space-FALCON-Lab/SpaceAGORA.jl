---
id: simulation.from_env__engine_env_haskey_with_env_fallback
label: _engine_env_haskey_with_env_fallback
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _engine_env_haskey_with_env_fallback
  lines:
  - 228
  - 228
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
  description: Return value of `_engine_env_haskey_with_env_fallback`.
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

# _engine_env_haskey_with_env_fallback

## Purpose
Presence test paired with `_engine_env_get_with_env_fallback`: true if the knob is in the active override dict or, failing that, in the process `ENV`.

## Design & Implementation
`_engine_env_haskey_with_env_fallback(name::String)::Bool`, `@inline`. It short-circuits to `true` when an override dict is active and contains `name`; otherwise it returns `haskey(ENV, name)`. This lets knobs outside the canonical override list keep working when set only in the shell environment while a `SimulationEngineConfig` scope is active.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_engine_env_haskey_with_env_fallback`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- [[simulation.solver_policy__solve_with_explicit_solver|_solve_with_explicit_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:327-327`
- [[simulation.solver_policy__solve_with_fixed_step_solver|_solve_with_fixed_step_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:384-384`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It cannot report `false` for a knob that is present in `ENV` but deliberately excluded from the config, so there is no way to mask an environment variable via the override mechanism. Shares the process-global `_engine_active_overrides_ref` and its task-safety caveat.

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 228.
