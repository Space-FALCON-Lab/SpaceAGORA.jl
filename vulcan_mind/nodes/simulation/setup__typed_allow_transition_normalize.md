---
id: simulation.setup__typed_allow_transition_normalize
label: _typed_allow_transition_normalize
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _typed_allow_transition_normalize
  lines:
  - 31
  - 31
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
  description: Return value of `_typed_allow_transition_normalize`. Returns `_engine_env_get("SPACEAGORA_ALLOW_TYPED_NORMALIZE",
    "0") == "1"`.
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

# _typed_allow_transition_normalize

## Purpose
Opt-in switch permitting the engine to normalise a typed state bundle when transitioning between checkpoint or resume phases, which is otherwise refused to protect layout invariants.

## Design & Implementation
Evaluates `_engine_env_get("SPACEAGORA_ALLOW_TYPED_NORMALIZE", "0") == "1"`. The default `"0"` means normalisation of typed states across transitions is disallowed; only the exact string `"1"` turns it on. Marked `@inline` with no arguments, no mutation and no exceptions.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_typed_allow_transition_normalize`. Returns `_engine_env_get("SPACEAGORA_ALLOW_TYPED_NORMALIZE", "0") == "1"`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_reporting_enforce_typed_normalize_policy__enforce_typed_normalize_policy_bang|_enforce_typed_normalize_policy!]] · `callees` → `callers` · call · `src/simulation/engine/reporting.jl:17-17`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/setup.jl:31-31`
<!-- vulcan:connections:end -->

## Limitations
Strict string equality means `"true"`, `"yes"`, or `" 1"` (with whitespace) all leave the feature off with no diagnostic. There is no logging when the override is active, so a run that silently normalised can be hard to distinguish from one that did not.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 31.
