---
id: simulation.solver_policy__retcode_is_stiff_symptom
label: _retcode_is_stiff_symptom
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _retcode_is_stiff_symptom
  lines:
  - 55
  - 55
inputs:
- id: retcode
  type: Any
  units: n/a
  required: true
  description: Positional argument `retcode`.
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
  description: Return value of `_retcode_is_stiff_symptom`.
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

# _retcode_is_stiff_symptom

## Purpose
Classifies a DiffEq return code as a symptom of stiffness so callers can decide whether an implicit fallback is warranted.

## Design & Implementation
Converts `retcode` to a `Symbol` once (zero-allocation for `Symbol` and `ReturnCode` inputs) and compares with `===` against `:Unstable`, `:DtLessThanMin`, `:MaxIters`, and `:InitialFailure`. Returns `true` on any match.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `retcode` | Any | n/a | yes | Positional argument `retcode`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_retcode_is_stiff_symptom`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The list is hard-coded; `ReturnCode.Failure` or `ConvergenceFailure` are not treated as stiff symptoms. `MaxIters` can also arise from a too-small `maxiters` setting rather than stiffness, producing a false positive. `String` inputs allocate a `Symbol`.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 55.
