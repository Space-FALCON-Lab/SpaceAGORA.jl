---
id: simulation.setup__reset_spice_rhs_memo_bang
label: _reset_spice_rhs_memo!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _reset_spice_rhs_memo!
  lines:
  - 1433
  - 1433
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  description: Return value of `_reset_spice_rhs_memo!`; mutates `p` in place. Returns
    `nothing`.
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

# _reset_spice_rhs_memo!

## Purpose
Invalidates the single-entry SPICE memo at run start so a body position memoised at a previous run's epoch is never served to the new one.

## Design & Implementation
Takes the memo's own `ReentrantLock`, sets its `et` to `NaN`, blanks `primary_body_name` and empties the position dictionary, then returns `nothing`. The `NaN` epoch guarantees the next lookup's equality test fails and repopulates. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_reset_spice_rhs_memo!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:203-203`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The memo is per `SharedBuffers`, so this only protects a reused buffer set; the reset is unnecessary but harmless for a freshly constructed one.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1433.
