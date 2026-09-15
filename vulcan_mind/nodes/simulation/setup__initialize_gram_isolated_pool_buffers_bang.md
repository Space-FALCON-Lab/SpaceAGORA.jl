---
id: simulation.setup__initialize_gram_isolated_pool_buffers_bang
label: _initialize_gram_isolated_pool_buffers!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_gram_isolated_pool_buffers!
  lines:
  - 1366
  - 1366
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
  description: Return value of `_initialize_gram_isolated_pool_buffers!`; mutates
    `p` in place. Returns `nothing`.
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

# _initialize_gram_isolated_pool_buffers!

## Purpose
Clears the isolated GRAM pool's model and lock vectors so a new run builds its pool fresh.

## Design & Implementation
Calls `empty!` on `gram_isolated_pool_models` and `gram_isolated_pool_locks`. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_gram_isolated_pool_buffers!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:194-194`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The pool is rebuilt on first batch use rather than here, so the first density callback of a run pays construction cost.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1366.
