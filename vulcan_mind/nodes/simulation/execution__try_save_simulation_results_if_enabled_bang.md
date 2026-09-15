---
id: simulation.execution__try_save_simulation_results_if_enabled_bang
label: _try_save_simulation_results_if_enabled!
kind: function
source:
  file: src/simulation/engine/execution.jl
  symbol: _try_save_simulation_results_if_enabled!
  lines:
  - 134
  - 134
inputs:
- id: args
  type: Vararg{Any}
  units: n/a
  required: false
  description: Positional argument `args` (variadic).
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
  description: Return value of `_try_save_simulation_results_if_enabled!`; mutates
    `args` in place. Returns `_save_simulation_results_if_enabled!(args...)` or `nothing`.
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

# _try_save_simulation_results_if_enabled!

## Purpose
Best-effort wrapper around `_save_simulation_results_if_enabled!` used on the failure paths of `run_simulation`. When an integration segment throws, the engine still wants to persist whatever was accumulated up to that point, but a second failure inside the save must not mask the original solver error.

## Design & Implementation
Accepts `args...` and forwards them unchanged to `_save_simulation_results_if_enabled!` inside a `try` block, returning its result (the CSV path or `nothing`). Any exception is caught, reported through `@warn "Failed to save partial simulation results after solve failure."` with `reason=sprint(showerror, err)`, and swallowed by returning `nothing`. `run_simulation` invokes it from the `catch` branches of both the checkpoint-loop segment solve and the single-shot solve, immediately before rethrowing the solver exception; the successful-completion path calls the unwrapped function instead.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vararg{Any} | n/a | no | Positional argument `args` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_try_save_simulation_results_if_enabled!`; mutates `args` in place. Returns `_save_simulation_results_if_enabled!(args...)` or `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:344-344`

**Downstream**

- `callees` → [[simulation.execution__save_simulation_results_if_enabled_bang|_save_simulation_results_if_enabled!]] · `callers` · call · `src/simulation/engine/execution.jl:136-136`
<!-- vulcan:connections:end -->

## Limitations
The splatted signature provides no arity or type checking, so an argument-count mistake surfaces only as a swallowed `MethodError` warning rather than a hard failure. All exception types are caught indiscriminately, including `InterruptException`, which can delay a user interrupt until the save attempt finishes. The warning does not include the partial file path, so a half-written CSV may be left on disk without being mentioned.

## Provenance
Mapped from `src/simulation/engine/execution.jl` line 134.
