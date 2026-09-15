---
id: simulation.setup__initialize_harmonics_workspace_buffers_bang
label: _initialize_harmonics_workspace_buffers!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_harmonics_workspace_buffers!
  lines:
  - 1372
  - 1372
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
  description: Return value of `_initialize_harmonics_workspace_buffers!`; mutates
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

# _initialize_harmonics_workspace_buffers!

## Purpose
Resets the per-satellite spherical-harmonics scratch workspace maps to `nothing` so each satellite lazily allocates workspaces sized to the degree and order actually requested in this run.

## Design & Implementation
Resizes `shared_buffers.harmonics_workspaces` to the spacecraft count if needed and fills with `nothing`. Each slot, once allocated by the harmonics effector, holds a `Dict` keyed by a degree-order hash so several harmonics models of different resolution can coexist per satellite. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_harmonics_workspace_buffers!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:195-195`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Clearing discards every allocated workspace map, including the recursion matrices for high-degree fields that are expensive to size; the first harmonics evaluation of every run re-allocates them.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1372.
