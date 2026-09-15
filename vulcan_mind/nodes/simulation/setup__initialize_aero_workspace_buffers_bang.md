---
id: simulation.setup__initialize_aero_workspace_buffers_bang
label: _initialize_aero_workspace_buffers!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_aero_workspace_buffers!
  lines:
  - 1392
  - 1392
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
  description: Return value of `_initialize_aero_workspace_buffers!`; mutates `p`
    in place. Returns `nothing`.
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

# _initialize_aero_workspace_buffers!

## Purpose
Resets the per-satellite aerodynamic scratch workspace slots to `nothing` at run start so they are lazily allocated on first use for the correct satellite count.

## Design & Implementation
Resizes `aero_workspaces` to the spacecraft count if needed and fills with `nothing`. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_aero_workspace_buffers!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:197-197`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Resizing then filling discards any previously allocated workspaces, so a campaign reusing `SharedBuffers` across runs re-allocates them on the first RHS call.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1392.
