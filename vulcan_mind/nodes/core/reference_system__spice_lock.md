---
id: core.reference_system__spice_lock
label: _spice_lock
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: _spice_lock
  lines:
  - 11
  - 11
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
  description: Return value of `_spice_lock`. Returns `getproperty(mod, :RuntimeServices).SPICE_LOCK`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# _spice_lock

## Purpose
Locates the process-wide SPICE lock from inside a file that is included into an arbitrary module, so frame transforms can serialise CSPICE calls without a hard dependency on where they were included.

## Design & Implementation
Starts at `@__MODULE__` and walks `parentmodule` upward until it finds a module defining `RuntimeServices`, returning that module's `SPICE_LOCK`. If it reaches the top without finding one it raises an `ErrorException`. The ancestry walk is what lets this file be included by both the main package and the sandbox test harness.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_spice_lock`. Returns `getproperty(mod, :RuntimeServices).SPICE_LOCK`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system__body_fixed_state_xform|_body_fixed_state_xform]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:48-48`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The walk runs on every call and is not `@inline`, so callers on a hot path should cache the result; the error path only fires if the file is included somewhere without the runtime services module, which is a build-structure mistake rather than a runtime condition.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 11.
