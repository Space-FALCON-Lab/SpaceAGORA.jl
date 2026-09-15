---
id: simulation.assembly__requires_entry_end_callback
label: _requires_entry_end_callback
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _requires_entry_end_callback
  lines:
  - 69
  - 69
inputs:
- id: effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `effectors`.
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
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
  description: Return value of `_requires_entry_end_callback`.
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

# _requires_entry_end_callback

## Purpose
Decides whether to install the callback that detects the end of an atmospheric entry pass, which requires both that entry targets were requested and that density is computable.

## Design & Implementation
Returns `_entry_target_count() > 0 && _requires_density_callback(effectors, args)`. The short-circuit order matters: the cheap-to-fail environment check runs first, so the effector scan is skipped entirely in the common case where no entry targets are configured. It is consulted twice during assembly — once by `_requires_staged_density_callback`, because entry-end detection reads staged density to compare altitude against the entry interface, and once by `get_callbacks` itself to install `get_entry_end_callback(num_sats, args)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effectors` | Tuple | n/a | yes | Positional argument `effectors`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_requires_entry_end_callback`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.assembly__requires_staged_density_callback|_requires_staged_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:42-42`
- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:172-172`

**Downstream**

- `callees` → [[simulation.assembly__entry_target_count|_entry_target_count]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:70-70`
- `callees` → [[simulation.assembly__requires_density_callback|_requires_density_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:70-70`
<!-- vulcan:connections:end -->

## Limitations
`_entry_target_count()` re-reads and re-parses the environment variable on each of those consultations, so a malformed value throws from whichever call site happens to run first, producing an error message that does not mention entry-end callbacks. Being gated on `_requires_density_callback`, an entry target configured against a `NoAtmosphereModel` run is silently ignored rather than reported as a contradictory configuration. In `get_callbacks` the installation is additionally suppressed whenever `backbone_mode` is active.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 69.
