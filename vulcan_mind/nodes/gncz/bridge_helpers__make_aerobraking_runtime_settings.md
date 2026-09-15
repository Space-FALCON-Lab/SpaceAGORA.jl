---
id: gncz.bridge_helpers__make_aerobraking_runtime_settings
label: _make_aerobraking_runtime_settings
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _make_aerobraking_runtime_settings
  lines:
  - 145
  - 163
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: GNC internal bridge namespace holding the typed argument accessors
    that each settings field is read through.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: settings
  type: NamedTuple
  units: mixed
  description: Flat aerobraking runtime settings record covering topography, interface
    altitudes, control mode, mass, thrust angle, and integrator selection.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncz
origin: agent
---

# _make_aerobraking_runtime_settings

## Purpose
`_make_aerobraking_runtime_settings` collects every aerobraking switch and scalar that the runtime needs into one flat named tuple, so downstream setup code reads a single settings object instead of probing the caller's argument structure field by field. It is the consolidation point of the bridge between loosely typed scenario arguments and the typed simulation runtime.

## Model & Assumptions
Each field is produced by a dedicated accessor that first looks for a required or optional field on the argument object and, where relevant, falls back to a mission-level value. Optional accessors carry explicit defaults, so a scenario that omits a switch still yields a fully populated settings record. The accessors assume that the argument object supports property queries and that mission values, when consulted, are already validated.

## Design & Implementation
The whole file is written defensively for repeated inclusion. Every definition is wrapped in a check that the name is not already defined in the enclosing module, and the functions are marked inline so the constructed named tuple is built without call overhead on the setup path. A companion routine builds the wider runtime context, embedding this settings record alongside the mission, phase indices, epoch, initial state, and atmosphere handles, and two thin wrappers extend that context with a control gain or a time switch.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | GNC internal bridge namespace holding the typed argument accessors that each settings field is read through. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `settings` | NamedTuple | mixed | — | Flat aerobraking runtime settings record covering topography, interface altitudes, control mode, mass, thrust angle, and integrator selection. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.bridge_helpers__make_aerobraking_runtime_context|_make_aerobraking_runtime_context]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:193-193`

**Downstream**

- `callees` → [[gnc.bridge_helpers__bridge_aerobraking_body_shape|_bridge_aerobraking_body_shape]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:150-150`
- `callees` → [[gnc.bridge_helpers__bridge_aerobraking_dry_mass|_bridge_aerobraking_dry_mass]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:156-156`
- `callees` → [[gnc.bridge_helpers__bridge_aerobraking_entry_interface_m|_bridge_aerobraking_entry_interface_m]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:148-148`
- `callees` → [[gnc.bridge_helpers__bridge_aerobraking_exit_interface_m|_bridge_aerobraking_exit_interface_m]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:149-149`
- `callees` → [[gnc.bridge_helpers__bridge_aerobraking_integrator_name|_bridge_aerobraking_integrator_name]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:159-159`
- `callees` → [[gnc.bridge_helpers__bridge_aerobraking_max_heat_rate|_bridge_aerobraking_max_heat_rate]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:152-152`
- `callees` → [[gnc.bridge_helpers__bridge_aerobraking_thrust_phi|_bridge_aerobraking_thrust_phi]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:157-157`
- `callees` → [[gnc.bridge_helpers__bridge_aerobraking_topography_enabled|_bridge_aerobraking_topography_enabled]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:147-147`
<!-- vulcan:connections:end -->

## Limitations
The named tuple is positional in its field names only, so adding a field changes the type seen by every consumer and forces recompilation. The redefinition guards mean that a stale definition loaded first wins, which can mask an edit during an interactive session.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl:1-233`.
