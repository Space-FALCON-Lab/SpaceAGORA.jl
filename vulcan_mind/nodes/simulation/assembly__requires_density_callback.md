---
id: simulation.assembly__requires_density_callback
label: _requires_density_callback
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _requires_density_callback
  lines:
  - 25
  - 25
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
  description: Return value of `_requires_density_callback`.
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

# _requires_density_callback

## Purpose
Backward-compatibility alias preserving the original predicate name for callers that only need to know whether density computation is possible at all, rather than whether it must be pre-staged.

## Design & Implementation
A one-line `@inline` forwarding to `_requires_density_for_rhs(effectors, args)`. The source comment records the reason for keeping it: existing thermal, drag-state, entry-end, and test call sites were written against this name before the staged and right-hand-side requirements were split apart, and rewriting them all was not worth the churn. Within this same file, `_requires_entry_end_callback`, `_requires_drag_state_callback`, and `_requires_thermal_callback` all still call through this alias rather than the new name.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effectors` | Tuple | n/a | yes | Positional argument `effectors`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_requires_density_callback`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- [[simulation.assembly__requires_drag_state_callback|_requires_drag_state_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:74-74`
- [[simulation.assembly__requires_entry_end_callback|_requires_entry_end_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:70-70`
- [[simulation.assembly__requires_thermal_callback|_requires_thermal_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:84-84`

**Downstream**

- `callees` → [[simulation.assembly__requires_density_for_rhs|_requires_density_for_rhs]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:26-26`
<!-- vulcan:connections:end -->

## Limitations
Two names for one predicate invites divergence: a future change to the staged-versus-inline split must remember that this alias is the entry point most callers still use. The name is now actively misleading, since answering `true` here does not mean a density callback is installed — that question belongs to `_requires_staged_density_callback`. Nothing marks it deprecated programmatically, so there is no compiler pressure to migrate.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 25.
