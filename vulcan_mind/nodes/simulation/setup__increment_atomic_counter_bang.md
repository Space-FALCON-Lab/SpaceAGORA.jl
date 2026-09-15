---
id: simulation.setup__increment_atomic_counter_bang
label: _increment_atomic_counter!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _increment_atomic_counter!
  lines:
  - 1541
  - 1541
inputs:
- id: counter
  type: Any
  units: n/a
  required: true
  description: Positional argument `counter`.
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
  description: Return value of `_increment_atomic_counter!`; mutates `counter` in
    place. Returns `nothing`.
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

# _increment_atomic_counter!

## Purpose
Bumps an optional SPICE call counter by one without a lock, tolerating a `nothing` counter so call sites need no conditional.

## Design & Implementation
Returns immediately for `nothing`, otherwise `Threads.atomic_add!(counter, 1)`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `counter` | Any | n/a | yes | Positional argument `counter`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_increment_atomic_counter!`; mutates `counter` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_setup_build_nbody_ephemeris_cache|_build_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1567-1567`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Accepting `nothing` means a typo in the counter field name at a call site silently disables counting instead of failing.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1541.
