---
id: simulation.assembly__requires_staged_density_callback
label: _requires_staged_density_callback
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _requires_staged_density_callback
  lines:
  - 33
  - 33
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
  description: Return value of `_requires_staged_density_callback`.
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

# _requires_staged_density_callback

## Purpose
Decides whether the density callback must actually be installed to pre-stage atmospheric samples into `shared_buffers` before each integrator step, as opposed to letting the right-hand side compute density inline.

## Design & Implementation
The source comment states the key insight: the right-hand side computes density inline through `sample_buffered_atmosphere` into `sample_atmosphere`, so it does not need the staged callback; only non-right-hand-side consumers do. The function first requires `_requires_density_for_rhs`, returning `false` otherwise. It then returns `true` for any of four reasons, in order: the explicit debug override `ParallelPolicy.parse_bool_env("SPACEAGORA_FORCE_DENSITY_CALLBACK", false)`; a thermal callback, which fires each step and reads `shared_buffers.densities`; a drag-state callback, which uses staged density to select the next tolerance set; or entry-end detection, which compares altitude against the entry interface using staged density. Falling through all four returns `false`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effectors` | Tuple | n/a | yes | Positional argument `effectors`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_requires_staged_density_callback`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:160-160`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:36-36`
- `callees` → [[simulation.assembly__requires_density_for_rhs|_requires_density_for_rhs]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:34-34`
- `callees` → [[simulation.assembly__requires_drag_state_callback|_requires_drag_state_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:40-40`
- `callees` → [[simulation.assembly__requires_entry_end_callback|_requires_entry_end_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:42-42`
- `callees` → [[simulation.assembly__requires_thermal_callback|_requires_thermal_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:38-38`
<!-- vulcan:connections:end -->

## Limitations
The three consumer predicates each independently re-evaluate `_requires_density_callback`, so the base condition is recomputed up to four times per call; the `@inline` annotation only helps when the compiler can fold it. The environment-variable override is read on every call rather than once at startup, so changing `SPACEAGORA_FORCE_DENSITY_CALLBACK` mid-process changes behaviour inconsistently. The comment names `get_callbacks` as the sole caller, an invariant nothing enforces, and a new non-right-hand-side consumer added elsewhere will silently read stale or unfilled buffers unless its predicate is added to this list.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 33.
