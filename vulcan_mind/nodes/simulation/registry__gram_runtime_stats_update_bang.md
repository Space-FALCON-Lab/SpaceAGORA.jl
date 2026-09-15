---
id: simulation.registry__gram_runtime_stats_update_bang
label: _gram_runtime_stats_update!
kind: function
source:
  file: src/simulation/callbacks/registry.jl
  symbol: _gram_runtime_stats_update!
  lines:
  - 76
  - 76
inputs:
- id: f
  type: Function
  units: n/a
  required: true
  description: Positional argument `f`.
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
  description: Return value of `_gram_runtime_stats_update!`; mutates `f` in place.
    Returns `nothing`.
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

# _gram_runtime_stats_update!

## Purpose
Applies a caller-supplied mutation closure to the global `GramRuntimeStats` record under the profiling lock, so density callbacks can bump counters from any thread without racing.

## Design & Implementation
Signature `_gram_runtime_stats_update!(f::Function)`. It acquires `_gram_runtime_stats_lock` via the `lock(...) do` form, invokes `f(_gram_runtime_stats[])` passing the mutable struct so `f` can do `s.cache_hits += 1` or accumulate `s.refresh_elapsed_s`, and returns `nothing`. Any exception raised by `f` propagates after the lock is released by the `do`-block unwinding. Callers are expected to check `_gram_runtime_stats_enabled()` first.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | Function | n/a | yes | Positional argument `f`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_gram_runtime_stats_update!`; mutates `f` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:113-113`
- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:282-282`
- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:177-177`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:282-282`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/simulation/callbacks/registry.jl:78-78`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/simulation/callbacks/registry.jl:78-78`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/simulation/callbacks/registry.jl:78-78`
<!-- vulcan:connections:end -->

## Limitations
Because the closure runs while a `ReentrantLock` is held, a long-running `f` serialises every other profiled density call in the process. The function does not itself check the profiling gate, so unconditional use adds lock traffic to the RHS. The return value of `f` is discarded.

## Provenance
Mapped from `src/simulation/callbacks/registry.jl` line 76.
