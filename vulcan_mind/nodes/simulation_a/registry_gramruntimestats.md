---
id: simulation_a.registry_gramruntimestats
label: GramRuntimeStats
kind: struct
source:
  file: src/simulation/callbacks/registry.jl
  symbol: GramRuntimeStats
  lines:
  - 40
  - 61
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: stats_events
  type: Function
  units: n/a
  required: true
  description: Mutating closures applied under `_gram_runtime_stats_lock` by `_gram_runtime_stats_update!`
    from the density callback, the track-cache query path and the refresh path.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: stats_snapshot
  type: NamedTuple
  units: n/a
  description: Immutable counter and error-magnitude snapshot returned by `_gram_runtime_stats_snapshot`
    for profiling reports.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation-a
origin: agent
---
# GramRuntimeStats

## Purpose
`GramRuntimeStats` is the mutable counter record behind the GRAM profiling facility. It accumulates how often the atmosphere was sampled, how the track cache performed, what a refresh cost, and how far the cached state drifted from the true state, giving a quantitative basis for deciding whether the track cache is worth enabling for a given mission profile.

## Model & Assumptions
The fields separate three concerns. Call accounting covers total density calls, calls made while the cache was enabled, hits, misses and direct calls, with misses further attributed to a time-window failure or a state-tolerance failure so the two distinct causes can be told apart. Refresh accounting records call count, total and maximum sample counts, failure count and cumulative elapsed seconds. Error accounting keeps sample count plus maximum and summed absolute altitude, latitude and longitude errors, from which a mean error is recoverable without storing per-sample history.

## Design & Implementation
The struct is declared with `Base.@kwdef` so every field has a zero default and a fresh instance is created by calling `GramRuntimeStats()`. A single instance lives in the module-level `Ref` `_gram_runtime_stats`, guarded by `_gram_runtime_stats_lock`, a `ReentrantLock`. Mutation always goes through `_gram_runtime_stats_update!`, which takes the lock and applies a caller-supplied closure, while `_gram_runtime_stats_reset!` replaces the record wholesale and `_gram_runtime_stats_snapshot` copies the fields into a `NamedTuple` under the same lock. Callers guard every update with `_gram_runtime_stats_enabled()`, so profiling costs a single boolean test when switched off.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `stats_events` | Function | n/a | yes | Mutating closures applied under `_gram_runtime_stats_lock` by `_gram_runtime_stats_update!` from the density callback, the track-cache query path and the refresh path. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `stats_snapshot` | NamedTuple | n/a | — | Immutable counter and error-magnitude snapshot returned by `_gram_runtime_stats_snapshot` for profiling reports. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/registry.jl`
- [[simulation.registry__gram_runtime_stats_reset_bang|_gram_runtime_stats_reset!]] · `callees` → `callers` · call · `src/simulation/callbacks/registry.jl:71-71`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
State is process-global rather than per-run, so concurrent solves in the same process share one record and must reset between runs to get meaningful figures. Counter updates serialise on a single lock, which is acceptable at the profiling sample rate but makes the facility unsuitable for always-on accounting in heavily threaded runs. Integer counters are 64-bit and will not realistically overflow, but the summed error fields lose precision over very long runs.

## Provenance
Mapped from `src/simulation/callbacks/registry.jl:40-61`.
