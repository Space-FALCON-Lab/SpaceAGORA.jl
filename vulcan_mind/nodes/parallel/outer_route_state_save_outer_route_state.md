---
id: parallel.outer_route_state_save_outer_route_state
label: save_outer_route_state
kind: function
source:
  file: src/parallel/routing/outer_route_state.jl
  symbol: save_outer_route_state
  lines:
  - 146
  - 146
inputs:
- id: state
  type: OuterRouteState
  units: n/a
  required: true
  description: Positional argument `state`.
- id: path
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `path`.
- id: metadata
  type: AbstractDict
  units: n/a
  required: false
  description: Keyword argument `metadata` (default `Dict{String, Any}()`).
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
  type: NamedTuple{(:path,
  units: n/a
  description: Return value of `save_outer_route_state`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# save_outer_route_state

## Purpose
Persists the adaptive outer-route history to a TOML file so route timings survive across Julia sessions. It writes a schema-versioned document with a UTC timestamp, caller-supplied metadata and one row per (signature, route) pair that has at least one sample.

## Design & Implementation
Under `state.lock`, signatures are visited in sorted order for deterministic output; for each, routes are emitted in the fixed order `(:none, :threads, :process)`, skipping entries that are not `OuterRouteStats` or have `samples == 0`. Each row is `Dict("signature", "route" => String(route), "stats" => _route_stats_payload(stats))`. Metadata values that are `Number` or `Bool` are kept, all others are stringified. The payload has `schema_version = 2`, `updated_utc = string(now(UTC))`, `metadata` and `history`. `mkpath(dirname(path))` is called before `TOML.print`. Returns the named tuple `(path, signatures, rows)` where `signatures` counts only signatures that produced at least one row.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `state` | OuterRouteState | n/a | yes | Positional argument `state`. |
| in | `path` | AbstractString | n/a | yes | Positional argument `path`. |
| in | `metadata` | AbstractDict | n/a | no | Keyword argument `metadata` (default `Dict{String, Any}()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NamedTuple{(:path, | n/a | — | Return value of `save_outer_route_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_state.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/parallel/routing/outer_route_state.jl:163-163`
- `callees` → [[parallel.outer_route_state__route_stats_payload|_route_stats_payload]] · `callers` · call · `src/parallel/routing/outer_route_state.jl:166-166`
<!-- vulcan:connections:end -->

## Limitations
The lock is held only while snapshotting rows, not during file I/O, so a concurrent `load_outer_route_state!` cannot interleave with the snapshot but the file may be written after further updates. Writing is not atomic (no temp-file rename), so a crash mid-write leaves a truncated TOML that `load_outer_route_state!` will fail to parse with an exception. Metadata keys are stringified with `String(k)`, which throws for non-string-convertible keys. Any route symbol other than the three listed is silently dropped.

## Provenance
Mapped from `src/parallel/routing/outer_route_state.jl` line 146.
