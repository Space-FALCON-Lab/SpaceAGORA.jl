---
id: parallel.outer_route_state_load_outer_route_state_bang
label: load_outer_route_state!
kind: function
source:
  file: src/parallel/routing/outer_route_state.jl
  symbol: load_outer_route_state!
  lines:
  - 198
  - 198
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
- id: replace
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `replace` (default `true`).
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
  description: Return value of `load_outer_route_state!`; mutates `state` in place.
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

# load_outer_route_state!

## Purpose
Restores adaptive outer-route history from a TOML file written by `save_outer_route_state`, either replacing or merging into the in-memory `OuterRouteState`. Returns how many rows and distinct signatures were loaded.

## Design & Implementation
Returns `(path, signatures=0, rows=0)` immediately when the file does not exist or when `history` is not a vector. Under `state.lock`, `replace=true` (default) first calls `empty!(state.history)`. Each row must be an `AbstractDict` with a non-empty stripped `signature`, a `route` in `(:none, :threads, :process)`, and a `stats` table that `_route_payload_stats` accepts; otherwise it is skipped. Stats are merged into `get!`-created buckets by adding `samples`, `successes`, `failures`, `elapsed_sum_s` and `elapsed_sq_sum_s`, so loading the same file twice with `replace=false` doubles every count. Loaded signatures are tracked in a `Set{String}` for the return value.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `state` | OuterRouteState | n/a | yes | Positional argument `state`. |
| in | `path` | AbstractString | n/a | yes | Positional argument `path`. |
| in | `replace` | Bool | n/a | no | Keyword argument `replace` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NamedTuple{(:path, | n/a | — | Return value of `load_outer_route_state!`; mutates `state` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_state.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/parallel/routing/outer_route_state.jl:234-234`
- `callees` → [[parallel.outer_route_state__route_payload_stats|_route_payload_stats]] · `callers` · call · `src/parallel/routing/outer_route_state.jl:219-219`
- `callees` → [[parallel.outer_route_state_outerroutestats|OuterRouteStats]] · `callers` · call · `src/parallel/routing/outer_route_state.jl:226-226`
<!-- vulcan:connections:end -->

## Limitations
`TOML.parsefile` errors (malformed file) propagate uncaught, unlike per-row problems which are silently skipped. `schema_version` is read from nowhere; a future incompatible schema would be merged blindly. Because merging is additive and the file is not deduplicated, repeated non-replacing loads inflate confidence in stale timings. `Symbol(String(get(row, "route", "")))` interns arbitrary strings from the file as symbols before validation, which is harmless but leaks symbol-table memory for adversarial inputs.

## Provenance
Mapped from `src/parallel/routing/outer_route_state.jl` line 198.
