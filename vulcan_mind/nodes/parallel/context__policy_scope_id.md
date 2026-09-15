---
id: parallel.context__policy_scope_id
label: _policy_scope_id
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _policy_scope_id
  lines:
  - 13
  - 13
inputs:
- id: ctx
  type: PolicyContext
  units: n/a
  required: true
  description: Positional argument `ctx`.
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
  type: UInt
  units: n/a
  description: Return value of `_policy_scope_id`.
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

# _policy_scope_id

## Purpose
Derives a stable `UInt` identifier from a `PolicyContext` instance so that dictionary keys for persistent pools can be built from `isbits` values rather than holding references to the context object itself.

## Design & Implementation
An `@inline` one-liner returning `UInt(objectid(ctx))`. `objectid` yields a hash-like identity for the specific heap object, unique among live objects. Wrapping in `UInt` makes the type explicit for the `Tuple{UInt, Symbol}` pool key. The function is pure and allocation-free.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ctx` | PolicyContext | n/a | yes | Positional argument `ctx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | UInt | n/a | — | Return value of `_policy_scope_id`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.context__active_policy_scope_id|_active_policy_scope_id]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:18-18`
- [[parcore.context_with_policy_context|with_policy_context]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:334-334`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`objectid` values can be reused after the object is garbage-collected, so a stale key from a scope that skipped cleanup could alias a new scope; the design relies on `with_policy_context` always running `_destroy_persistent_foreach_scope!`. The id is not stable across processes or Julia sessions and must not be persisted. Because `PolicyContext` currently carries no fields, two contexts are distinguishable only by identity.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 13.
