---
id: cli.spaceagora_cli__normalize_example_name
label: _normalize_example_name
kind: function
source:
  file: src/cli/spaceagora_cli.jl
  symbol: _normalize_example_name
  lines:
  - 22
  - 22
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
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
  type: String
  units: n/a
  description: Return value of `_normalize_example_name`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- cli
charts:
- cli
origin: agent
---

# _normalize_example_name

## Purpose
Canonicalises a user-supplied example name from the `spaceagora run --example=` flag into a `.jl` filename so it can be joined onto `EXAMPLES_DIR` by `_resolve_example_path`.

## Design & Implementation
Two `@inline` methods. The `String` method strips surrounding whitespace with `strip(name)`, returns the token unchanged if `endswith(token, ".jl")`, and otherwise appends `".jl"` via `string(token, ".jl")`. The `AbstractString` method converts to `String` and forwards, so `SubString` results from `split` are accepted. Return type is annotated `String`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_normalize_example_name`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.spaceagora_cli__resolve_example_path|_resolve_example_path]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:34-34`
- [[module.cli|SpaceAGORACLI]] · `api` → `module_api` · call · `src/cli/spaceagora_cli.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the exact lowercase suffix `.jl` is recognised; `Example.JL` gets a second extension appended. An empty or whitespace-only name becomes the literal `".jl"`, which then fails lookup with a confusing message. Nested paths such as `sub/dir/name` are accepted verbatim and simply have `.jl` appended.

## Provenance
Mapped from `src/cli/spaceagora_cli.jl` line 22.
