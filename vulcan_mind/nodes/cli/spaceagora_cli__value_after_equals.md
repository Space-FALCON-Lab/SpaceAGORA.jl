---
id: cli.spaceagora_cli__value_after_equals
label: _value_after_equals
kind: function
source:
  file: src/cli/spaceagora_cli.jl
  symbol: _value_after_equals
  lines:
  - 18
  - 18
inputs:
- id: arg
  type: String
  units: n/a
  required: true
  description: Positional argument `arg`.
- id: prefix
  type: String
  units: n/a
  required: true
  description: Positional argument `prefix`.
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
  type: Any
  units: n/a
  description: Return value of `_value_after_equals`. Returns `split(arg, "=", limit=2)[2]`.
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

# _value_after_equals

## Purpose
Extracts the value portion of a `--flag=value` command-line token for the CLI parsers, returning everything after the first `=`.

## Design & Implementation
`_value_after_equals(arg::String, prefix::String)` calls `split(arg, "=", limit=2)[2]`, so a value containing further `=` characters is preserved intact (for example a Windows path with query-like content). The `prefix` argument is accepted for call-site symmetry with `_starts_with` but is not used in the computation. The result is a `SubString{String}`; downstream code wraps it with `abspath`, `strip`, `lowercase` or `String` as needed.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `arg` | String | n/a | yes | Positional argument `arg`. |
| in | `prefix` | String | n/a | yes | Positional argument `prefix`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_value_after_equals`. Returns `split(arg, "=", limit=2)[2]`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.spaceagora_cli__run_benchmark|_run_benchmark]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:156-156`
- [[cli.spaceagora_cli__run_example|_run_example]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:80-80`
- [[cli.spaceagora_cli__run_telemetry|_run_telemetry]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:115-115`
- [[module.cli|SpaceAGORACLI]] · `api` → `module_api` · call · `src/cli/spaceagora_cli.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
If the token has no `=`, `split` returns a one-element vector and indexing `[2]` throws a `BoundsError` rather than a descriptive `ArgumentError`; callers guard this only by checking the prefix first. The unused `prefix` parameter is misleading. The returned `SubString` keeps the whole original token alive.

## Provenance
Mapped from `src/cli/spaceagora_cli.jl` line 18.
