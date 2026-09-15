---
id: cli.spaceagora_cli__starts_with
label: _starts_with
kind: function
source:
  file: src/cli/spaceagora_cli.jl
  symbol: _starts_with
  lines:
  - 14
  - 14
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
  description: Return value of `_starts_with`. Returns `startswith(arg, prefix)`.
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

# _starts_with

## Purpose
Thin `@inline` wrapper around `Base.startswith` used by the CLI argument parsers to test whether a token begins with a `--flag=` prefix before extracting its value.

## Design & Implementation
`_starts_with(arg::String, prefix::String)` returns `startswith(arg, prefix)`. It exists so the parsing loops in `_run_example`, `_run_telemetry` and `_run_benchmark` read uniformly alongside `_value_after_equals`, and it constrains both arguments to `String`, giving a clearer method error if a `SubString` slips in.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `arg` | String | n/a | yes | Positional argument `arg`. |
| in | `prefix` | String | n/a | yes | Positional argument `prefix`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_starts_with`. Returns `startswith(arg, prefix)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.spaceagora_cli__run_benchmark|_run_benchmark]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:155-155`
- [[cli.spaceagora_cli__run_example|_run_example]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:79-79`
- [[cli.spaceagora_cli__run_telemetry|_run_telemetry]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:114-114`
- [[module.cli|SpaceAGORACLI]] · `api` → `module_api` · call · `src/cli/spaceagora_cli.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It adds no behaviour over `startswith` and does not accept `AbstractString`, so callers must convert substrings first. Matching is case-sensitive and exact, so `--Example=` is not recognised.

## Provenance
Mapped from `src/cli/spaceagora_cli.jl` line 14.
