---
id: cli.spaceagora_cli__resolve_example_path
label: _resolve_example_path
kind: function
source:
  file: src/cli/spaceagora_cli.jl
  symbol: _resolve_example_path
  lines:
  - 29
  - 29
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
  description: Return value of `_resolve_example_path`.
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

# _resolve_example_path

## Purpose
Turns the `--example=` argument into an absolute script path, accepting either an existing file path or a bare example name that lives under the repository `examples/` directory.

## Design & Implementation
The `String` method first computes `direct = abspath(name)` and returns it if `isfile(direct)`. Otherwise it builds `candidate = joinpath(EXAMPLES_DIR, _normalize_example_name(name))` where `EXAMPLES_DIR` is `<repo>/examples`, and throws `ArgumentError("Example not found: <name>. Checked <candidate>.")` if that file does not exist. An `AbstractString` method converts to `String` and forwards. The resolved path is passed to `_run_subprocess` as the script to execute under the `.AGORA` project.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_resolve_example_path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.spaceagora_cli__run_example|_run_example]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:92-92`
- [[module.cli|SpaceAGORACLI]] · `api` → `module_api` · call · `src/cli/spaceagora_cli.jl`

**Downstream**

- `callees` → [[cli.spaceagora_cli__normalize_example_name|_normalize_example_name]] · `callers` · call · `src/cli/spaceagora_cli.jl:34-34`
<!-- vulcan:connections:end -->

## Limitations
Because the direct-path check runs first and is relative to the current working directory, a file in the cwd that happens to share a name with an example shadows the repository example. Only a single fallback directory is searched; subdirectories of `examples/` are not walked. The error message reports only the last candidate checked.

## Provenance
Mapped from `src/cli/spaceagora_cli.jl` line 29.
