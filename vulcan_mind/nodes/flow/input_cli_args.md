---
id: input.cli_args
label: CLI arguments
kind: external
inputs: []
outputs:
- id: argv
  type: Vector{String}
  units: n/a
  description: Subcommand and --key=value options as typed at the shell.
tags:
- master-flow
charts:
- master
origin: agent
---

# CLI arguments

## Purpose
The command line a user types to drive SpaceAGORA without writing Julia: a subcommand such as `run`, `check-assets` or `verify`, followed by `--key=value` options.

## Design & Implementation
Read by `run_cli` in `src/cli/spaceagora_cli.jl`, which dispatches on the first token and parses the rest into a `SimulationConfiguration` or a verification request. Options mirror the fields of the configuration structs so anything expressible in a script is expressible here.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `argv` | Vector{String} | n/a | — | Subcommand and --key=value options as typed at the shell. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `argv` → [[flow.configure|Configure a run]] · `argv` · dataflow · `src/cli/spaceagora_cli.jl`
<!-- vulcan:connections:end -->

## Limitations
Only the documented subcommands are recognised; an unknown token exits non-zero with a usage message rather than falling through to a default run.
