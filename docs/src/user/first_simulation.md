# First Simulation

Use this page when you are ready to run a fuller packaged scenario after the
initial no-GRAM smoke path.

This page is for users who already have the repository environment instantiated
and want a concrete next command.

Shortest successful command:

```text
julia --project=. examples/AGORA_Earth_NoGRAM.jl
```

`Earth_Thruster_Test.jl`, which this page used to recommend, builds its planet
with `Earth("", SPICE_PATH)` and needs the SPICE kernels shipped in the
`data/GRAMSuite.jl` submodule; on a fresh clone it stops with "Required SPICE
kernel not found". Run it after [GRAMSuite Setup](gramsuite_setup.md).

What to read next:

- [Verification Study](verification_study.md)
- [CLI](../cli.md)
- [Simulation Outputs](outputs.md)
- [Examples Catalog](examples_catalog.md)
- [Concepts](concepts.md)

## Two practical ways to run a first scenario

### Example script

Use a repository-owned example when you want the smallest amount of setup:

```text
julia --project=. examples/AGORA_Earth_NoGRAM.jl
```

It is the quickstart's planet and ephemerides (`make_no_gram_planet(:earth)`,
`SimpleEphemeridesModel()`) with a fuller mission configuration, and it writes
the same three result files under `output/` (no plots).

### CLI wrapper

Use the CLI when you want a stable packaged command surface:

```text
julia --project=. src/cli/main.jl run --example=AGORA_Earth_NoGRAM.jl --output-dir=output/cli_run
```

## When to pick a different path

Choose [Verification Study](verification_study.md) instead when your goal is a
known study workflow with enforcement and report outputs.

Choose [Simulation Outputs](outputs.md) when the run completed and you need to
interpret the CSV, Feather, or manifest files.

Choose [Assets & Modes](../assets.md) instead when you need to decide whether a
machine is ready for GRAM/SPICE-backed runs.
