# First Simulation

Use this page when you are ready to run a fuller packaged scenario after the
initial no-GRAM smoke path.

This page is for users who already have the repository environment instantiated
and want a concrete next command.

Shortest successful command:

```bash
julia --project=.AGORA examples/Earth_Thruster_Test.jl
```

What to read next:

- [Verification Study](verification_study.md)
- [CLI](../cli.md)
- [Concepts](concepts.md)

## Two practical ways to run a first scenario

### Example script

Use a repository-owned example when you want the smallest amount of setup:

```bash
julia --project=.AGORA examples/Earth_Thruster_Test.jl
```

### CLI wrapper

Use the CLI when you want a stable packaged command surface:

```bash
./bin/spaceagora run --example=AGORA_Earth_NoGRAM.jl --output-dir=output/cli_run
```

## When to pick a different path

Choose [Verification Study](verification_study.md) instead when your goal is a
known study workflow with enforcement and report outputs.

Choose [Assets & Modes](../assets.md) instead when you need to decide whether a
machine is ready for GRAM/SPICE-backed runs.
