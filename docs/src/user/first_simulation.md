# First Simulation

Use this page to run the quickstart's Earth scenario without generating plots,
or to choose between a script and the command-line launcher.

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

It uses the same 12-hour mission configuration as the quickstart, including
`make_no_gram_planet(:earth)` and `SimpleEphemeridesModel()`. It writes the same
three result files under `output/`, without generating plots.

### CLI wrapper

Use the CLI to keep this run's files in their own directory:

```text
julia --project=. src/cli/main.jl run --example=AGORA_Earth_NoGRAM.jl --output-dir=output/cli_run
```

### Your own script

Every type needed to set up a run is exported from the root module, so
`using SpaceAGORA` is the only import. `make_example_config` and
`make_three_body_spacecraft` build a complete configuration from a planet, a
spacecraft, and an initial condition. `make_example_config` defaults to SPICE
ephemerides, so pass `SimpleEphemeridesModel()` for a run without SPICE kernels:

```julia
using SpaceAGORA

planet = make_no_gram_planet(:earth)
spacecraft = make_three_body_spacecraft(
    bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.01, 2.0, 1.0),
    bus_mass=500.0, panel_mass_each=10.0, panel_offset_y=2.0,
    ic=InitialCondition(ra=planet.Rp_e + 500e3, rp=planet.Rp_e + 500e3,
                        i=45.0, ω=0.0, Ω=0.0, ν=0.0),
)
config = make_example_config(
    planet=planet, spacecraft=spacecraft, mission_time=300.0,
    initial_time=InitialTime(year=2024, month=1, day=1),
    dynamic_effectors=(InverseSquaredGravityModel(),),
    ephemerides_model=SimpleEphemeridesModel(), results=false,
)
run_simulation(config)
```

## When to pick a different path

Choose [Verification Study](verification_study.md) instead when your goal is a
known study workflow with enforcement and report outputs.

Choose [Simulation Outputs](outputs.md) when the run completed and you need to
interpret the CSV, Feather, or manifest files.

Choose [Assets & Modes](../assets.md) instead when you need to decide whether a
machine is ready for GRAM/SPICE-backed runs.
