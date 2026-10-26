# SpaceAGORA Documentation

This site follows `main`. For a released version, read the `docs/` folder at
its [release tag](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/tags).
The [migration guide](reference/migrating_extension_hooks.md) lists renamed and removed hooks.

SpaceAGORA is a Julia toolkit for spacecraft orbit and atmospheric-flight
simulation. Use this site to run a first example, configure a study, understand
saved results or extend a model through the supported API.

## Start in under a minute

Use this path if you want the fastest successful run on a machine that does not
have GRAM or SPICE configured yet. SpaceAGORA requires Julia 1.12; other
releases are refused (see [Julia version](user/installation_environment.md#Julia-version)).
Terminal, Linux/macOS; the first line creates the checkout and the rest run
inside it:

```text
GIT_LFS_SKIP_SMUDGE=1 git clone --filter=blob:none https://github.com/Space-FALCON-Lab/SpaceAGORA.jl
cd SpaceAGORA.jl
julia --project=. -e "using Pkg; Pkg.instantiate(julia_version_strict=true)"
julia --startup-file=no examples/setup.jl
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

On Windows set the variable first (PowerShell `$env:GIT_LFS_SKIP_SMUDGE = "1"`,
Command Prompt `set GIT_LFS_SKIP_SMUDGE=1`) and run the same `git clone` line
without the prefix; the [Quickstart](user/quickstart.md) shows both forms and
a checkpoint that confirms the setup before the first run.

SpaceAGORA has two working modes. The no-GRAM path, built on
`NoAtmosphereModel()` and `SimpleEphemeridesModel()`, covers first runs, smoke
runs and lightweight studies. The higher-fidelity path needs GRAM and SPICE
assets on the local machine. You do not need it for a first successful run.

## Choose your route

### Run something now

- [Quickstart](user/quickstart.md): install, check the setup and run the first example.
- [Your First Script](user/first_script.md): write your own scenario from that example.

### Build your own simulation

- [Configuring a Simulation](user/simulation_configuration.md) and [Gravity and Forces](user/gravity_forces.md)
- [Adding a Force or Torque](user/custom_effector.md), [Adding an Actuator](user/custom_actuator.md) and [The Integrated State](user/integrated_state.md)
- [Stopping on a Condition](user/stop_conditions.md) and [Mission Workflows](user/mission_workflows.md)
- [Simulation Outputs](user/outputs.md)

### Guidance and control with a frozen atmosphere

Develop a guidance, control or autonomy algorithm against a reproducible Mars
atmosphere without native GRAM.

- [Odyssey Guidance with a Frozen Atmosphere](tutorials/odyssey_surrogate.md)
- [GRAM Grid Surrogates](user/gram_grid_surrogates.md)

### High-fidelity GRAM and SPICE

- [Assets & Modes](assets.md)
- [GRAMSuite Setup](user/gramsuite_setup.md)
- [Verification Study](user/verification_study.md)

### Optional packages

Motion planning, contact scenes, energy-depletion guidance and the interactive
viewer each need their own package installation; each page says which.

- [Motion Planning with HYPR](user/hypr_planning.md)
- [Orbital Scenes with MuJoCo](user/mujoco_scenes.md)
- [Energy Depletion Guidance](user/edg_guidance.md)
- [Interactive Visualization](user/visualization.md)

## How a run fits together

![Run flow: scenario configuration contains the spacecraft collection; individual spacecraft contain rigid, articulated or compliant assemblies. Environment, spacecraft and GNC supply the repeated dynamics and solver loop; callbacks update commands and save values for outputs.](assets/architecture_run_flow.svg)

[How a Simulation Fits Together](user/architecture.md) explains this flow and where each part lives.

## How these pages are organized

- **Get Started**: installation and the first run.
- **Build a Simulation**: configuring a scenario, gravity and forces, stopping conditions, missions, solver and parallel settings.
- **Vehicles and Models**: articulated and compliant vehicles, mesh aerodynamics, plumes and laser interlinks.
- **Atmosphere and Data**: atmosphere models, GRAM grid surrogates, assets and native GRAM.
- **Guidance and Planning**: energy-depletion guidance and motion planning.
- **Tutorials**: complete worked exercises.
- **Results**: output files, columns and visualization.
- **How It Works**: the concepts and the state the solver integrates.
- **Studies**: verification and benchmark workflows.
- **Help**: troubleshooting, recipes and the examples catalog.
- **Extend SpaceAGORA** and **API Reference**: extension contracts, the public API and the CLI.

To cite SpaceAGORA, see [Publications and Citing](publications.md).
