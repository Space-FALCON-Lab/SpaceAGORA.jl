# SpaceAGORA Documentation

SpaceAGORA is a Julia toolkit for simulating spacecraft orbits and atmospheric
flight. These pages cover running a first example, configuring your own study,
reading the saved results and extending a model through the supported API.

The site follows `main`. For a released version, read the `docs/` folder at
its [release tag](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/tags); the
[migration guide](reference/migrating_extension_hooks.md) lists renamed and
removed hooks.

## A first run without GRAM

This route works on a machine that has neither GRAM nor SPICE set up.
SpaceAGORA requires Julia 1.12, and other releases are refused (see
[Julia version](user/installation_environment.md#Julia-version)). In a
terminal on Linux or macOS, the first line creates the checkout and the others
run inside it:

```text
GIT_LFS_SKIP_SMUDGE=1 git clone --filter=blob:none https://github.com/Space-FALCON-Lab/SpaceAGORA.jl
cd SpaceAGORA.jl
julia --project=. -e "using Pkg; Pkg.instantiate(julia_version_strict=true)"
julia --startup-file=no examples/setup.jl
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

On Windows, set the variable first (PowerShell `$env:GIT_LFS_SKIP_SMUDGE = "1"`,
Command Prompt `set GIT_LFS_SKIP_SMUDGE=1`) and run the same `git clone` line
without the prefix. The [Quickstart](user/quickstart.md) shows both forms and a
check that confirms the setup before the first run.

SpaceAGORA has two working modes. The no-GRAM mode, built on
`NoAtmosphereModel()` and `SimpleEphemeridesModel()`, is enough for first runs,
smoke runs and lightweight studies. The high-fidelity mode needs GRAM and SPICE
assets on the local machine; [Assets & Modes](assets.md) lists them.

## Where to go next

To run something now, follow the [Quickstart](user/quickstart.md), then turn
its example into your own scenario with [Your First Script](user/first_script.md).

To build your own simulation:

- [Configuring a Simulation](user/simulation_configuration.md) and [Gravity and Forces](user/gravity_forces.md)
- [Adding a Force or Torque](user/custom_effector.md), [Adding an Actuator](user/custom_actuator.md) and [The Integrated State](user/integrated_state.md)
- [Stopping on a Condition](user/stop_conditions.md) and [Mission Workflows](user/mission_workflows.md)
- [Simulation Outputs](user/outputs.md)

To develop guidance, control or autonomy against a reproducible Mars
atmosphere without native GRAM, work through
[Odyssey Guidance with a Frozen Atmosphere](tutorials/odyssey_surrogate.md) and
read [GRAM Grid Surrogates](user/gram_grid_surrogates.md).

For high-fidelity GRAM and SPICE runs, start with [Assets & Modes](assets.md),
then [GRAMSuite Setup](user/gramsuite_setup.md) and the
[Verification Study](user/verification_study.md).

Motion planning, contact scenes, energy-depletion guidance and the interactive
viewer each need a separate package installation, described on their pages:
[Motion Planning with HYPR](user/hypr_planning.md),
[Orbital Scenes with MuJoCo](user/mujoco_scenes.md),
[Energy Depletion Guidance](user/edg_guidance.md) and
[Interactive Visualization](user/visualization.md).

## How a run fits together

![Run flow: scenario configuration contains the spacecraft collection; individual spacecraft contain rigid, articulated or compliant assemblies. Environment, spacecraft and GNC supply the repeated dynamics and solver loop; callbacks update commands and save values for outputs.](assets/architecture_run_flow.svg)

[How a Simulation Fits Together](user/architecture.md) explains this flow and
where each part lives in the source.

## How the site is organized

The sidebar groups the pages as follows.

- **Get Started**: installation and the first run.
- **Build a Simulation**: configuration, gravity and forces, custom effectors and actuators, power, stop conditions, mission workflows, solver and parallel settings.
- **Vehicles and Models**: articulated and compliant vehicles, mesh aerodynamics, constellations, plume interaction and interlinks.
- **Atmosphere and Data**: atmosphere models, GRAM grid surrogates, assets, native GRAM and lunar imagery.
- **Guidance and Planning**: writing GNC algorithms, feedback controllers and estimators, planners, energy-depletion guidance, HYPR motion planning and RPO planners.
- **Tutorials**: complete worked exercises.
- **Results**: output files and columns, the interactive viewer and MuJoCo scenes.
- **How It Works**: concepts, how a simulation fits together and the integrated state.
- **Studies**: the verification study and the benchmark workflows.
- **Help**: troubleshooting, recipes, the examples catalog and upgrade notes.
- **Extend SpaceAGORA** and **API Reference**: extension contracts, the public API and the CLI.

To cite SpaceAGORA, see [Publications and Citing](publications.md).
