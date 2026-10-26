# SpaceAGORA.jl

[![Julia CI](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/actions/workflows/julia-ci.yml/badge.svg?branch=main)](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/actions/workflows/julia-ci.yml)
[![Documentation](https://img.shields.io/badge/docs-online-blue.svg)](https://space-falcon-lab.github.io/SpaceAGORA.jl/)
[![Julia 1.12](https://img.shields.io/badge/julia-1.12-9558b2.svg)](https://julialang.org/downloads/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

SpaceAGORA is a Julia simulator for spacecraft that fly in or near a planet's
atmosphere. It models aerobraking, entry, aerocapture and orbital flight around
Earth, Mars, Venus, Titan and the Moon, with guidance and control in the loop.
It is developed by the Space-FALCON Lab at the University of Michigan.

![A simulated Mars aerobraking pass in the SpaceAGORA viewer: the Mars Odyssey spacecraft model coloured by its heating, with the Martian limb behind it and a panel listing altitude, density, heat rate and drag.](docs/src/assets/figures/viewer_heating_closeup.jpg)

*A Mars aerobraking pass in the interactive viewer. Spacecraft model: NASA/JPL/Eyes on the Solar System.*

**[Documentation](https://space-falcon-lab.github.io/SpaceAGORA.jl/)** ·
[Quickstart](https://space-falcon-lab.github.io/SpaceAGORA.jl/user/quickstart/) ·
[Tutorial](https://space-falcon-lab.github.io/SpaceAGORA.jl/tutorials/odyssey_surrogate/) ·
[Publications](https://space-falcon-lab.github.io/SpaceAGORA.jl/publications/)

## Get started

SpaceAGORA needs Julia 1.12. In a terminal on Linux or macOS:

```text
GIT_LFS_SKIP_SMUDGE=1 git clone --filter=blob:none https://github.com/Space-FALCON-Lab/SpaceAGORA.jl
cd SpaceAGORA.jl
julia --project=. -e "using Pkg; Pkg.instantiate(julia_version_strict=true)"
julia --startup-file=no examples/setup.jl
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

The third line installs SpaceAGORA. The fourth prepares the examples
environment, which adds the packages the examples use for plots and analysis,
starting from the package versions SpaceAGORA's manifest records. The example
flies an Earth orbit without atmosphere, GRAM or SPICE data, and writes its
results and plots under `output/`. The clone options skip large data files that
the first run does not need. On Windows, set `GIT_LFS_SKIP_SMUDGE` in its own
command first; the
[Quickstart](https://space-falcon-lab.github.io/SpaceAGORA.jl/user/quickstart/)
shows both Windows shells, a setup check and what the run prints.

For a short check run through the command-line tool, from the same directory
(it runs the example in the examples environment):

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Basic_Quickstart.jl --smoke
```

## What you can simulate

- **Atmospheric flight**: aerobraking, entry and aerocapture with analytic,
  NRLMSISE-00, frozen GRAM grid or native NASA GRAM atmospheres, including heat
  rates and heat loads.
- **Gravity and perturbations**: point mass, J2, spherical harmonics, third
  bodies, solar radiation pressure and gravity-gradient torque.
- **Vehicles**: rigid, articulated and compliant spacecraft carrying thrusters,
  reaction wheels, torque rods, panel hinges, laser terminals or actuator kinds
  of your own.
- **Guidance and control**: your own controllers in the loop, energy-depletion
  aerobraking guidance, rendezvous planning and powered lunar descent with
  plume-surface interaction.
- **Missions and campaigns**: operational segments, phases and events within a
  run, and Monte Carlo and constellation ensembles in parallel.
- **Results**: CSV and Feather tables with a manifest, plots, an interactive 3D
  viewer and a packaged verification study against flight telemetry.

Optional companion packages add HYPR motion planning, MuJoCo contact scenes,
energy-depletion guidance and the 3D viewer. Each installs into its own Julia
environment, as the
[installation guide](https://space-falcon-lab.github.io/SpaceAGORA.jl/user/installation_environment/)
describes.

## Learn more

- [Configuring a Simulation](https://space-falcon-lab.github.io/SpaceAGORA.jl/user/simulation_configuration/):
  build your own scenario.
- [Odyssey guidance with a frozen atmosphere](https://space-falcon-lab.github.io/SpaceAGORA.jl/tutorials/odyssey_surrogate/):
  a complete aerobraking exercise without licensed data.
- [Interactive Visualization](https://space-falcon-lab.github.io/SpaceAGORA.jl/user/visualization/):
  open a run in the 3D viewer.
- [GRAMSuite setup](https://space-falcon-lab.github.io/SpaceAGORA.jl/user/gramsuite_setup/):
  the higher-fidelity path with NASA GRAM and SPICE.
- [Extensibility](https://space-falcon-lab.github.io/SpaceAGORA.jl/extensibility/):
  add your own forces, atmospheres, actuators and controllers.

## Citing

If you use SpaceAGORA, please cite it with the metadata in
[CITATION.cff](CITATION.cff); GitHub's "Cite this repository" button turns it
into BibTeX. Papers built on SpaceAGORA are listed on
[Publications](https://space-falcon-lab.github.io/SpaceAGORA.jl/publications/).

## Contributing and license

The [contributor guide](CONTRIBUTING.md) covers repository layout, ownership
rules and release checks. SpaceAGORA is released under the [MIT License](LICENSE).
