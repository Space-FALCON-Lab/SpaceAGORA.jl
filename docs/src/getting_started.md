# Getting Started

## Clone and instantiate

```bash
git clone https://github.com/Space-FALCON-Lab/SpaceAGORA.jl
cd SpaceAGORA.jl
julia --project=.AGORA -e 'using Pkg; Pkg.instantiate()'
```

## GRAM and data prerequisites

SpaceAGORA can run with local fallback atmosphere models, but high-fidelity GRAM-based studies require a local GRAM installation under `data/GRAMSuite.jl/GRAM Suite 2.0` and the expected SPICE assets.

For the current operational setup, use the repository root `data/GRAMSuite.jl/GRAM Suite 2.0` folder structure and keep host-specific GRAM build outputs local to the machine where GRAM is compiled.

## Run a minimal example

```bash
julia --project=.AGORA src/examples/Earth_Thruster_Test.jl
```

## Run the telemetry verification study

```bash
julia --project=.AGORA benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true
```

## Build the documentation locally

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```
