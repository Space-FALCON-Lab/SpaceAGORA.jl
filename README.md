# SpaceAGORA.jl

SpaceAGORA.jl is a Julia toolkit for spacecraft orbit and atmospheric-flight
simulation, including aerobraking, entry, and aerocapture workflows. The
supported package surface is the root `SpaceAGORA` API together with the
package-owned CLI.

## Install

Use the committed root environment as the canonical committed execution environment
for examples, tests, benchmarks, and most local runs. In normal
repository usage there is no bootstrap step beyond instantiating that committed
environment:

```bash
git clone https://github.com/Space-FALCON-Lab/SpaceAGORA.jl
cd SpaceAGORA.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## First no-GRAM run

The fastest first run does not require GRAM or SPICE:

```bash
julia --project=. examples/AGORA_Earth_NoGRAM.jl
```

SpaceAGORA also supports higher-fidelity GRAM/SPICE-backed runs when local
assets are available under `data/GRAMSuite.jl/GRAM Suite 2.0`.

## Docs and CLI

- Docs site: [space-falcon-lab.github.io/SpaceAGORA.jl](https://space-falcon-lab.github.io/SpaceAGORA.jl)
- CLI wrapper: `./bin/spaceagora`
- Asset check: `./bin/spaceagora assets check`
- Verification study: `./bin/spaceagora telemetry quick --output-dir=output/telemetry_cli --enforce=1`

## For contributors

- User docs live in `docs/src/` and are built with `julia --project=docs docs/make.jl`
- Generated API pages come from `docs/public_api_symbols.jl`
- Architecture and quality references live under the docs Maintainer Guide

Generated reports, docs builds, and other local outputs stay in ignored paths
such as `output/`, `docs/build/`, `docs/site/`, and `docs/src/generated/`.
