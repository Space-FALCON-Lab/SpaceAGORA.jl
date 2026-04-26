# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

SpaceAGORA.jl is a Julia package for high-fidelity spacecraft orbit and atmospheric-flight simulation, including aerobraking, entry, and aerocapture scenarios. It wraps SciML/OrdinaryDiffEq for integration, supports optional GRAM-backed atmospheric modeling and SPICE-backed ephemerides, and includes a telemetry verification pipeline and CLI.

## Development environment

Use the base directory project environment (`--project=.`). The `.AGORA/` and `.SpaceAGORA/` subdirectory environments have been removed.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Common commands

```bash
# Run all tests
julia --project=. -e 'include("test/runtests.jl")'

# Run a single test suite file
julia --project=. -e 'include("test/suites/03_persistence_units_and_rotational_tests.jl")'

# Run a CI gate directly
julia --project=. test/ci_architecture_contract_gate.jl

# Run no-GRAM examples (no local GRAM/SPICE required)
julia --project=. examples/AGORA_Earth_NoGRAM.jl
julia --project=. examples/AGORA_Mars_NoGRAM.jl

# CLI
./bin/spaceagora run --example=AGORA_Earth_NoGRAM.jl --output-dir=output/cli_run
./bin/spaceagora telemetry quick --output-dir=output/telemetry_cli --enforce=1
./bin/spaceagora benchmark runtime-analysis smoke --output-dir=output/perf_cli
./bin/spaceagora assets check

# Rebuild ignored outputs
julia --project=. scripts/regenerate_ignored_outputs.jl runtime-analysis quick
julia --project=. scripts/regenerate_ignored_outputs.jl telemetry quick
julia --project=. scripts/regenerate_ignored_outputs.jl docs
```

## Architecture

The codebase is organized as layered aggregator modules with no behavior ownership in the aggregator files themselves.

### Module hierarchy

```
src/SpaceAGORA.jl                  ← public package facade; re-exports and forwards
  src/core/simulation_model.jl     ← composition root; @reexports all physical models,
                                     types, hooks, IO, callbacks, parallel policy
  src/simulation/engine/
    simulation_engine.jl           ← engine aggregator
    public_api.jl                  ← typed entrypoint + env override boundary
    execution.jl                   ← run_simulation lifecycle (setup → solve → persist)
    setup.jl                       ← buffer/cache initialization
    dynamics_rhs.jl                ← spacecraft_dynamics! ODE RHS, build_initial_conditions
    solver_policy.jl               ← explicit / IMEX / multirate solver selection
    persistence.jl                 ← CSV/bundle output forwarding
    resume_checkpoint.jl           ← checkpoint segmentation and reload
  src/analysis/verification/
    telemetry_verification.jl      ← verification pipeline aggregator
  src/cli/spaceagora_cli.jl        ← CLI dispatch
  src/parallel/
    policy/                        ← inner (per-satellite) parallel policy
    routing/                       ← outer route/profile adaptive policy
```

### Key ownership rules (enforced by CI gates)

- Benchmarks and studies → `benchmarks/` only (never under `src/`)
- Runnable examples → `examples/` only
- Dynamics → `src/dynamics/translational/`, `rotational/`, `coupled/`
- Mission configuration → `src/mission/initial_conditions.jl` (initial state) and `src/mission/operations/` (plans, schedules, aerobraking policy)
- Vehicle structure → `src/vehicle/spacecraft/` (composition) and `src/vehicle/structure/` (mass/inertia/geometry)
- Simulation execution → `src/simulation/engine/` (legacy `events/`, `execution/`, `solver_orchestration/` are retired)

### Extension hooks (stable public API)

Custom models are wired in by dispatching on abstract types exported from `SimulationModel`:

| Hook | Abstract type | Purpose |
|---|---|---|
| `calcForceTorque(model, x, p, i)` | `AbstractForceTorqueModel` | translational + rotational wrench |
| `getDensity(model, h, lat, lon, el_time, wind, p)` | `AbstractDensityModel` | scalar atmosphere query |
| `getDensityBatch!(rhos, Ts, winds, ...)` | `AbstractDensityModel` | optional batch atmosphere |
| `calcControlEffect!(model, u, p, t, i)` | `AbstractControlEffectorModel` | update control state |
| `calcControlForceTorque(model, u, p, i, t)` | `AbstractControlEffectorModel` | control force/torque |
| `calcControlMassFlowRate(model, u, p, i, t)` | `AbstractControlEffectorModel` | propellant consumption |

### Simulation call flow

`run_simulation(SimulationConfiguration)` → `execution.jl`:
1. Validate + initialize buffers (heat rates, density, ephemeris caches, memo mode)
2. `build_initial_conditions` → `ODEParams{N}` + `SharedBuffers{N}`
3. Build callback graph via `SimulationModel.get_callbacks(...)`
4. `_build_typed_solver_problem` → `_solve_with_solver_policy` (explicit / IMEX / multirate)
5. Optional checkpoint segmentation loop
6. `_build_results_dataframe` → write CSV (+ optional bundle)

### Data types

- **Input contract**: `SimulationConfiguration` in `src/core/state/simulation_configuration.jl`
- **Runtime mutable state**: `ODEParams{N}` in `src/core/types/runtime_types.jl` carrying `SharedBuffers{N}`
- **Callback output**: `SavedValues(Float64, Dict{Symbol,Any})`
- **GNC commands**: typed structs in `src/gnc/command_types.jl`

## No-GRAM mode

For onboarding or CI without a local GRAM install, use the built-in no-GRAM builders:

```julia
using SpaceAGORA
env = make_no_gram_environment(planet=:earth, atmosphere=:none, EI_km=120.0, wind=false, topography=false)
```

High-fidelity GRAM/SPICE mode expects assets at `data/GRAMSuite.jl/GRAM Suite 2.0`. `Build/` inside that tree is host-specific generated output.

## Test structure

`test/runtests.jl` orchestrates eight numbered suite files (`test/suites/01_*.jl` through `08_*.jl`). The `test/ci_*.jl` files are standalone CI contract/architecture gates that can be run individually. The `test/coverage_*.jl` files are targeted coverage probes.

`run_simulation` uses `isolate_state=true` by default (deep-copies configuration); use `isolate_state=false` only as an expert performance lever.
