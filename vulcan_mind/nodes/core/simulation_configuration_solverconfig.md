---
id: core.simulation_configuration_solverconfig
label: SolverConfig
kind: struct
source:
  file: src/core/state/simulation_configuration.jl
  symbol: SolverConfig
  lines:
  - 73
  - 73
inputs:
- id: solver_mode
  type: Symbol
  units: n/a
  required: false
  description: Field `solver_mode` (default `:tsit5`).
- id: maxiters
  type: Union{Nothing, Int}
  units: n/a
  required: false
  description: Field `maxiters` (default `nothing`).
- id: symplectic_dt_s
  type: Union{Nothing, Float64}
  units: n/a
  required: false
  description: Field `symplectic_dt_s` (default `nothing`).
- id: gravity_backbone_dt_s
  type: Union{Nothing, Float64}
  units: n/a
  required: false
  description: Field `gravity_backbone_dt_s` (default `nothing`).
- id: split_imex_solver
  type: Symbol
  units: n/a
  required: false
  description: Field `split_imex_solver` (default `:kencarp4`).
- id: multirate_slow_dt_s
  type: Union{Nothing, Float64}
  units: n/a
  required: false
  description: Field `multirate_slow_dt_s` (default `nothing`).
- id: multirate_fast_substeps
  type: Int
  units: n/a
  required: false
  description: Field `multirate_fast_substeps` (default `8`).
- id: multirate_slow_solver
  type: Symbol
  units: n/a
  required: false
  description: Field `multirate_slow_solver` (default `:tsit5`).
- id: multirate_fast_solver
  type: Symbol
  units: n/a
  required: false
  description: Field `multirate_fast_solver` (default `:auto_stiff`).
- id: auto_stiff_gravity_tsit5
  type: Bool
  units: n/a
  required: false
  description: Field `auto_stiff_gravity_tsit5` (default `true`).
- id: auto_stiff_switch_max
  type: Int
  units: n/a
  required: false
  description: Field `auto_stiff_switch_max` (default `50`).
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: SolverConfig
  units: n/a
  description: Constructed `SolverConfig` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# SolverConfig

## Purpose
`SolverConfig` is the typed, immutable record of solver selection and stepping policy that can be pinned on a `SimulationConfiguration` to override the `SPACEAGORA_SOLVER_*` environment variables. When left as `nothing`, `run_simulation` derives the same fields from the environment at call time.

## Design & Implementation
A `Base.@kwdef struct` with eleven fields: `solver_mode::Symbol = :tsit5`, `maxiters::Union{Nothing,Int}`, `symplectic_dt_s` and `gravity_backbone_dt_s` (`Union{Nothing,Float64}`, seconds, for fixed-step modes), `split_imex_solver::Symbol = :kencarp4` for the atmosphere-implicit IMEX partition, `multirate_slow_dt_s`, `multirate_fast_substeps::Int = 8`, `multirate_slow_solver = :tsit5`, `multirate_fast_solver = :auto_stiff`, `auto_stiff_gravity_tsit5::Bool = true`, and `auto_stiff_switch_max::Int = 50` bounding the number of stiff/non-stiff switches. Each field corresponds one-to-one to an environment variable of the same name, and the docstring notes that `gravity_backbone_split` is a fixed-step symplectic gravity backbone rather than a whole-system symplectic solve.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `solver_mode` | Symbol | n/a | no | Field `solver_mode` (default `:tsit5`). |
| in | `maxiters` | Union{Nothing, Int} | n/a | no | Field `maxiters` (default `nothing`). |
| in | `symplectic_dt_s` | Union{Nothing, Float64} | n/a | no | Field `symplectic_dt_s` (default `nothing`). |
| in | `gravity_backbone_dt_s` | Union{Nothing, Float64} | n/a | no | Field `gravity_backbone_dt_s` (default `nothing`). |
| in | `split_imex_solver` | Symbol | n/a | no | Field `split_imex_solver` (default `:kencarp4`). |
| in | `multirate_slow_dt_s` | Union{Nothing, Float64} | n/a | no | Field `multirate_slow_dt_s` (default `nothing`). |
| in | `multirate_fast_substeps` | Int | n/a | no | Field `multirate_fast_substeps` (default `8`). |
| in | `multirate_slow_solver` | Symbol | n/a | no | Field `multirate_slow_solver` (default `:tsit5`). |
| in | `multirate_fast_solver` | Symbol | n/a | no | Field `multirate_fast_solver` (default `:auto_stiff`). |
| in | `auto_stiff_gravity_tsit5` | Bool | n/a | no | Field `auto_stiff_gravity_tsit5` (default `true`). |
| in | `auto_stiff_switch_max` | Int | n/a | no | Field `auto_stiff_switch_max` (default `50`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SolverConfig | n/a | — | Constructed `SolverConfig` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:128-128`
- [[simulation.from_env__solver_config_from_env|_solver_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:128-128`
- [[simulation.simulation_engine_config|SimulationEngineConfig]] · `callees` → `callers` · call · `src/simulation/engine/config/simulation_engine_config.jl:10-10`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No field validation is performed: negative `dt_s` values, `multirate_fast_substeps <= 0` or unknown `solver_mode` symbols are accepted here and only fail inside the engine. `nothing` defaults mean the effective values depend on the engine's environment fallback, so a partially specified `SolverConfig` mixes explicit and environment-derived settings. The struct is immutable, so changing one field requires reconstruction.

## Provenance
Mapped from `src/core/state/simulation_configuration.jl` line 73.
