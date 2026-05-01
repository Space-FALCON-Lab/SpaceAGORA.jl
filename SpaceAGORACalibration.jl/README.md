# SpaceAGORACalibration.jl

`SpaceAGORACalibration.jl` is a standalone orchestration package for calibration workflows.
It does not implement spacecraft physics. It drives calibration through a backend interface.

## Architecture

- `Spec`: versioned calibration campaign specification and validation.
- `ParamSpace`: candidate representation and random/LHS sampling/perturbation.
- `Backend`: evaluation interface (`AbstractBackend`) with mock and command-backed implementations.
- `Objective`: score helpers and robust aggregation.
- `GlobalSearch`: surrogate-assisted Bayesian global search (EI/LCB).
- `LocalRefine`: trust-region local refinement with restart support.
- `Robustness`: uncertainty Monte Carlo and robust ranking (`weighted`, `cvar`, `p95`).
- `Store`: immutable run store (`spec.toml`, `evaluations.arrow`, `state.toml`, `best_manifest.toml`).
- `Engine`: restartable lifecycle orchestration.
- `Report`: markdown run summary artifact.

## Spec Contract

The versioned TOML spec includes:

- `schema_version`
- `manifest_paths`
- `scenario_weights`
- parameter list with bounds/types and manifest target mappings
- budgets and robustness settings (`initial_design`, `parallel_evaluations`, `global_acquisition`, `bo_*`, `local_refine_*`, `robustness_*`, `robust_*`)
- objective settings (`objective_huber_delta`, `objective_lambda_fail`, `objective_lambda_time`, `objective_runtime_budget_s`, `objective_telemetry_noise_frac`)
- uncertainty settings (`uncertainty_atm_scale`, `uncertainty_ic_scale`)
- seed and verification command settings

## Run Store Contract

Each run is created under `runs/<run_id>/` and contains:

- `spec.toml`
- `evaluations.arrow` (append-only ledger)
- `state.toml` (checkpoint)
- `best_manifest.toml`
- `report.md`

`run_id` is deterministic from spec + git SHA + dependency lock hash + telemetry hashes + seeds.

## Lifecycle

`prepare -> global_search_quick -> local_refine_full -> robustness_validation -> promote`

Stages are restartable from `state.toml` and `evaluations.arrow`.

## Usage

```text
julia --project=SpaceAGORACalibration.jl -e "using Pkg; Pkg.instantiate()"

# Mock smoke run (package-only)
julia --project=SpaceAGORACalibration.jl SpaceAGORACalibration.jl/bin/run_calibration.jl

# Telemetry calibration driver (process-isolated verify_telemetry backend)
julia --project=SpaceAGORACalibration.jl \
  SpaceAGORACalibration.jl/bin/run_telemetry_tuner.jl \
  --spec=SpaceAGORACalibration.jl/examples/telemetry_hybrid_spec.toml \
  --project=.AGORA \
  --verification-script=src/analysis/verification/TelemetryVerification.jl \
  --profile=quick \
  --plots=0

# Earth-focused telemetry calibration
julia --project=SpaceAGORACalibration.jl \
  SpaceAGORACalibration.jl/bin/run_telemetry_tuner.jl \
  --spec=SpaceAGORACalibration.jl/examples/telemetry_earth_spec.toml \
  --project=.AGORA \
  --verification-script=src/analysis/verification/TelemetryVerification.jl \
  --profile=quick \
  --plots=0

# Earth-only offline study
julia --project=SpaceAGORACalibration.jl \
  SpaceAGORACalibration.jl/bin/run_earth_offline_study.jl \
  --profile=quick \
  --plots=0
```

`CommandBackend` computes the robust scalar objective:
`J(θ) = Σ_s w_s Σ_e [Huber(rmse_km/limit_rmse) + 0.5*Huber(max_abs_km/limit_abs)] + λfail*I(run_failed) + λtime*max(0, runtime/t_budget - 1)`.

Local refinement perturbs continuous parameters only; robustness uses Monte Carlo uncertainty draws and robust ranking (`mean + α*p95 + β*fail_rate`) before full-profile promote confirmation.

Candidate evaluations can run concurrently via `budgets.parallel_evaluations` (override with `SPACEAGORA_CALIBRATION_PARALLEL_EVALUATIONS`).
