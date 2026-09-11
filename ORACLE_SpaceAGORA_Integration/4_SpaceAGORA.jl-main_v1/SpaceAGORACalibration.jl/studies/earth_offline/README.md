# Earth Offline Study

This study is an Earth-only telemetry calibration campaign designed for reproducible offline runs.

## Files

- `spec_earth_offline.toml`: immutable study spec (search space, budgets, robustness).
- `spec_earth_offline_smoke.toml`: reduced-budget smoke spec for fast validation.
- `manifest_earth_offline.toml`: Earth-only verification manifest.

## Run (from repo root)

PowerShell:

```powershell
$env:JULIA_PKG_OFFLINE = "true"
julia --project=. SpaceAGORACalibration.jl/bin/run_earth_offline_study.jl --profile=quick --plots=0
```

`cmd.exe`:

```bat
set JULIA_PKG_OFFLINE=true
julia --project=. SpaceAGORACalibration.jl/bin/run_earth_offline_study.jl --profile=quick --plots=0
```

POSIX shells:

```bash
export JULIA_PKG_OFFLINE=true
julia --project=. SpaceAGORACalibration.jl/bin/run_earth_offline_study.jl --profile=quick --plots=0
```

## Smoke run (fast sanity)

PowerShell:

```powershell
$env:JULIA_PKG_OFFLINE = "true"
julia --project=. SpaceAGORACalibration.jl/bin/run_earth_offline_study.jl --spec=SpaceAGORACalibration.jl/studies/earth_offline/spec_earth_offline_smoke.toml --profile=quick --plots=0
```

`cmd.exe`:

```bat
set JULIA_PKG_OFFLINE=true
julia --project=. SpaceAGORACalibration.jl/bin/run_earth_offline_study.jl --spec=SpaceAGORACalibration.jl/studies/earth_offline/spec_earth_offline_smoke.toml --profile=quick --plots=0
```

POSIX shells:

```bash
export JULIA_PKG_OFFLINE=true
julia --project=. SpaceAGORACalibration.jl/bin/run_earth_offline_study.jl --spec=SpaceAGORACalibration.jl/studies/earth_offline/spec_earth_offline_smoke.toml --profile=quick --plots=0
```

The engine stage map remains:

- `global_search_quick` uses quick profile.
- `local_refine_full`, `robustness_validation`, and `promote` use full profile.

## Resume behavior

Run the same command again to resume from checkpoint (`state.toml`) for the same deterministic `run_id`.

## Artifacts

Outputs are written under:

- `output/calibration/earth_offline/runs/<run_id>/`
- `spec.toml`, `evaluations.arrow`, `state.toml`, `best_manifest.toml`, `report.md`

This study enables candidate-level outer parallelism with `budgets.parallel_evaluations` in the spec.
You can override it with the `SPACEAGORA_CALIBRATION_PARALLEL_EVALUATIONS`
environment variable in your shell before launching the study.

Replica artifacts during robustness stage are isolated as:

- `candidate_XXXX_robustness_validation_rNNN/`

## Optional overrides

- `--spec=...`: use a copy/tweak of this study spec.
- `--manifest=...`: override manifest path for what-if runs.
- `--enforce=1`: fail verification process when tolerance checks fail.
