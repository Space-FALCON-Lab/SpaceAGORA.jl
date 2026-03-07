# SpaceAGORA `src/` Canonical Owner Audit

This audit records the post-cleanup ownership boundary for runnable examples, plotting scripts, and retired helper scaffolding.

## Cleanup Status

1. Runnable example entrypoints are owned by top-level `examples/`.
2. Plotting/report scripts are owned by top-level `scripts/plotting/`.
3. Shared example helper builders live in `src/analysis/verification/telemetry_verification/example_support.jl`.
4. `src/core/utils/typed_example_utils.jl` is retired and must not exist.
5. `src/examples/`, `src/analysis/plotting/`, `src/analysis/reports/`, `src/parallel/state/`, `src/parallel/tuning/`, and `src/vehicle/sensors/` are retired and must not exist.

## Canonical Owners

1. Example bootstrap: `examples/common.jl`
2. Example entrypoints: `examples/*.jl`
3. Plotting scripts: `scripts/plotting/plot_data.jl` and `scripts/plotting/telemetry_orbit_accuracy_plots.jl`
4. Example helper builders: `src/analysis/verification/telemetry_verification/example_support.jl`
5. Telemetry verification package surface: `src/analysis/verification/telemetry_verification.jl`

## Verification Notes

1. No example file should include package internals by relative path.
2. Benchmark plotting launchers should forward to `scripts/plotting/`.
3. Telemetry verification reporting should resolve plotting from `scripts/plotting/`.
