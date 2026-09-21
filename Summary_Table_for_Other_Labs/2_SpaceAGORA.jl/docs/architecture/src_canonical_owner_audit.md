# SpaceAGORA `src/` Canonical Owner Audit

This audit records the post-cleanup ownership boundary for runnable examples,
plotting scripts, and helper ownership.

## Cleanup Status

1. Runnable example entrypoints are owned by top-level `examples/`.
2. Plotting/report scripts are owned by top-level `scripts/plotting/`.
3. Shared example helper builders live in `src/analysis/verification/telemetry_verification/example_support.jl`.
4. Shared runtime lock ownership lives in `src/simulation/runtime_services.jl`.
5. Unfinished subsystem scaffolds are quarantined under top-level
   `experimental/*`, not under `src/*`.
6. Only the canonical roots above are valid ownership locations for runnable
   examples and plotting/report entrypoints.
7. Legacy or forbidden paths are enforced by CI gates rather than repeated in
   this audit.

## Canonical Owners

1. Example bootstrap: `examples/common.jl`
2. Example entrypoints: `examples/*.jl`
3. Plotting scripts: `scripts/plotting/plot_data.jl` and `scripts/plotting/telemetry_orbit_accuracy_plots.jl`
4. Example helper builders: `src/analysis/verification/telemetry_verification/example_support.jl`
5. Runtime serialization locks: `src/simulation/runtime_services.jl`
6. Telemetry verification package surface: `src/analysis/verification/telemetry_verification.jl`

## Verification Notes

1. No example file should include package internals by relative path.
2. Benchmark plotting launchers should forward to `scripts/plotting/`.
3. Telemetry verification reporting should resolve plotting from `scripts/plotting/`.
4. Experimental scaffolds must remain outside the package load graph until they
   own real runtime behavior.
