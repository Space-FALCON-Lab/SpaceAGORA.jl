# SpaceAGORA `src/` Canonical Owner Audit

This audit records the ownership locations for runnable examples, maintained
plotting/report entrypoints and shared helpers, together with retained legacy
plotting code awaiting a capability decision.

## Cleanup Status

1. Runnable example entrypoints are owned by top-level `examples/`.
2. Plotting/report scripts are owned by top-level `scripts/plotting/`.
3. Shared example helper builders live in `src/analysis/verification/telemetry_verification/example_support.jl`.
4. Shared runtime lock ownership lives in `src/simulation/runtime_services.jl`.
5. Unfinished subsystem scaffolds are quarantined under top-level
   `experimental/*`, not under `src/*`.
6. Only the canonical roots above are valid ownership locations for runnable
   examples and plotting/report entrypoints.
7. Forbidden paths are enforced by CI gates. Retained legacy plotting code is
   described separately below.

## Canonical Owners

1. Example bootstrap: `examples/common.jl`
2. Example entrypoints: `examples/*.jl`
3. Maintained telemetry plotting: `scripts/plotting/telemetry_orbit_accuracy_plots.jl`,
   invoked by `src/analysis/verification/telemetry_verification/reporting.jl`.
4. Example helper builders: `src/analysis/verification/telemetry_verification/example_support.jl`
5. Runtime serialization locks: `src/simulation/runtime_services.jl`
6. Telemetry verification package surface: `src/analysis/verification/telemetry_verification.jl`

## Retained legacy plotting

`scripts/plotting/plot_data.jl` retains the older plotting dispatcher. It
expects a caller-provided `SimulationModel` module, dictionary arguments and
nested solution structures. No tracked Julia caller was located in this source
review; current workflow support has not been established.

Preserve the file pending capability and ownership review. Modern examples
provide some overlapping plots, but replacement coverage has not been established
for costate/switching diagnostics, closed-form comparisons, per-link histories,
attitude and reaction-wheel traces, or torque/inertia histories. Similar charts
or available saved fields alone do not establish equivalent inputs and outputs.

The telemetry caller above establishes a maintained source route. Runtime and
plotting-equivalence acceptance remain separate from this source review.

## Verification Notes

1. No example file should include package internals by relative path.
2. Benchmark plotting launchers should forward to `scripts/plotting/`.
3. Telemetry verification reporting should resolve plotting from `scripts/plotting/`.
4. Experimental scaffolds must remain outside the package load graph until they
   own real runtime behavior.
