# SpaceAGORA `src/` Canonical Owner Audit

This audit records the ownership locations for runnable examples, maintained
plotting/report entrypoints and shared helpers, together with retained legacy
plotting code awaiting a capability decision.

## Cleanup Status

1. Runnable example entrypoints are owned by top-level `examples/`.
2. Plotting/report scripts are owned by top-level `scripts/plotting/`.
3. Shared example helper builders live in `src/interfaces/analysis/verification/telemetry_verification/example_support.jl`.
4. Shared runtime lock ownership lives in `src/orchestration/simulation/runtime_services.jl`.
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
   invoked by `src/interfaces/analysis/verification/telemetry_verification/reporting.jl`.
4. Example helper builders: `src/interfaces/analysis/verification/telemetry_verification/example_support.jl`
5. Runtime serialization locks: `src/orchestration/simulation/runtime_services.jl`
6. Telemetry verification package surface: `src/interfaces/analysis/verification/telemetry_verification.jl`

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

## Simulation preparation and viewer payload phases

The execution coordinator retains the public `run_simulation` entrypoint and
its checkpoint/solve/save/resume sequence. Ordered private helpers in
`src/orchestration/simulation/engine/setup.jl` own model validation,
runtime/cache initialization, ephemerides preparation and execution-policy
preparation. The policy phase still calibrates before starting in-run width
identification and resets density histories after all pre-solve probes.
Callback creation, resume restoration and output handling retain their existing
positions in `execution.jl`; the helpers introduce no new lifecycle framework.

`build_viewer_frames` retains its public signature in the scene module.
Private helpers in `viewer_payload.jl` discover available frame blocks, select
rows and precision under the byte budget, and pack the selected blocks.
The layout carries existing counts and offsets; selection preserves mission
boundary rows and the Float32/Float64 policy. Packing retains existing units,
field names, optional-block rules and byte encodings. These are responsibility
boundaries within existing modules, not new public namespaces.
