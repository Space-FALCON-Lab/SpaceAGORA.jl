# Canonical owner audit: examples, plotting and shared helpers

Runnable examples, plotting and report entry points, and the helpers they share
each have one home. This page names those homes and how each rule is checked.
The [topology contract](canonical_topology_contract.md) and the
[completeness contract](src_completeness_contract.md) cover the rest of the
source tree.

| Responsibility | Owner |
|---|---|
| Example bootstrap | `examples/common.jl` |
| Example entry points | `examples/*.jl` |
| Maintained telemetry plotting | `scripts/plotting/telemetry_orbit_accuracy_plots.jl` |
| Telemetry verification package surface | `src/interfaces/analysis/verification/telemetry_verification.jl` |
| Example helper builders | `src/interfaces/analysis/verification/telemetry_verification/example_support.jl` |
| Runtime serialization locks (`SPICE_LOCK`, `GRAM_LOCK`) | `src/orchestration/simulation/runtime_services.jl` |
| Unfinished subsystem scaffolds | top-level `experimental/` |

## Rules and checks

- Examples and plotting or report scripts have no owner under `src/`.
  `test/gates/ci_path_policy_gate.jl` fails if `src/examples/`,
  `src/interfaces/analysis/plotting/` or `src/interfaces/analysis/reports/`
  exists.
- Telemetry verification resolves its plots from `scripts/plotting/`:
  `src/interfaces/analysis/verification/telemetry_verification/reporting.jl`
  runs `scripts/plotting/telemetry_orbit_accuracy_plots.jl` in a separate Julia
  process and refuses if the script is missing.
- Benchmark plotting launchers forward to `scripts/plotting/`.
  `benchmarks/studies/telemetry_orbit_accuracy_plots.jl` includes the script
  above and fails if it is missing; the path-policy gate requires the launcher.
- `test/gates/ci_architecture_contract_gate.jl` fails unless `RuntimeServices`
  defines `SPICE_LOCK` and `GRAM_LOCK` and `SimulationModel` does not.
- Experimental scaffolds stay outside the package load graph until they own real
  runtime behavior; `test/gates/ci_scaffold_interface_gate.jl` loads them on its
  own and checks that each one refuses with `Not implemented: `.
- An example does not include package source files by relative path; it loads
  SpaceAGORA through `examples/common.jl`. No gate checks this.

## Retained legacy plotting

`scripts/plotting/plot_data.jl` keeps the older plotting dispatcher. It expects
a caller-provided `SimulationModel` module, dictionary arguments and nested
solution structures. No maintained workflow calls it; unit tests load
individual functions from it to pin the saved-data contracts described in
[legacy plotting contracts](../quality/legacy_plotting_contracts.md).

Keep the file until its capabilities and ownership have been reviewed. Modern
examples produce some overlapping plots, but no replacement has been shown for
the costate and switching diagnostics, closed-form comparisons, per-link
histories, attitude and reaction-wheel traces, or torque and inertia histories.
A similar chart, or a saved field that happens to exist, does not show that the
inputs and outputs are equivalent.

## Simulation preparation and viewer payload phases

`run_simulation` stays the public entry point, and
`src/orchestration/simulation/engine/execution.jl` keeps the order of a run:
checkpoint, solve, save and resume. It prepares a run by calling four private
helpers in `src/orchestration/simulation/engine/setup.jl`, in this order:
`_validate_run_models!` (model and configuration checks),
`_prepare_run_runtime!` (run buffers and caches), `_prepare_run_ephemerides!`
(the ephemeris caches) and `_prepare_run_execution_policy!` (execution policy).
The last one calibrates the right-hand-side plan before it starts the in-run
width identification, and resets the density models' histories after every
pre-solve probe. Callback creation, the restore of a resumed run and output
handling stay in `execution.jl`. The helpers add no lifecycle framework of
their own.

`build_viewer_frames` keeps its public signature in the scene module. Three
private helpers in `src/interfaces/analysis/visualization/scene/viewer_payload.jl`
do its work: `_viewer_frame_layout` finds the frame blocks a bundle holds,
`_select_viewer_frame_rows` chooses the rows and the position precision within
the byte budget, keeping mission-boundary rows and the Float32/Float64 rule, and
`_pack_viewer_frame_blocks` packs the selected blocks with the payload's units,
field names, optional-block rules and byte encodings. These are boundaries
inside existing modules, not new public namespaces.
