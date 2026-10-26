# SpaceAGORA Verification Contract

## Objective
This contract defines merge and release quality gates for SpaceAGORA architecture and runtime changes.

## Required PR Gates
1. `tests` (`test/runtests.jl`)
2. `coverage-quality-gate` (`test/gates/ci_coverage_quality_gate.jl`)
3. `ai-review-artifact-gate` (`test/gates/ci_ai_review_artifact_gate.jl`)
4. The contract gates run by `test/contracts/pr_runtests.jl`:
   - `test/gates/ci_p1_findings_gate.jl`
   - `test/gates/ci_naming_contract_gate.jl`
   - `test/gates/ci_architecture_contract_gate.jl`
   - `test/gates/ci_path_policy_gate.jl`
   - `test/gates/ci_thin_entry_files_gate.jl`
   - `test/gates/ci_no_artifact_files_gate.jl`
   - `test/gates/ci_stale_deps_gate.jl` (Aqua stale-dependency check: every `[deps]` entry is loaded by the package; the packages only the examples and tests still load are an explicit ignore list, and the test and docs tools live in `test/Project.toml` and `docs/Project.toml`)
   - `test/gates/ci_tracked_feather_registry_gate.jl` (every tracked feather is named in `test/telemetry_benchmark_manifest.toml` or lies under an in-tree required, non lab-private `data/assets_manifest.toml` asset)
   - `test/gates/ci_vehicle_structure_boundary_gate.jl`
   - `test/gates/ci_src_completeness_contract_gate.jl`
   - `test/gates/ci_no_legacy_include_chains_gate.jl`
   - `test/gates/ci_gnc_aerobraking_boundary_gate.jl`
   - `test/gates/ci_no_guidance_control_cross_include_gate.jl`
   - `test/gates/ci_rotational_ownership_gate.jl`
   - `test/gates/ci_scaffold_interface_gate.jl`
   - `test/gates/ci_typed_config_equivalence_gate.jl`
   - `test/gates/ci_rhs_calibration_gate.jl`
   - `test/gates/ci_rhs_parallel_route_gate.jl`
   - `test/gates/ci_docs_contract_gate.jl`
   - `test/gates/ci_public_api_surface_gate.jl`
   - `test/gates/ci_hpc_extensibility_docs_gate.jl` (the extension templates, starters and HPC scripts exist, the extensibility and HPC pages keep their markers, and every snippet region of the scripts in `examples/gnc_algorithm_env/` and `examples/walkthrough_env/` appears unchanged as a code block of the page that shows it (`gnc_feedback.md`, `gnc_planning.md`, `custom_effector.md`, `custom_actuator.md`), with each script listed and run by its environment's `smoke.jl`)
   - `test/gates/ci_extension_hook_reference_gate.jl` (every exported function in the hook sections of `docs/public_api_symbols.jl` is a row of `docs/extension_hooks.jl` or a listed non-hook; each row's signature is its docstring's first line, its arity and fallback agree with the method table, its page names it, deprecated rows are consistent and current rows use `spacecraft_index`; listed known gaps fail once they no longer apply)
   - `test/gates/ci_gnc_typed_command_boundary_gate.jl` (GNC reaches actuators through the typed command boundary; algorithm directories stay in scan scope)
   - `test/gates/ci_gnc_random_stream_gate.jl` (GNC and the in-repo HYPR adapter draw random numbers only from an explicit stream argument: no `default_rng()`, `seed!` or stream-less `rand`/`randn`)
   - `test/gates/ci_project_compat_gate.jl` (every non-standard-library `[deps]` entry has a `[compat]` bound)
   - `test/gates/ci_sysimage_package_lists_gate.jl` (every package named in the sysimage builder's `FULL_PACKAGES` and `CI_PACKAGES` is `:SpaceAGORA` or a root `[deps]` entry; parses `build/build_sysimage.jl` and `Project.toml`, builds nothing)
   - `test/gates/ci_runtime_legacy_surface_gate.jl` (the active runtime does not regain prohibited legacy entry surfaces)
   - `test/gates/ci_translational_ownership_gate.jl` (translational equations stay with their owner directory)
   - `test/gates/ci_runtime_any_hotpath_gate.jl` (no new explicit `Vector{Any}` or `Dict{Symbol, Any}` in runtime hot paths, including whitespace and qualified spellings; each inventoried heterogeneous boundary must occur exactly once in its declared owner)
   - `test/gates/ci_ppc_mode_profile_parity_gate.jl` (the parallelization benchmark's evaluated mode values, scheduler and declared profile identity agree with the shipped routing profiles it names)
   - `test/gates/ci_runtime_analysis_copy_overhead_gate.jl` (synthetic rows check actual summary/report copy-overhead values, successful-run denominators, failed/empty/zero cases, ordering and input preservation; no simulation or performance measurement)
   - `test/gates/ci_removed_hook_names_gate.jl` (no method, declaration or other use of an extension hook name removed in 0.5.0, and no five-argument `wrench` method or call, in `src/`, `ext/`, `templates/` or `examples/`, found by parsing; the contributor starters do not name a removed hook in any file)

The paper workload coverage gate belongs to
SpaceAGORAParallelizationStudies.
It checks that the research workload names every case and mode its paper phases
run. The core mode-profile parity and runtime-analysis gates above remain
required.

Additional required smoke checks for this migration track:
1. `test/smoke/ci_clean_depot_smoke.jl`
2. `test/smoke/ci_threaded_smoke.jl`
3. `test/smoke/ci_examples_regression.jl`
4. `benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true`

Telemetry threshold failures are blocking in this cleanup track.

## Coverage Policy
Coverage enforcement is implemented by `test/gates/ci_coverage_quality_gate.jl`.
Lines inside PrecompileTools `@setup_workload` / `@compile_workload` macro calls
run only while the package precompiles, never under `--code-coverage=user`, so
the gate drops them from the executable-line count (found by parsing each source
file; a file that fails to parse fails the gate) and prints each affected file
as a `precompile_only_excluded` line; functions called from those blocks still count.

Main thresholds:
1. Main overall coverage: `>= 90.0%`
2. Main per-file coverage: `>= 80.0%` (unless explicit override exists)

Critical-file thresholds (`>= 90.0%`):
1. `src/orchestration/simulation/engine/execution.jl`
2. `src/models/gnc/control/maneuvers/propulsive_maneuvers.jl`
3. `src/foundation/core/interfaces/reference_system.jl`

Main per-file overrides:
1. `src/orchestration/simulation/engine/adapters/from_env.jl` => `>= 70.0%`
2. `src/orchestration/simulation/engine/dynamics_rhs.jl` => `>= 70.0%`

## Native GRAM probe reporting

The threaded coverage child prints `spaceagora_native_gram_probes=completed`
only after both native construction and density-worker probes return successfully.
Its parent forwards that status even when the child passes. When native setup
is unavailable, a developer run prints `spaceagora_native_gram_probes=skipped`.

The native-enabled coverage-quality job sets
`SPACEAGORA_REQUIRE_NATIVE_GRAM_PROBES=1`. A skip, missing or malformed status,
failed child, or disabled coverage then fails the test. This makes native execution
visible without depending on a fixed assertion count. A completion marker covers
these named probes, not every native-GRAM scenario or a performance comparison.

## Architecture/Dependency Policy
1. `SpaceAGORA.run_simulation` forwards to `SimulationEngine.run_simulation`.
2. `TelemetryVerification` calls `SimulationEngine.run_simulation`.
3. `TelemetryVerification` does not include retired simulation-execution wrappers directly.
4. Engine internals do not read `ENV` directly outside adapter files.
5. Retired wrapper files remain absent (benchmark/study wrappers under `test/`).

Items 1-4 are enforced by `test/gates/ci_architecture_contract_gate.jl`; item 5 by `test/gates/ci_path_policy_gate.jl`.

## P1 Findings Policy
Enforced by `test/gates/ci_p1_findings_gate.jl`.

1. New unallowlisted P1 markers fail CI.
2. Allowlist format: `path::exact_line`.
3. Optional metadata fields: `owner`, `opened`, `expires`.
4. Expired allowlist entries fail CI.

## AI Artifact Policy
Enforced by `test/gates/ci_ai_review_artifact_gate.jl`.

Required path: `test/ai_reviews/PR_<number>.md`.
Required sections:
1. `Scope`
2. `Changed Files`
3. `Findings`
4. `P1 Assessment`
5. `Tests Added/Updated`
6. `Residual Risk`

## Nightly Policy
Nightly stress checks remain required for release readiness and run through `.github/workflows/nightly-stress.yml`.

## Windows Install Check
`.github/workflows/windows-install.yml` runs the documented install path and the README quickstart smoke on `windows-latest` with Julia 1.12, weekly, on manual dispatch and on pull requests that change the install path or the pages that document its commands. It is not a required check. A non-blocking job in the same workflow records how Julia 1.13 is refused or fails on that path.

## Contributor CI Workflow
`.github/workflows/contributor-ci.yml` is a reusable workflow (`on: workflow_call`) for repositories built on SpaceAGORA. The caller passes the same full commit SHA as the `uses:` ref and as `spaceagora-rev`, and records it in `SpaceAGORASource.toml` beside its `Project.toml`, which has no `[sources]` entry for SpaceAGORA. The workflow refuses a mismatch, tests the repository in a separate environment built by `scripts/contributor_environment.jl` with that core developed by path, and certifies interface conformance only. The `contributor-kit` job in `julia-ci.yml` runs it on each starter of `templates/starters/` (project, package and companion), generated from the commit under test by `scripts/new_contributor_repository.jl`, with this checkout as the core. It is not a required check. `.github/workflows/contributor-kit-demonstration.yml` runs the same workflow on the kit's demonstration repository, `Space-FALCON-Lab/SpaceAGORAKitDemo.jl` (input `repository`), every night against `main` and in each release pull request (a branch named `release-*`) against its merge commit, after `scripts/contributor_environment.jl repin` moves the demonstration's pin to the commit under test. It is not a required check either; a failure is resolved before the release is tagged.
