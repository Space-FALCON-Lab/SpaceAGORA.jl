# Verification contract

What a pull request must pass before it merges into `main`, what a release needs
on top of that, and the coverage floors. The jobs live in
`.github/workflows/julia-ci.yml`. Which of them are required is set by the
repository's `tests` ruleset on GitHub, not by a file in this checkout; the table
below mirrors that ruleset.

## Required checks

| Check | What it runs |
| --- | --- |
| `tests` | Collects the `tests (…)` shards of `test/runtests.jl` (run with `--project=test`) and checks that every shard ran every planned item exactly once (`.github/scripts/ci_shard_verify.py`). It also fails unless the optional companion installation jobs (`optional-hypr-installation`, `optional-mujoco-installation`, `optional-visualizer-installation`, `optional-edg-installation`, `optional-interlink-installation`) pass. |
| `coverage-quality-gate` | Merges the coverage shards (`test/runtests.jl` under `--code-coverage=user`) with the companion installation coverage and runs `test/gates/ci_coverage_quality_gate.jl`. See [Coverage floors](#coverage-floors). |
| `p1-findings-gate` | The contract suite `test/contracts/pr_runtests.jl`, listed in the next section. |
| `ai-review-artifact-gate` | `test/gates/ci_ai_review_artifact_gate.jl`, on pull requests only. See [Review record](#review-record). |
| `clean-depot-smoke` | Instantiates and precompiles the package in an empty Julia depot, then runs `test/smoke/ci_clean_depot_smoke.jl` with `--depwarn=error`. |
| `threaded-smoke` | `test/smoke/ci_threaded_smoke.jl` with four threads and `--depwarn=error`. |
| `readme-command-smoke` | The README's instantiate and examples-setup commands, the quickstart check `test/smoke/ci_readme_quickstart.py`, and `examples/Earth_Thruster_Test.jl`. |
| `example-regression` | `test/integration/examples/runtests.jl`, which runs `test/smoke/ci_examples_regression.jl`, with native GRAM. |
| `example-suite-smoke` | Collects the example smoke shards and checks every example ran once. |
| `telemetry-quick-regression` | Collects the telemetry shards, each running `benchmarks/studies/telemetry_orbit_accuracy_study.jl quick --enforce=true` on its scenarios, and checks every scenario was graded once (`.github/scripts/ci_telemetry_verify.py`). With `--enforce=true`, any row that misses its `tolerances_quick` thresholds in `test/telemetry_benchmark_manifest.toml` fails the run. |
| `documentation-release-gate` | `docs/make.jl --no-deploy`, then `test/gates/ci_public_api_surface_gate.jl` and `test/gates/ci_hpc_extensibility_docs_gate.jl`. |
| `docs` | The documentation build and deploy in `.github/workflows/docs.yml`. |

The `contributor-kit` job, the Windows install workflow and the kit demonstration
workflow also run on some pull requests but are not required; they are described
at the end of this page.

## Contract gates on every pull request

`test/contracts/pr_runtests.jl` includes these gates, in this order. The
`p1-findings-gate` job runs it with the test environment, and
`test/contracts/nightly_runtests.jl` includes the same file, so the nightly run
checks the same list.

| Gate | What it checks |
| --- | --- |
| `test/gates/ci_p1_findings_gate.jl` | No P1 marker outside the allowlist. See [P1 findings](#p1-findings). |
| `test/gates/ci_naming_contract_gate.jl` | `lower_snake_case.jl` file names in the scanned folders; the retired module names `ref_sys` and `PhysicalModel` stay gone. See the [API and naming contract](api_naming_contract.md). |
| `test/gates/ci_architecture_contract_gate.jl` | Module load order (`RuntimeServices`, then `SimulationModel`, then `SimulationEngine`), lock ownership, the typed `run_simulation` boundary, `TelemetryVerification` calling `SimulationEngine.run_simulation` without including engine files, no `ENV` reads in the engine outside `adapters/from_env.jl`, and the split runtime-analysis and parallelization benchmark files. |
| `test/gates/ci_path_policy_gate.jl` | The four source groups and their domain folders, retired paths staying absent (including the old benchmark and study wrappers under `test/`), required owner paths, and forbidden path tokens. |
| `test/gates/ci_thin_entry_files_gate.jl` | Each listed entry file stays at or under 80 non-comment lines and defines no functions or structs. |
| `test/gates/ci_no_artifact_files_gate.jl` | No tracked `*.cov`, `.DS_Store`, `output/`, `docs/build/`, `docs/site/` or `docs/src/generated/` files, and no absolute home-directory paths (under macOS `Users`, Linux `home` or Windows `Users` folders) in tracked text files. |
| `test/gates/ci_tracked_feather_registry_gate.jl` | Every tracked feather file is named in `test/telemetry_benchmark_manifest.toml` or lies under an asset in `data/assets_manifest.toml` that is required and not lab-private. |
| `test/gates/ci_stale_deps_gate.jl` | Aqua's stale-dependency check: every `[deps]` entry of `Project.toml` is loaded by the package. No dependency is exempt; a package that only the tests, the docs or the examples use belongs in `test/Project.toml`, `docs/Project.toml` or `examples/Project.toml`. Aqua comes from the test environment, so the gate runs with `--project=test`. |
| `test/gates/ci_vehicle_structure_boundary_gate.jl` | The structure owner files under `src/models/vehicle/structure/` exist, and nothing under `src/models/vehicle/spacecraft/` defines a structure-owner function such as `get_COM` or `get_inertia_tensor`. |
| `test/gates/ci_src_completeness_contract_gate.jl` | No silent `*_stub() = nothing` stubs or placeholder tokens in `src/`; each listed aggregator carries the `# Canonical aggregator: no behavior ownership.` header and defines no functions or types. |
| `test/gates/ci_no_legacy_include_chains_gate.jl` | No `__legacy_` token and no three-level `"..", "..", ".."` include in the scanned source folders, apart from a listed set of composition owners. |
| `test/gates/ci_gnc_aerobraking_boundary_gate.jl` | Aerobraking guidance does not include control source, and GNC source does not use `DynamicEffectors` directly. |
| `test/gates/ci_no_guidance_control_cross_include_gate.jl` | Guidance, planning, companion adapters, contracts and the analysis planning folder do not `include` control files. |
| `test/gates/ci_rotational_ownership_gate.jl` | The rotational equations appear only under `src/models/dynamics/rotational/`. |
| `test/gates/ci_scaffold_interface_gate.jl` | The experimental resources, constellation-network and estimation scaffolds throw `ErrorException` messages that start with `Not implemented: `. |
| `test/gates/ci_typed_config_equivalence_gate.jl` | A typed `SimulationEngineConfig` produces the expected environment overrides and solver policy. |
| `test/gates/ci_rhs_calibration_gate.jl` | The RHS calibration functions in `src/orchestration/simulation/engine/rhs_calibration.jl` exist and are wired into setup and execution. |
| `test/gates/ci_rhs_parallel_route_gate.jl` | The RHS parallel-route helpers, routing predicates and policy counters exist and are used where routing is decided. |
| `test/gates/ci_docs_contract_gate.jl` | Markers in the pull request template, the documentation policy and `docs/make.jl`; contributor material stays outside `docs/src/`; the listed contributor contracts exist, this page included. |
| `test/gates/ci_public_api_surface_gate.jl` | The exports of `SpaceAGORA` equal the list in `docs/public_api_symbols.jl` (the `SpaceAGORA.Conformance` exports are compared separately), and internal modules such as `SimulationModel` and `SimulationEngine` are not exported. |
| `test/gates/ci_hpc_extensibility_docs_gate.jl` | The extension templates, starters and HPC scripts exist; the extensibility and HPC pages keep their markers; and every snippet region of the scripts in `examples/gnc_algorithm_env/` and `examples/walkthrough_env/` appears unchanged as a code block of the page that shows it (`gnc_feedback.md`, `gnc_planning.md`, `custom_effector.md`, `custom_actuator.md`), with each script listed and run by its environment's `smoke.jl`. |
| `test/gates/ci_extension_hook_reference_gate.jl` | Every exported function in the hook sections of `docs/public_api_symbols.jl` is a row of `docs/extension_hooks.jl` or a listed non-hook. Each row's signature is its docstring's first line, its arity and fallback agree with the method table, its page names it, deprecated rows are consistent, and current rows use `spacecraft_index`. Listed known gaps fail once they no longer apply. |
| `test/gates/ci_gnc_typed_command_boundary_gate.jl` | GNC reaches actuators through the typed command boundary (`PropulsiveManeuverCommand`, `AerobrakingControlCommand`), and the algorithm folders stay in the scan. |
| `test/gates/ci_gnc_random_stream_gate.jl` | Code under `src/models/gnc/`, the in-repository HYPR adapter included, draws random numbers only from an explicit stream argument: no `default_rng()`, `seed!` or stream-less `rand`/`randn`. |
| `test/gates/ci_project_compat_gate.jl` | Every non-standard-library `[deps]` entry has a `[compat]` bound. |
| `test/gates/ci_sysimage_package_lists_gate.jl` | Every package named in the `FULL_PACKAGES` and `CI_PACKAGES` lists of `build/build_sysimage.jl` is `:SpaceAGORA` or a root `[deps]` entry. It parses the two files and builds nothing. |
| `test/gates/ci_runtime_legacy_surface_gate.jl` | The active runtime does not regain prohibited legacy entry points, such as dictionary-typed `args` or the old `gravity_*` functions. |
| `test/gates/ci_translational_ownership_gate.jl` | The translational equations appear only under `src/models/dynamics/translational/`. |
| `test/gates/ci_runtime_any_hotpath_gate.jl` | No new explicit `Vector{Any}` or `Dict{Symbol, Any}` in the runtime hot paths, whatever the spacing or qualification; each listed heterogeneous boundary occurs exactly once, in its declared owner. |
| `test/gates/ci_ppc_mode_profile_parity_gate.jl` | The parallelization benchmark's base modes, scheduler and profile identities agree with the shipped routing profiles they name. |
| `test/gates/ci_runtime_analysis_copy_overhead_gate.jl` | On synthetic rows: the runtime-analysis summary and report copy-overhead values, successful-run denominators, the failed, empty and zero cases, ordering and input preservation. It runs no simulation and measures nothing. |
| `test/gates/ci_removed_hook_names_gate.jl` | No method, declaration or other use of an extension hook name removed in 0.5.0, and no five-argument `wrench` method or call, in `src/`, `ext/`, `templates/` or `examples/` (found by parsing); no file of the contributor starters names a removed hook. |

The workload-coverage check for the parallelization paper, which confirms that
the research workload names every case and mode its phases run, lives with
SpaceAGORAParallelizationStudies, not here. The mode-profile parity and
runtime-analysis gates above stay required in this repository.

## Coverage floors

`test/gates/ci_coverage_quality_gate.jl` reads the `.cov` files under `src/` and
the pinned EDG source `packages/gnc/SpaceAGORAEDG.jl/src/`, keeps those written
within `SPACEAGORA_COVERAGE_WINDOW_SECONDS` (default 900) of the newest one, and
fails on any of these:

- overall coverage below 90.0% (`MIN_MAIN_OVERALL`, overridable with
  `SPACEAGORA_COVERAGE_MIN_OVERALL`);
- any file below 80.0% (`MIN_MAIN_FILE`, overridable with
  `SPACEAGORA_COVERAGE_MIN_FILE`), except the per-file overrides below;
- any critical file below 90.0%, or a critical file with no coverage data;
- an EDG source file with no coverage data.

| File | Floor |
| --- | --- |
| `src/orchestration/simulation/engine/execution.jl` | 90.0% (critical) |
| `src/models/gnc/control/maneuvers/propulsive_maneuvers.jl` | 90.0% (critical) |
| `src/foundation/core/interfaces/reference_system.jl` | 90.0% (critical) |
| `src/orchestration/simulation/engine/adapters/from_env.jl` | 70.0% |
| `src/orchestration/simulation/engine/dynamics_rhs.jl` | 70.0% |

The gate refuses to start if any file named in those two tables no longer exists,
so a moved file cannot fall back to the default floor.

Lines inside PrecompileTools `@setup_workload` and `@compile_workload` calls run
only while the package precompiles, never under `--code-coverage=user`, so the
gate drops them from the executable-line count. It finds them by parsing each
source file; a file that does not parse fails the gate. Each affected file is
printed as a `precompile_only_excluded` line. Functions called from those blocks
still count. `test/unit/quality/coverage_gate_precompile_exclusion_tests.jl`
tests this exclusion.

## Native GRAM probe reporting

The coverage shards run the threaded probe child `test/probes/coverage_threaded_probes.jl`
from the "Coverage Threaded Probe Driver" testset in
`test/suites/02_callbacks_parallel_and_smoke_tests.jl`. The protocol is in
`test/helpers/native_probe_reporting.jl`. The child prints
`spaceagora_native_gram_probes=completed` only after both the native
construction probes and the density-service probes return. When native GRAM is
not set up, it prints `spaceagora_native_gram_probes=skipped`. The parent
forwards the status even when the child passes.

The coverage shards set `SPACEAGORA_REQUIRE_NATIVE_GRAM_PROBES=1` (the only
accepted values are `0` and `1`). With it set, a skip, a missing or malformed
status, a failed child, or a run without coverage enabled fails the test. This
makes native execution visible without relying on an assertion count. The
completion marker covers these named probes only, not every native GRAM scenario.

## P1 findings

`test/gates/ci_p1_findings_gate.jl` scans `.jl`, `.md`, `.txt`, `.yml`, `.yaml`
and `.toml` files under `src/` and `test/`, plus the root `README.md`, for lines
containing `[P1]`, `priority = 1` or `severity: 1` (also `severity=P1`; the last
two ignore case).

- A marker not in `test/p1_findings_allowlist.txt` fails the gate.
- An allowlist entry is `path::exact_line`, where `exact_line` is the line with
  surrounding whitespace removed.
- An entry may end with `# owner=<name> opened=YYYY-MM-DD expires=YYYY-MM-DD`;
  each field is optional.
- An entry past its `expires` date fails the gate. Entries that no longer match
  any line are printed but do not fail it.

## Review record

`test/gates/ci_ai_review_artifact_gate.jl` applies to a pull request that adds,
modifies or renames any `src/**/*.jl` file, compared with `origin/<base>`
(falling back to `HEAD~1` when the base cannot be read). Such a pull request
needs `test/ai_reviews/PR_<number>.md` with these `##` sections, none of them
empty:

1. `Scope`
2. `Changed Files`
3. `Findings`
4. `P1 Assessment`
5. `Tests Added/Updated`
6. `Residual Risk`

`Changed Files` must name every changed `src/*.jl` file by its full path. A pull
request that changes no file under `src/` needs no record.

## Nightly stress and releases

`.github/workflows/nightly-stress.yml` runs every day at 07:00 UTC, on manual
dispatch, and on pull requests that change the dynamics, environment, GNC,
simulation, verification or telemetry paths listed in the workflow. It runs the
contract suite through `test/contracts/nightly_runtests.jl`, the threaded smoke,
the example regression, the stress suite `test/stress/runtests.jl` and the
full-profile telemetry regression, one job per scenario.

Pushing a `v*` tag runs `.github/workflows/release-tag-gate.yml`, which fails
unless the latest successful nightly-stress run on `main` finished within the
last 36 hours (`SPACEAGORA_NIGHTLY_MAX_AGE_HOURS`).

## Windows install check

`.github/workflows/windows-install.yml` runs the documented install path and the
README quickstart smoke on `windows-latest` with Julia 1.12, weekly, on manual
dispatch, and on pull requests that change the install path or the pages that
document its commands. It is not a required check. A non-blocking job in the same
workflow records how Julia 1.13 is refused or fails on that path.

## Contributor CI workflow

`.github/workflows/contributor-ci.yml` is a reusable workflow (`on: workflow_call`)
for repositories built on SpaceAGORA. The caller passes the same full commit SHA
as the `uses:` ref and as `spaceagora-rev`, and records it in
`SpaceAGORASource.toml` beside its `Project.toml`, which has no `[sources]` entry
for SpaceAGORA. The workflow refuses a mismatch. It tests the repository in a
separate environment built by `scripts/contributor_environment.jl`, with that core
commit developed by path, and certifies interface conformance only.

The `contributor-kit` job in `julia-ci.yml` runs the workflow on each starter in
`templates/starters/` (project, package and companion), generated from the commit
under test by `scripts/new_contributor_repository.jl`, with this checkout as the
core. It is not a required check.

`.github/workflows/contributor-kit-demonstration.yml` runs the same workflow on
the kit's demonstration repository, `Space-FALCON-Lab/SpaceAGORAKitDemo.jl`
(input `repository`): every night against `main`, and in each release pull
request (a branch named `release-*`) against its merge commit. It first moves the
demonstration's pin to the commit under test with
`scripts/contributor_environment.jl repin`. It is not a required check either,
but a failure is resolved before the release is tagged.
