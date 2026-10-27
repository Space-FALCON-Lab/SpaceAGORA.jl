# Tests

## Prepare

```bash
julia --startup-file=no test/setup.jl
julia --startup-file=no examples/setup.jl
```

The test environment (`test/Project.toml`) holds the test tools, Aqua among
them, and every package the harness imports, with SpaceAGORA from this
checkout; the root project does not include Aqua. Suites 01 and 03 load the
package source into sandbox modules, so the test environment also declares all
of SpaceAGORA's own dependencies, and `test/setup.jl` refuses a test project
that misses any. `test/setup.jl` regenerates the ignored `test/Manifest.toml`
from the root manifest every time, so the tests use the versions the root
project pins; run it again after the root manifest changes. `Pkg.test()` from
the root project uses the same test project.

The harness requires Aqua and stops if it would test another checkout's
SpaceAGORA. The tests also start example and command-line child processes,
which run in the examples environment (`examples/Project.toml`, prepared by
`examples/setup.jl`). The EDG and HYPR providers stay in separate
environments: point `SPACEAGORA_EDG_TEST_ENV` and `SPACEAGORA_HYPR_TEST_ENV` at
environments prepared with `scripts/optional_packages.jl setup edg <directory>`
and `scripts/optional_packages.jl setup hypr <directory>`. The older
`scripts/setup_edg.jl` and `scripts/setup_hypr.jl` commands are compatibility
shortcuts to those companion-owned installers.

## Run

```bash
julia --project=test test/runtests.jl
julia --project=test test/unit/runtests.jl
julia --project=test test/contracts/pr_runtests.jl
julia --project=test test/contracts/nightly_runtests.jl
julia --project=test test/contracts/runtests.jl
julia --threads=2 --project=. test/smoke/runtests.jl
julia --project=. test/stress/runtests.jl
julia --project=. test/coverage/runtests.jl
julia --startup-file=no test/grid_atmosphere/setup.jl
julia --startup-file=no --threads=4 --project=test/grid_atmosphere test/grid_atmosphere/runtests.jl
```

The smoke runner needs at least two threads because it includes the threaded
smoke. `test/coverage/runtests.jl` grades coverage data, so run the tests with
`--code-coverage=user` first. CI runs most smoke files one at a time rather
than through `test/smoke/runtests.jl`.

## Layout

- `runtests.jl` is the default entry point. It hands over to
  `integration/runtests.jl`, which loads the package (`using SpaceAGORA`),
  binds `SimulationModel`, `SimulationEngine`, `TelemetryVerification`, the
  parallel modules and a few frame helpers as aliases, and includes the
  numbered suites in `suites/` (00 to 10) into one shared scope with mock
  models and builder helpers. Because everything comes from the loaded package,
  the process holds one copy of each module. The harness includes no file from
  `src/` except `src/orchestration/mission/operations/maneuver_plans.jl`, which
  the package does not load. Suites 01 and 03 deliberately load source files
  into throwaway sandbox modules to test standalone loading; keep those. The
  reference-system sandboxes use `helpers/reference_system_sandbox.jl`, which
  loads the Geodesy and OrbitalGeometry sources before the frame file without
  importing SpaceAGORA or `RuntimeServices`.
- `unit/` holds standalone test files (each starts with `using SpaceAGORA`)
  run by `unit/runtests.jl`. They cover parallel execution, simulation,
  dynamics, environment, GNC, I/O, vehicles, visualization, quality checks,
  gates, the CLI, missions, conformance and robotics. The default chain also
  runs the whole tree once, as a subprocess started from
  `suites/09_probe_drivers.jl` with the coverage flag passed on, and suite 10
  runs most of `unit/parallel/` again as separate subprocesses.
- `integration/` holds that shared harness, the `examples/` regression entry
  point, and the Python interlink installation checks
  (`interlink_python.jl`, `interlink_python_setup.jl`) that CI's
  `optional-interlink-installation` job runs.
- `gnc/aerobraking/` holds the EDG guidance and control tests that suite 00
  includes.
- `smoke/` holds environment and startup smokes (clean depot, threaded, no-GRAM
  and the examples suite), which `smoke/runtests.jl` lists, and entry points
  that CI workflows run on their own: `installation_contract.jl`,
  `optional_packages_tests.jl`, `contributor_environment_tests.jl`,
  `new_contributor_repository_tests.jl`, `gram_dependency_integration_tests.jl`,
  `gram_package_readiness.jl`, the surrogate preset, workflow and CLI tests,
  `ci_examples_regression.jl` and `ci_readme_quickstart.py`.
- `optional/` holds the companion suites `hypr/`, `edg/` and `visualizer/`,
  each run by its installation job in `julia-ci.yml`.
- `fixtures/` holds tracked fixtures: the Earth regression CSV that suite 04
  reads, the conformance instances, and the small research-derived fixtures in
  `research_contracts/`. The `contributor-kit` CI job tests the starters in
  `templates/starters/`, generated from the commit under test, through the
  reusable contributor workflow.
- `contracts/` chooses which gates run for each tier (`pr_runtests.jl`,
  `nightly_runtests.jl`, `runtests.jl`). `edg_ownership_checks.jl` there is a
  shared check module that some gates include.
- `gates/` holds the `ci_*_gate.jl` implementations and two shared check files.
  `contracts/*.jl` runs most gates; `coverage/runtests.jl` and CI run the
  coverage gate, and CI alone runs the review-record gate.
- `stress/` holds long determinism, Monte Carlo and flake checks, run nightly.
- `coverage/` holds the coverage quality gate and two runtime-analysis gates.
  The probe files whose coverage it measures are in `probes/`.
- `probes/` holds standalone, coverage-targeted probe files:
  - suite 05 includes four of them directly, among them
    `coverage_r6_routing_probes.jl`, whose calibration probes need that suite's
    multi-satellite fixtures;
  - suite 09 starts most of the rest as subprocesses on every run, not only
    under coverage;
  - suite 02 includes `native_probe_reporting_tests.jl`, and starts
    `coverage_threaded_probes.jl` only under `--code-coverage=user`; that file
    in turn includes the GRAM construction and density-service probes;
  - `gram_perturbation_walk_probes.jl` has no driver; its header says how to
    run it.

  Every probe starts the same way as the harness (`using SpaceAGORA` plus
  module aliases), so a probe subprocess uses the precompiled package instead
  of compiling `src/` again. GRAM-backed checks need the
  `SpaceAGORAGRAMSuiteExt` extension, which needs a `data/GRAMSuite.jl`
  checkout with the hooks the extension expects (CI uses the submodule). Local
  runs report missing native prerequisites as skipped. The coverage job with
  native GRAM sets `SPACEAGORA_REQUIRE_NATIVE_GRAM_PROBES=1` and fails unless
  the construction and density-service probes both finish and the child process
  exits successfully.
- `helpers/` holds shared companion loaders, native-probe reporting,
  reference-system sandboxes and reusable scenarios; see
  [the helper notes](helpers/README.md).
- `scripts/` holds Python regression tests of repository scripts, run by CI's
  `plan` job.
- Top-level files: `gmat_scenario_matrix.jl` (the GMAT, Basilisk and CYGNSS
  comparisons, which skip when their data is missing),
  `telemetry_benchmark_manifest.toml` (the verification study's scenarios),
  `p1_findings_allowlist.txt`, and `ai_review_template.md` with the review
  records in `ai_reviews/`.

## Grid adapter suite

`grid_atmosphere/` tests the native-free GRAM grid adapter with ordinary
package imports. Synthetic fixtures cover scalar and batch queries, strict
domain checks, buffered sampling, ENU wind rotation, ownership, serialization,
threaded reads and full-grid schema compatibility. It needs the matching
GRAMSuite pure-grid API and four Julia threads, builds no native model, and
needs no GRAM data, shared library or SPICE kernels. Its dependencies must
already be installed (by `test/grid_atmosphere/setup.jl`); the suite itself
runs no package operations. `EXPECTED_GRAMSUITE_ROOT` can pin the resolved
wrapper folder when you test a separately staged package. CI runs it as its own
step ("Native-free GRAM grid adapter") in the `test-shards` job and again under
coverage in `coverage-shards`; neither the default nor the unit harness
includes it.

An optional regression checks 654 saved Odyssey passage points against an
existing fine-grid interpolation CSV. Set all six variables:
`SPACEAGORA_TEST_GRID_FILE`, `SPACEAGORA_TEST_GRID_FILE_SHA256`,
`SPACEAGORA_TEST_GRID_POINTS`, `SPACEAGORA_TEST_GRID_POINTS_SHA256`,
`SPACEAGORA_TEST_GRID_REFERENCE` and `SPACEAGORA_TEST_GRID_REFERENCE_SHA256`,
which point at the fine `.jls` payload, the passage-point CSV and the
`grid_fine` evaluation CSV and give their hashes. Each input is checked before
and after use. The test suite does not distribute these inputs, and no private
path is built into the tests. This is a regression against kept results, not a
new native comparison or a claim of physical accuracy.
