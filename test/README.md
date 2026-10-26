# Test Layout

Test folders follow their execution and dependency boundaries. Shared drivers
select the suites; helpers supply reusable setup and scenarios.

- `test/runtests.jl`
  Default package test entrypoint. Delegates to `test/integration/runtests.jl`,
  which raw-includes the numbered `test/suites/` legacy suites into
  its own shared scope (mock models, builder helpers) — this is still the
  canonical default-entrypoint path, not yet split up.
  The harness loads the package (`using SpaceAGORA`) and binds
  `SimulationModel`, `SimulationEngine`, `TelemetryVerification`, the parallel
  modules and a handful of frame helpers as aliases into it, so there is one
  copy of every module in the process. It does not include `src/`; the one
  exception is `src/orchestration/mission/operations/maneuver_plans.jl`, which is not part
  of the package. Suites 01 and 03 deliberately raw-include source files into
  throwaway sandbox modules to test standalone loading — keep those.
  Reference-system sandboxes use `test/helpers/reference_system_sandbox.jl`:
  it raw-loads the canonical Geodesy and OrbitalGeometry numeric sources before
  the frame leaf, without importing SpaceAGORA or RuntimeServices. This preserves
  includer-local `args` and the cold missing-lock probe while keeping dependency
  orchestration out of the internal `reference_system.jl` leaf.
- `test/unit/`
  A smaller, standalone (`using SpaceAGORA`, no shared-scope dependency) set of
  domain tests, runnable on their own via
  `julia --project=test test/unit/runtests.jl`. The default
  `test/runtests.jl` chain also runs the whole tree once, as a subprocess
  dispatched from `test/suites/09_probe_drivers.jl` with the coverage flag
  forwarded, so their line data reaches the coverage gate.
- `test/integration/`
  The shared integration harness (mock models, builder helpers and `test/suites/`
  includes), plus the `examples/` regression entrypoint. Domain folders are added
  when they contain tests, rather than reserved in advance.
- `test/smoke/`
  Environment and startup smokes such as clean-depot, threaded, and no-GRAM checks,
  plus the workflow entry points that check installation and surrogate contracts
  (`installation_contract.jl`, `optional_packages_tests.jl`, `contributor_environment_tests.jl`,
  `gram_dependency_integration_tests.jl`, the surrogate preset, workflow and CLI
  tests). `runtests.jl` lists its files explicitly; the entry points are run by
  the workflows that name them.
- `test/optional/`
  The companion suites: `hypr/`, `edg/` and `visualizer/`, each driven by its
  installation job in `julia-ci.yml`.
- `test/fixtures/`
  Tracked fixtures such as the Earth regression CSV read by suite 04. (The
  `contributor-kit` job tests the starters of `templates/starters/`, generated
  from the commit under test, through the reusable contributor workflow.)
- `test/grid_atmosphere/`
  A standalone native-free grid-adapter suite with normal package imports.
  Synthetic fixtures cover scalar/batch queries, strict domain checks, buffered
  sampling, ENU wind rotation, ownership, serialization, threaded reads and
  full-grid schema compatibility. It requires the matching GRAMSuite pure-grid
  API and four Julia threads. The `tests-matrix` CI job runs it once as a separate
  bounded step; it is not also included in the default or unit harness.
- `test/contracts/`
  Orchestration only (`pr_runtests.jl`, `nightly_runtests.jl`, `runtests.jl`);
  the architecture/API-surface/boundary/naming/docs/policy gate implementations
  themselves live in `test/gates/`.
- `test/gates/`
  The `ci_*_gate.jl` implementations wired in by `test/contracts/*.jl`.
- `test/stress/`
  Long-running determinism, Monte Carlo, and flake-resistance checks.
- `test/coverage/`
  Coverage quality gate plus two runtime-analysis gates. The probe suites it
  measures coverage of live in `test/probes/`, not here.
- `test/probes/`
  Standalone (`using SpaceAGORA`) coverage-targeted probe files. Four are
  raw-included directly from `test/suites/05_thruster_control_and_quality_tests.jl`
  (`coverage_r6_routing_probes.jl` among them, because its calibration probes
  need that suite's multi-satellite fixtures);
  most of the rest are dispatched as subprocesses from
  `test/suites/09_probe_drivers.jl`, every run (not just under coverage).
  `coverage_threaded_probes.jl` is the one exception — its
  `test/suites/02_callbacks_parallel_and_smoke_tests.jl` driver only
  dispatches it when running with `--code-coverage=user`.
  Every probe bootstraps the same way as the harness (`using SpaceAGORA` plus
  module aliases), so a probe subprocess starts from the precompiled package
  instead of recompiling `src/`. GRAM-backed probe checks need the
  `SpaceAGORAGRAMSuiteExt` extension to load, which requires a vendored
  `data/GRAMSuite.jl` checkout that provides the hooks the extension expects
  (CI uses the dev submodule). Optional local runs report unavailable native
  prerequisites as skipped. The native-enabled coverage job sets
  `SPACEAGORA_REQUIRE_NATIVE_GRAM_PROBES=1` and fails unless construction and
  density-service probes both finish and the child exits successfully.
- `test/helpers/`
  Shared companion loaders, native-probe reporting, reference-system sandboxes
  and reusable scenarios live here. See [the helper guide](helpers/README.md).

Recommended commands:

```bash
julia --startup-file=no test/setup.jl
julia --startup-file=no examples/setup.jl
julia --project=test test/runtests.jl
julia --project=. test/smoke/runtests.jl
julia --project=test test/contracts/pr_runtests.jl
julia --project=test test/contracts/nightly_runtests.jl
julia --project=test test/contracts/runtests.jl
julia --project=. test/stress/runtests.jl
julia --project=. test/coverage/runtests.jl
julia --startup-file=no --threads=4 --project=. test/grid_atmosphere/runtests.jl
```

The test environment (`test/Project.toml`) holds the test tools, Aqua among
them, and every package the harness imports, with SpaceAGORA from this checkout;
the root project does not include Aqua. Suites 01 and 03 load the package source
into sandbox modules, so the test environment declares all of SpaceAGORA's own
dependencies too, and `test/setup.jl` refuses one that misses any. The tests
also start example and CLI children, which run in the examples environment
(`examples/Project.toml`, prepared by `examples/setup.jl`), as do the example
smoke runs of `test/smoke/runtests.jl`.

`test/setup.jl` generates the ignored `test/Manifest.toml` from the root
manifest on every run, so the tests use the package versions the root project
pins. Run it again after the root manifest changes. `Pkg.test()` from the root
project uses the same test project. The harness requires Aqua and stops when it
would test another checkout's SpaceAGORA. The EDG and HYPR providers stay
separate environments: `SPACEAGORA_EDG_TEST_ENV` and `SPACEAGORA_HYPR_TEST_ENV`
name prepared ones (`scripts/setup_edg.jl`, `scripts/setup_hypr.jl`).

The grid suite constructs no native model and needs no GRAM data, shared library
or SPICE kernels. `EXPECTED_GRAMSUITE_ROOT` can pin the resolved wrapper directory
when testing a separately staged package. Package dependencies must already be
installed; the suite performs no package operations.

The optional retained Odyssey fine-grid regression checks 654 saved passage
points against an existing interpolation CSV. Supply all six variables:
`SPACEAGORA_TEST_GRID_FILE`, `SPACEAGORA_TEST_GRID_FILE_SHA256`,
`SPACEAGORA_TEST_GRID_POINTS`, `SPACEAGORA_TEST_GRID_POINTS_SHA256`,
`SPACEAGORA_TEST_GRID_REFERENCE`, and `SPACEAGORA_TEST_GRID_REFERENCE_SHA256`.
Paths refer respectively to the retained fine `.jls` payload, passage point CSV,
and `grid_fine` evaluation CSV. Each input is checked before and after use.
These inputs are optional and are not distributed by this test suite. This is a
software regression against retained results, not a fresh native comparison or
a claim of physical accuracy. No private workspace path is built into the tests.

Migration notes:

- The current legacy runtime harness now lives at `test/integration/runtests.jl`.
- `ci_*.jl` gate implementations live in `test/gates/`; `test/contracts/*.jl`
  only orchestrates which subset runs for a given tier.
- `test/unit/runtests.jl` covers a growing but still partial slice of domains
  (RPO, robotics) — the bulk of coverage remains in the numbered
  `test/suites/` legacy suites, not yet split into `test/unit/`.
