# Test Layout

Current layout (flat, CI-oriented — a prior purpose-oriented restructure was
reverted to keep this repo mergeable with its upstream):

- `test/runtests.jl`
  Default package test entrypoint. Delegates to `test/integration/runtests.jl`,
  which raw-includes the numbered `test/suites/01..09_*.jl` legacy suites into
  its own shared scope (mock models, builder helpers) — this is still the
  canonical default-entrypoint path, not yet split up.
  The harness loads the package (`using SpaceAGORA`) and binds
  `SimulationModel`, `SimulationEngine`, `TelemetryVerification`, the parallel
  modules and a handful of frame helpers as aliases into it, so there is one
  copy of every module in the process. It does not include `src/`; the one
  exception is `src/mission/operations/maneuver_plans.jl`, which is not part
  of the package. Suites 01 and 03 deliberately raw-include source files into
  throwaway sandbox modules to test standalone loading — keep those.
- `test/unit/`
  A smaller, standalone (`using SpaceAGORA`, no shared-scope dependency) set of
  domain tests, runnable on their own via `test/unit/runtests.jl`. The default
  `test/runtests.jl` chain also runs the whole tree once, as a subprocess
  dispatched from `test/suites/09_probe_drivers.jl` with the coverage flag
  forwarded, so their line data reaches the coverage gate.
- `test/integration/`
  Current home of the legacy end-to-end harness (mock models, builder helpers,
  `test/suites/` includes) plus example, persistence, CLI, and telemetry-facing
  tests.
- `test/smoke/`
  Environment and startup smokes such as clean-depot, threaded, and no-GRAM checks.
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
  Shared native-probe reporting and parent verification live here.

Recommended commands:

```bash
julia --project=. test/runtests.jl
julia --project=. test/smoke/runtests.jl
julia --project=. test/contracts/pr_runtests.jl
julia --project=. test/contracts/nightly_runtests.jl
julia --project=. test/contracts/runtests.jl
julia --project=. test/stress/runtests.jl
julia --project=. test/coverage/runtests.jl
julia --startup-file=no --threads=4 --project=. test/grid_atmosphere/runtests.jl
```

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
