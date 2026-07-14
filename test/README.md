# Test Layout

This directory is being migrated from a flat, CI-oriented layout to a purpose-oriented layout.

Current intent:

- `test/runtests.jl`
  Default package test entrypoint. For now it delegates to `test/integration/runtests.jl` so existing behavior stays unchanged.
- `test/unit/`
  Future home for fast, isolated logic tests after the legacy suites are split up.
- `test/integration/`
  Current home of the legacy end-to-end harness plus example, persistence, CLI, and telemetry-facing tests.
- `test/smoke/`
  Environment and startup smokes such as clean-depot, threaded, and no-GRAM checks.
- `test/contracts/`
  Architecture, API-surface, boundary, naming, docs, and repository policy gates.
- `test/stress/`
  Long-running determinism, Monte Carlo, and flake-resistance checks.
- `test/coverage/`
  Coverage quality gates and targeted probe suites.
- `test/helpers/`
  Shared harness code and fixtures once the large legacy test file is decomposed.

Recommended commands:

```bash
julia --project=. test/runtests.jl
julia --project=. test/smoke/runtests.jl
julia --project=. test/contracts/pr_runtests.jl
julia --project=. test/contracts/nightly_runtests.jl
julia --project=. test/contracts/runtests.jl
julia --project=. test/stress/runtests.jl
julia --project=. test/coverage/runtests.jl
```

Migration notes:

- The current legacy runtime harness now lives at `test/integration/runtests.jl`.
- Most standalone `ci_*.jl` files have not been renamed yet; the new aggregators include them from their current paths.
- `test/unit/runtests.jl` is intentionally light until the large numbered suites are split into smaller domain files.
- See `test/TEST_RESTRUCTURE_PLAN.md` for the active, multi-session plan and progress log for this migration.

## Adding a new test

**Function-level unit test** (fast, no SPICE/GRAM data, calls a function
directly rather than running a full `run_simulation`):

- Add it under `test/unit/<domain>/`, mirroring the `src/<domain>/` path of
  the code under test (e.g. a test for `src/vehicle/actuators/foo.jl` goes in
  `test/unit/vehicle/actuators/foo_tests.jl`).
- Reference the module function explicitly in your test file, e.g.
  `SimulationEngine.some_internal_fn(...)` — don't add new names to the
  shared alias block in `test/integration/runtests.jl`. That block only
  exists for the not-yet-migrated `test/suites/NN_*.jl` legacy files; it is
  being phased out, not extended.
- If you need a mock model or builder function, check `test/helpers/` first
  (`mock_models.jl`, `spacecraft_builders.jl`, `config_builders.jl`,
  `simulation_run_helpers.jl`, `numeric_helpers.jl`, `test_constants.jl`,
  `persistence_fixtures.jl`) before writing a new one.

**Reproducible propagation (golden) test** — a full simulation run compared
against a checked-in fixture within tolerance:

- Add a new scenario under `test/golden/<scenario_name>/` following the
  pattern established there (config, fixture CSV, provenance metadata) and
  use the shared harness rather than hand-writing a bespoke comparison
  testset. See `test/golden/README.md` once Phase 1 lands, or
  `test/TEST_RESTRUCTURE_PLAN.md` in the meantime.

**Everything else** (architecture/naming/boundary policy gates, coverage
probes, stress/Monte Carlo, environment smokes) has its own existing bucket
above — add alongside the existing files of that kind rather than in
`test/unit/` or `test/golden/`.
