# Test Layout

Two extensible pillars, plus supporting policy/quality infrastructure:

- **Function-level unit tests** — fast, call functions/APIs directly, no full
  `run_simulation`. Live under `test/unit/`, mirroring `src/`.
- **Reproducible propagation ("golden") tests** — a full `run_simulation`
  compared against a checked-in fixture within tolerance. Live under
  `test/golden/`.

Everything else (architecture/policy gates, coverage quality, stress/Monte
Carlo, environment smokes) has its own bucket below.

- `test/runtests.jl`
  Default package test entrypoint (what CI runs). Runs `test/integration/runtests.jl` (bootstrap trigger; legacy suites are gone, see below), `test/unit/runtests.jl`, then `test/golden/runtests.jl` (PR tier only by default).
- `test/unit/`
  Function-level unit tests, organized by domain (`dynamics/`, `simulation_engine/`, `vehicle/`, `io/`, `environment/`, `gnc/`, `parallel/`, `mission/`, `core/`, `analysis/`, `cli/`, `coverage_probes/`, plus `rpo_port_tests.jl` and `robotics/`). See "Adding a new test" below.
- `test/golden/`
  Reproducible-propagation regression tests. See `test/golden/README.md`.
- `test/integration/`
  Bootstrap trigger point (`runtests.jl` just includes `test/helpers/bootstrap.jl`) and the future home for genuine end-to-end integration tests (examples, CLI, persistence, telemetry — see the subdirectory READMEs). The old numbered-suite legacy harness that used to live here was fully migrated out and deleted (see `test/TEST_RESTRUCTURE_PLAN.md`).
- `test/smoke/`
  Environment and startup smokes: clean-depot, threaded, no-GRAM, examples-suite.
- `test/contracts/`
  Architecture, API-surface, boundary, naming, docs, and repository policy gates (the `ci_*_gate.jl` files), split into `pr_runtests.jl` (fast, every PR), `nightly_runtests.jl` (slower, nightly), and `runtests.jl` (full superset for local runs).
- `test/stress/`
  Long-running determinism, Monte Carlo, and flake-resistance checks.
- `test/coverage/`
  Coverage quality gate plus the coverage-probe files' non-domain-specific gates. (The probe test content itself lives at `test/unit/coverage_probes/`.)
- `test/helpers/`
  Shared harness code: `bootstrap.jl` (module/engine setup, GRAM extension re-patching — included by `test/integration/`, `test/unit/`, and `test/golden/` runners) plus mock models, builders, and the golden-test harness. See `test/helpers/README.md`.

Recommended commands:

```bash
julia --project=. test/runtests.jl                                          # default: unit + golden (PR tier)
julia --project=. test/unit/runtests.jl
julia --project=. test/golden/runtests.jl
SPACEAGORA_GOLDEN_TIER=all julia --project=. test/golden/runtests.jl        # include nightly-tier golden scenarios
julia --project=. test/smoke/runtests.jl
julia --project=. test/contracts/pr_runtests.jl
julia --project=. test/contracts/nightly_runtests.jl
julia --project=. test/contracts/runtests.jl                                # full superset, local use
julia --project=. test/stress/runtests.jl
julia --project=. test/coverage/runtests.jl                                 # needs --code-coverage=user data first
```

See `test/TEST_RESTRUCTURE_PLAN.md` for the history of how this layout came
to be and a log of what changed in each phase.

## Adding a new test

**Function-level unit test:**

- Add it under `test/unit/<domain>/`, mirroring the `src/<domain>/` path of
  the code under test (e.g. a test for `src/vehicle/actuators/foo.jl` goes in
  `test/unit/vehicle/actuators/foo_tests.jl`).
- Reference module functions explicitly, e.g. `SimulationEngine.some_internal_fn(...)`.
  Don't rely on unqualified names beyond what's already exported/aliased by
  `test/helpers/bootstrap.jl` — new unit test files should be explicit about
  what they depend on rather than accumulating more implicit global aliases.
- If you need a mock model or builder function, check `test/helpers/` first
  (`mock_models.jl`, `spacecraft_builders.jl`, `config_builders.jl`,
  `simulation_run_helpers.jl`, `numeric_helpers.jl`, `test_constants.jl`,
  `persistence_fixtures.jl`) before writing a new one.
- Wire your new file into the matching `@testset` block in `test/unit/runtests.jl`.

**Reproducible propagation (golden) test:**

- Add a new scenario under `test/golden/<scenario_name>/` following the
  pattern established there (config, fixture CSV, provenance metadata) and
  use the shared harness (`run_and_compare_golden`) rather than hand-writing
  a bespoke comparison testset. See `test/golden/README.md`.

**Architecture/policy gate:**

- Add a new `test/contracts/ci_<name>_gate.jl` file, self-contained (own
  `using`s, own `REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))`) like
  its siblings — don't depend on `test/helpers/bootstrap.jl` unless you have
  a specific reason to (see the comment atop `architecture_and_export_contracts.jl`
  for why that one is the exception, and why mixing `using SpaceAGORA` gates
  with `bootstrap.jl`-based ones in the same process breaks).
- Wire it into `test/contracts/pr_runtests.jl` (every PR) or
  `test/contracts/nightly_runtests.jl` (slower/nightly-only) as appropriate.

**Everything else** (coverage probes, stress/Monte Carlo, environment
smokes) has its own existing bucket above — add alongside the existing files
of that kind.
