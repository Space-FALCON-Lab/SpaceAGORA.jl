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
