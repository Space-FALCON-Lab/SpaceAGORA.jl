# Test Layout

This directory is being migrated from a flat, CI-oriented layout to a purpose-oriented layout.

Current intent:

- `test/runtests.jl`
  Default package test entrypoint. Delegates to `test/integration/runtests.jl`,
  which raw-includes the numbered `test/suites/01..09_*.jl` legacy suites into
  its own shared scope (mock models, builder helpers) — this is still the
  canonical default-entrypoint path, not yet split up.
- `test/unit/`
  A smaller, standalone (`using SpaceAGORA`, no shared-scope dependency) set of
  domain tests, separate from and not included by the default `test/runtests.jl`
  chain — run via `test/unit/runtests.jl`.
- `test/integration/`
  Current home of the legacy end-to-end harness (mock models, builder helpers,
  `test/suites/` includes) plus example, persistence, CLI, and telemetry-facing
  tests.
- `test/smoke/`
  Environment and startup smokes such as clean-depot, threaded, and no-GRAM checks.
- `test/contracts/`
  Orchestration only (`pr_runtests.jl`, `nightly_runtests.jl`, `runtests.jl`);
  the architecture/API-surface/boundary/naming/docs/policy gate implementations
  themselves live in `test/gates/`.
- `test/gates/`
  The `ci_*_gate.jl` implementations wired in by `test/contracts/*.jl`.
- `test/stress/`
  Long-running determinism, Monte Carlo, and flake-resistance checks.
- `test/coverage/`
  Coverage quality gates and targeted probe suites.
- `test/helpers/`
  Currently unused (placeholder); no shared harness code lives here.

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
- `ci_*.jl` gate implementations live in `test/gates/`; `test/contracts/*.jl`
  only orchestrates which subset runs for a given tier.
- `test/unit/runtests.jl` covers a growing but still partial slice of domains
  (RPO, robotics) — the bulk of coverage remains in the numbered
  `test/suites/` legacy suites, not yet split into `test/unit/`.
