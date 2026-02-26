# SpaceAGORA Verification Contract

## Objective
This contract defines the minimum quality bar for merging changes into SpaceAGORA.  
The target is zero known critical defects (P1), deterministic CI behavior, and sustained nightly stability.

## CI Merge Gates
The following checks are required for pull requests:

1. `tests` in `.github/workflows/julia-ci.yml`
2. `coverage-quality-gate`
3. `p1-findings-gate`
4. `ai-review-artifact-gate`

Nightly checks are required for release readiness and are defined in `.github/workflows/nightly-stress.yml`.

## Branch Protection Rollout
Branch protection on `main` must require these checks:

1. `tests`
2. `coverage-quality-gate`
3. `p1-findings-gate`
4. `ai-review-artifact-gate`

Nightly rollout policy:

1. `nightly-stress` is informational from March 2, 2026 through March 8, 2026.
2. Starting March 9, 2026, release tagging requires a recent successful nightly run.
3. Release-tag enforcement is implemented in `.github/workflows/release-tag-gate.yml`.

## Coverage Policy
Coverage rules are enforced by `test/ci_coverage_quality_gate.jl`.

- Main overall coverage: `>= 90.0%`
- Main per-file coverage: `>= 80.0%`
- Critical-file coverage: `>= 90.0%` for:
  - `src/simulation/run_simulation.jl`
  - `src/control/Propulsive_maneuvers.jl`
  - `src/utils/Closed_form_solution.jl`
  - `src/utils/Save_results.jl`

### Legacy Scope Handling
Legacy code can be handled with explicit policy exceptions only:

- Fully excluded from main coverage threshold:
  - `src/simulation/Complete_passage.jl`
- Legacy per-file override:
  - `src/control/heatload_control/Second_tsw_calcs.jl` with minimum `50.0%`

Excluded/overridden files must stay explicitly documented in this contract and in the gate script.

## P1 Findings Policy
P1 findings are enforced by `test/ci_p1_findings_gate.jl`.

- New unallowlisted P1 markers fail CI.
- Allowlist entries must use `path::exact_line` format.
- Optional metadata may be appended in an inline comment:
  - `owner=<value>`
  - `opened=<YYYY-MM-DD>`
  - `expires=<YYYY-MM-DD>`
- If `expires` is present and in the past, CI fails.

## AI Verification Policy
AI review artifacts are enforced by `test/ci_ai_review_artifact_gate.jl`.

- Required artifact path: `test/ai_reviews/PR_<number>.md`
- Required sections:
  - `Scope`
  - `Changed Files`
  - `Findings`
  - `P1 Assessment`
  - `Tests Added/Updated`
  - `Residual Risk`
- Every changed `src/**/*.jl` file in the PR must be listed in `Changed Files`.

## Nightly Stress Policy
Nightly stress checks run on a schedule in `.github/workflows/nightly-stress.yml` and include:

1. Full test suite with coverage
2. Coverage quality gate
3. P1 findings gate
4. Threaded smoke stress
5. Example regression stress
6. Nightly Monte Carlo stress
7. Flake guard

For pull request runs of the nightly workflow, the AI artifact gate is also required.
