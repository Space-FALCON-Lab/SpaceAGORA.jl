# Test Suite Restructure Plan

Status doc for an in-progress, multi-session effort to restructure `test/` into
two clear, extensible pillars: **function-level unit tests** and
**reproducible propagation (golden) tests**. This file is the source of truth
for what's done and what's next — update it as phases complete instead of
relying on conversation history.

## Why

The current suite (~90 files, ~21.7k lines) is a single monolith:
`test/runtests.jl` -> `test/integration/runtests.jl` (mocks, builders, and a
~90-function alias block exposing `SimulationEngine` internals) -> 8 numbered
`test/suites/NN_*.jl` files (~18.6k lines) that freely mix unit-level physics
checks, regression/golden checks, contract checks, and CLI/smoke checks with
no directory boundary. A second, disconnected axis of ~40 root-level
`ci_*_gate.jl` policy/architecture checks was never relocated into
`test/contracts/` despite that folder existing for exactly that purpose.
Coverage is chased via dedicated `coverage_*_probes.jl` files written to hit a
numeric threshold rather than to exercise real behavior. Reproducible
propagation testing exists as exactly one instance (the AGORA Earth golden
regression, buried at `test/suites/04_solver_env_and_regression_tests.jl:2520`)
with no pattern for adding a second scenario.

`test/README.md` already describes a target layout (`unit/`, `integration/`,
`contracts/`, `smoke/`, `stress/`, `coverage/`) but the migration stalled:
`test/unit/` is a placeholder and `test/integration/{cli,persistence,simulation,telemetry}/`
are README-only stubs.

## Decisions made

- **Scope**: full cleanup, including relocating `ci_*_gate.jl` into
  `test/contracts/` and auditing the coverage-probe files for redundancy.
- **Cadence**: incremental, domain by domain, each as its own reviewable PR —
  not a big-bang rewrite.
- **Golden scenario coverage**: formalize the existing AGORA Earth case, plus
  add a GRAM/Mars-backed scenario and a multi-satellite constellation
  scenario, to prove the harness generalizes before calling it done.
- **Internals access in new tests**: new `test/unit/` files reference module
  functions explicitly per file (e.g. `SimulationEngine.foo`) instead of
  relying on the shared ~90-name alias block. That block stays in
  `test/integration/runtests.jl` for the *not-yet-migrated* legacy suites
  until they're migrated away from it in Phase 2 — don't remove it early or
  suites 01-08 break.

## Phases

### Phase 0 — Scaffolding (no test migration yet)
- Extract mock models/builder functions out of `test/integration/runtests.jl`
  into `test/helpers/*.jl`, included in the same order so behavior is
  unchanged (Julia `include()` at top level means file location doesn't
  affect evaluation order/semantics, only inclusion order does).
- Stand up `test/unit/<domain>/` directories mirroring `src/` for the domains
  Phase 2 will target first.
- Document the internals-access convention (see above) in `test/README.md`.
- Stand up `test/golden/<scenario>/` layout convention + shared harness.

### Phase 1 — Golden/propagation pattern, proven with 3 scenarios
- `test/golden/agora_earth/`: migrate the existing regression fixture + test,
  add provenance metadata (git SHA, Julia version, generation date).
- `test/golden/gram_mars/` (or similar): new GRAM/Mars-backed scenario.
- `test/golden/constellation/` (or similar): new multi-satellite scenario.
- Shared `run_and_compare_golden(scenario)` harness in `test/helpers/` so
  adding a 4th scenario is: write a config, run the regen script once, add
  one entry — not hand-roll a testset.
- `scripts/regenerate_golden.jl` tool + documented regen policy.
- Fast (non-GRAM) golden cases run in the PR tier; GRAM/SPICE-heavy ones in
  the nightly tier, mirroring the existing contracts pr/nightly split.

### Phase 2 (future) — Incremental unit-test migration, domain by domain
Order (largest/most tangled first): suite 04 (solver/env/regression, 2696
lines) -> suite 05 (thruster/control, 933 lines) -> suite 03
(persistence/rotational, 1353 lines) -> suites 01, 02, 06, 07, 08. Each PR:
move testsets into `test/unit/<domain>/`, verify green, delete the migrated
block from the suite file. `test/suites/` gets deleted once empty.

### Phase 3 (future) — Contract-gate relocation
Pure move of the ~40 root `ci_*_gate.jl` files into `test/contracts/`,
updating include paths in `test/contracts/{pr,nightly,}runtests.jl`. No
behavior change.

### Phase 4 (future) — Coverage-probe audit
For each of `coverage_parallel_telemetry_probes.jl`,
`coverage_runtime_boundary_probes.jl`, `coverage_targeted_90_probes.jl`: check
whether Phase 2's unit tests now naturally cover the same lines; retire
redundant probes, relocate genuinely-useful ones into `test/unit/` next to the
code they cover.

### Phase 5 (future) — Wire-up and docs
Update `test/runtests.jl` and CI workflows for a fast (unit + fast golden)
vs. slow (propagation-heavy + contracts + stress) split. Rewrite
`test/README.md` as the definitive "how to add a test" doc.

## Progress log

- 2026-07-13: Plan written. Starting Phase 0 and Phase 1.
- 2026-07-13: **Phase 0 done.** Extracted mock models/builders/numeric
  helpers/sandbox modules out of `test/integration/runtests.jl` into 8 files
  under `test/helpers/`; split the remaining setup into a shared,
  double-include-guarded `test/helpers/bootstrap.jl` used by both
  `test/integration/runtests.jl` (now a thin suite-runner) and
  `test/golden/runtests.jl`. Stood up `test/unit/<domain>/` placeholders and
  documented the new convention in `test/README.md` / `test/helpers/README.md`.
  Verified with a full `test/runtests.jl` run: 4088 pass / 1 broken (the one
  broken test is pre-existing and unrelated) across all 8 suites — the
  extraction changed no behavior.
- 2026-07-13: **Phase 1 done.** Built `test/helpers/golden_harness.jl`
  (`run_and_compare_golden`, generalized to N spacecraft; PR/nightly tier
  gating via `metadata.toml` + `SPACEAGORA_GOLDEN_TIER`) and
  `scripts/regenerate_golden.jl`. Three scenarios now live under
  `test/golden/`: `agora_earth` (migrated from the old buried suite-04 test),
  `constellation` (new, 3 satellites via `build_config_multi`), and
  `gram_mars` (new, GRAM-backed Mars drag case). All three pass end-to-end
  (372/372 assertions with `SPACEAGORA_GOLDEN_TIER=all`).
  Along the way, found and fixed a real gap in `test/helpers/bootstrap.jl`'s
  GRAM re-patch block: it mirrored most of
  `ext/SpaceAGORAGRAMSuiteExt.jl.__init__` but was missing the Mars/non-Earth
  ephemeris-state bypass (`GRAMSuite._GRAM_EPHEMERIS_STATE_FN[]` /
  `_GRAM_DEFAULT_LOCK_HOOK[]`) — because the test harness raw-`include`s
  `SpaceAGORA.jl` rather than `using` it as a package, Julia's extension
  auto-load never triggers, so this hook has to be reinstalled by hand just
  like the rest of that block already was. No prior test exercised a live
  Mars+GRAM query, so this had never surfaced before; see memory
  `project_gram_mars_isolated_cspice_bug`. Also removed the old golden
  testset from `test/suites/04_solver_env_and_regression_tests.jl` and the
  superseded flat `test/golden/agora_earth_regression.csv`.
  Wired `test/golden/runtests.jl` into the default `test/runtests.jl`.
  Confirmed with a final default `julia --project=. test/runtests.jl` run:
  legacy suites 4088/4088 pass (the earlier "1 broken" was a pre-existing,
  unrelated flake — it doesn't appear in this run; not something introduced
  here, see `ci_flake_guard.jl`), golden 316/316 pass with `gram_mars`
  correctly `@test_skip`'d by the default PR-tier gate. **Phase 1 is closed
  out.**
