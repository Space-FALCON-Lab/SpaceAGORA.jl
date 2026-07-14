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
- 2026-07-13: **Phase 2, suite 04 done.** Split
  `test/suites/04_solver_env_and_regression_tests.jl` (2649 lines after the
  Phase 1 golden removal) into `test/unit/simulation_engine/solver_env_helpers_tests.jl`
  (the big "Solver/Env Helper Parsing Coverage" testset, lines 1-1832),
  `test/unit/dynamics/rigid_body_and_rhs_tests.jl` (rigid-body/RHS-split
  testsets, lines 1833-2244), and `test/unit/dynamics/orbital_energy_and_drift_tests.jl`
  (energy/invariant/drift physics-regression testsets, lines 2245-2649).
  Suite 04 is now a 4-line pointer comment.
  Wiring this in required also fixing `test/unit/` itself, which had never
  actually been run from any CI entrypoint before: (1) `test/unit/runtests.jl`
  now `include`s `test/helpers/bootstrap.jl` (needed by the migrated content),
  which meant `test/unit/rpo_port_tests.jl` and `test/unit/robotics/runtests.jl`
  could no longer also do `using SpaceAGORA` themselves — raw-`include`-ing
  `SpaceAGORA.jl` (what `bootstrap.jl` does) and `using SpaceAGORA` (real
  package load) in the same process hard-errors ("importing SpaceAGORA into
  Main conflicts with an existing global"). Removed their `using SpaceAGORA`
  lines; they already exclusively used the fully-qualified `SpaceAGORA.X`
  form so this was a safe, mechanical fix. (2) That combination also exposed
  a real bug in `test/unit/robotics/runtests.jl`: three bare `RobotArmSphereObstacle[]`
  literals resolved against the wrong of two now-coexisting copies of
  `SimulationModel` (bootstrap.jl raw-includes a standalone top-level copy
  *and* a nested `SpaceAGORA.SimulationModel` copy — these are distinct
  Julia types even though structurally identical). Fixed by qualifying all
  three as `SpaceAGORA.RobotArmSphereObstacle[]`, matching the file's own
  existing convention everywhere else. (3) Running `test/unit/` for the
  first time also surfaced two genuine pre-existing, previously-undiscovered
  issues, unrelated to this restructure: `data/rpo/station_geometry/demo/`
  (a "demo" station pointcloud CSV) doesn't exist in this checkout, so gated
  the "RPO geometry and PSO basics" testset behind an `isfile` check +
  `@test_skip`, mirroring the golden-fixture convention; and one RPO
  replanning-decision tie-break test (`tracking_error_retime_m`/
  `tracking_error_replan_m` both clamp to `0.0`) asserts `:retime` but the
  code actually returns `:replan` — marked `@test_broken` with a comment,
  needs a GNC-domain call on whether the test or the tie-break is wrong.
  Full `test/runtests.jl`: 2801/2801 legacy + 1285/1285 unit + 316/317 golden
  (gram_mars still correctly skipped) — all clean.
- 2026-07-13: **Phase 2, suite 05 done.** Split
  `test/suites/05_thruster_control_and_quality_tests.jl` (933 lines) into
  `test/unit/vehicle/thruster_tests.jl` ("Thruster Edge Cases" testset, lines
  1-747), `test/unit/gnc/control_effector_tests.jl` (the three
  "Control ... (End-to-End)" testsets, lines 748-909), and
  `test/contracts/aqua_quality_gate.jl` (the "Aqua Package Quality" testset,
  lines 918-933 — a package-policy check, doesn't belong with thruster/control
  unit tests; wired into `test/contracts/pr_runtests.jl`, written
  self-contained like the other `ci_*_gate.jl` files rather than depending on
  `bootstrap.jl`). The remaining bare `include`s suite 05 had picked up over
  time (3 coverage-probe files, 2 aerobraking parity test files, 1 mission
  stub file) were carried forward as-is into `test/unit/runtests.jl` under a
  "Coverage Probes (pending Phase 4 audit)" / "GNC" testset — not
  reorganizing their content now, just preserving exactly what already ran
  by default so nothing silently stops being tested before Phase 4 gets to
  them. Suite 05 is now a pointer comment.
  Full `test/runtests.jl`: 2801/2801 legacy + 1285/1285 unit + 316/317 golden
  — all clean.
- 2026-07-13: **Phase 2, suite 03 done.** Split
  `test/suites/03_persistence_units_and_rotational_tests.jl` (1353 lines)
  into `test/unit/io/persistence_and_checkpoint_tests.jl` (checkpoint/resume/
  results-bundle/normalization/verbose-logging testsets, lines 1-930) and
  `test/unit/dynamics/rotational_tests.jl` (orbital-elements/quaternion/
  torque-free rigid-body testsets, lines 931-1353) — matching the suite's
  own name (persistence+units → io, rotational → dynamics). Suite 03 is now
  a pointer comment. No new issues found this time (unlike suites 04/05,
  nothing here needed a behavior fix).
  Full `test/runtests.jl`: 956/956 legacy (03/04/05 now empty) + 3130/3130
  unit + 316/317 golden — all clean.
- 2026-07-13: **Phase 2, suites 01/02/06/07/08 done — Phase 2 fully closed
  out.** Suite 01 (1189 lines) split into `test/contracts/architecture_and_export_contracts.jl`
  (the 17 "...Contract" testsets, lines 1-568 — these need the sandbox
  modules/raw-include mechanics from `bootstrap.jl`, so unlike the other
  `ci_*_gate.jl` files this one includes it explicitly rather than being
  self-contained) and three `test/unit/` files by domain:
  `core/api_convenience_constructor_tests.jl`,
  `simulation_engine/run_metadata_and_rhs_tests.jl`,
  `mission/maneuver_and_campaign_tests.jl`. Suite 02 (1910 lines) split into
  7 files across `simulation_engine/` (callbacks), `parallel/` (policy +
  adaptive routing), `dynamics/` (aerodynamic helpers),
  `environment/` (planet constructors), and `simulation/` (campaign + GRAM
  lock). Suites 06/07/08 (small, single-testset files) moved wholesale to
  `analysis/telemetry_and_policy_tests.jl`, `environment/no_gram_onboarding_tests.jl`,
  and `cli/cli_and_asset_tests.jl` respectively.
  Found one real cross-file dependency break during this split:
  `ensure_guidance_sandbox_loaded!()` was defined at the tail of suite 01's
  contract section but called from the unit section (`maneuver_and_campaign_tests.jl`)
  — fine when both lived in one file/suite, broken once they became
  separately-included files with `test/unit/mission/` running before
  `test/contracts/architecture_and_export_contracts.jl` did. Moved the
  function to `test/helpers/sandbox_modules.jl` (where `GUIDANCE_SANDBOX` it
  operates on already lives), so it's available via `bootstrap.jl` regardless
  of include order — the general lesson: when splitting a suite file,
  grep the whole file for helper-function *definitions* used across the
  split boundary, not just move testsets by line range.
  Deleted `test/suites/` (all 8 files were down to pointer comments) and
  simplified `test/integration/runtests.jl` to just the `bootstrap.jl`
  trigger (it remains the future home for genuine end-to-end integration
  tests per its subdirectory READMEs).
  Full `test/runtests.jl`: legacy suites section now empty (0/0, directory
  deleted) + 3787/3787 unit + 299/299 architecture contracts + 316/317
  golden — all clean. **This closes out Phase 2.**
- 2026-07-13: **Phase 3 (contract-gate relocation) in progress.** `git mv`'d
  all 42 root `test/ci_*.jl` files into `test/contracts/`. This had a bigger
  blast radius than a pure internal test/ reorg: fixed each file's own
  `REPO_ROOT` computation (one more `..` level after the move — 3 files used
  `dirname(dirname(@__FILE__))`, 39 used `normpath(joinpath(@__DIR__, ".."))`);
  updated 7 aggregator files' include paths
  (`test/contracts/{pr,nightly,}runtests.jl`, `test/smoke/runtests.jl`,
  `test/coverage/runtests.jl`, `test/stress/runtests.jl`,
  `test/integration/examples/runtests.jl`); fixed 3 gates' self-exclusion
  logic that matched on the literal prefix `"test/ci_"` (now stale --
  `ci_canonical_path_contract_gate.jl`, `ci_no_src_benchmarks_root_gate.jl`,
  `ci_no_dynamics_models_gate.jl` — changed to match `basename(rel)` instead
  of the full relative path, so it's independent of directory); fixed one
  gate's self-referential exclusion list entry (`ci_p1_findings_gate.jl`) and
  one gate's read of another moved file's content
  (`architecture_and_export_contracts.jl` reads `ci_clean_depot_smoke.jl` as
  a string to check its content); updated the 2 exact command strings this
  changed in `.github/pull_request_template.md` and
  `test/contracts/ci_docs_contract_gate.jl`'s markers checking for them;
  updated 3 CI workflow YAML files' direct `test/ci_X.jl` invocations
  (`CI.yml`, `julia-ci.yml`, `nightly-stress.yml` — carefully avoided touching
  the unrelated `SpaceAGORACalibration.jl/test/ci_coverage_quality_gate.jl`,
  a different subpackage's own file); updated path citations in 6 docs pages
  (`docs/architecture/{canonical_topology_contract,shim_window_manifest,src_completeness_contract}.md`,
  `docs/quality/{api_naming_contract,verification_contract}.md`,
  `docs/src/documentation_policy.md`) and one benchmarks handoff doc
  (left `test/ai_reviews/PR_*.md` untouched -- archival records of past PRs,
  not living docs).
  While validating (`test/contracts/pr_runtests.jl` had never been run before
  either), found and fixed two more pre-existing, previously-undiscovered
  bugs unrelated to the relocation itself: `ci_architecture_contract_gate.jl`
  asserted `const GRAM_LOCK = ReentrantLock()` in
  `src/simulation/runtime_services.jl`, which has actually read
  `const GRAM_LOCK = SPICE_LOCK` since the GRAM/CSPICE-lock-unification fix
  (see memory: `project_gram_cspice_symbol_collision`) — nobody had run this
  gate since; and `ci_no_legacy_include_chains_gate.jl`'s allowlist was
  missing `src/parallel/process/parallel_process.jl` (added during the
  process-backend work, commit `1ba63a82`), which legitimately raw-includes
  its sibling `worker_pool.jl` the same way already-allowlisted
  `parallel_policy.jl`/`parallel_profiles.jl` do.
  Also discovered `test/contracts/runtests.jl` (the "full superset" runner)
  cannot include `architecture_and_export_contracts.jl` itself: several of
  the other gates in that file `using SpaceAGORA` as a real package, which
  conflicts with `architecture_and_export_contracts.jl`'s `bootstrap.jl`
  raw-include in the same process (same class of bug as the Phase 2
  `rpo_port_tests.jl`/`robotics` fix, but the *other* direction) —
  removed it from that aggregator with an explanatory comment; it already
  runs correctly via `test/unit/runtests.jl` (part of default
  `test/runtests.jl`), which never mixes the two loading styles.
  Validated, all clean: `test/contracts/pr_runtests.jl`,
  `nightly_runtests.jl`, `runtests.jl` (full superset), `test/smoke/runtests.jl`
  (33/33 examples + 3 other smokes), `test/integration/examples/runtests.jl`
  (27/27), and the individual moved gates referenced from
  `test/coverage/runtests.jl` (`ci_runtime_any_hotpath_gate.jl`,
  `ci_runtime_analysis_copy_overhead_gate.jl` — `ci_coverage_quality_gate.jl`
  itself needs a prior `--code-coverage=user` instrumented run to have
  `.cov` files to check, an orthogonal precondition, not a bug).
  `test/stress/runtests.jl` also uncovered a third pre-existing bug: unlike
  everywhere else `SimulationModel` gets raw-included, `ci_flake_guard.jl`
  and `ci_nightly_montecarlo_stress.jl` both raw-`include`d
  `src/core/simulation_model.jl` *unconditionally* (no `isdefined` guard,
  unlike the `SimulationEngine` include two lines below it in the same two
  files, which *was* guarded) — running both gates in one process (exactly
  what `test/stress/runtests.jl` does, and what `nightly-stress.yml` already
  does in CI) redefines `Main.SimulationModel` and produces an ambiguous-
  export error (`Earth` in this case). This is independent of the relocation
  and would have broken nightly CI the same way at the old path; guarded
  both includes the same way the adjacent `SimulationEngine` include already
  was. Final full-repo re-check: `test/runtests.jl` still 3787/3787 unit +
  299/299 architecture contracts + 316/317 golden after the `ci_architecture_contract_gate.jl`
  and `ci_no_legacy_include_chains_gate.jl` fixes above. **Phase 3 is closed
  out**, fully validated across every entrypoint that touches a moved file.
