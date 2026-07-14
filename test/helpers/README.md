Shared fixtures and harness code, extracted out of the legacy `test/integration/runtests.jl` so they can be reused by both the not-yet-migrated `test/suites/NN_*.jl` files and new `test/unit/` / `test/golden/` tests.

- `mock_models.jl` — mock density/effector/force-torque/planet models used to exercise `SimulationModel` interfaces without real physics.
- `test_constants.jl` — shared `EARTH`, `SPICE_PATH`, and small numeric helpers used across test files.
- `numeric_helpers.jl` — orbital-mechanics and harmonics numeric helpers (specific energy/angular momentum, angle unwrapping, regression slope, ICGEM/harmonics CSV writers).
- `spacecraft_builders.jl` — `SpacecraftModel` construction helpers (single/multi-link, AGORA Earth reference spacecraft).
- `config_builders.jl` — `SimulationConfiguration` builders (`build_config`, `build_config_multi`).
- `simulation_run_helpers.jl` — helpers that actually invoke `run_simulation` and read back results (`run_case`, `run_case_silent`, `run_case_capture_stdout`, `interp_linear`).
- `persistence_fixtures.jl` — `Solution` seeding helper for persistence/CSV-save tests.
- `sandbox_modules.jl` — throwaway modules used by architecture/contract tests that need an isolated module to include code into.
- `golden_harness.jl` — `run_and_compare_golden` and the PR/nightly tier gating (`golden_scenario_tier`, `golden_scenario_enabled`) used by `test/golden/`.
- `bootstrap.jl` — includes all of the above (in dependency order) plus the one-time module/engine setup (`SimulationModel`, `SimulationEngine`, GRAM extension re-patching, the legacy internals-alias block). Both `test/integration/runtests.jl` and `test/golden/runtests.jl` `include` this; it's guarded against double-inclusion in the same process.

New unit tests should `include(joinpath(REPO_ROOT, "test", "helpers", "<file>.jl"))` (or the equivalent relative path) for whichever of these they need, rather than redefining a mock model or builder inline. See `test/README.md` for the full convention.
