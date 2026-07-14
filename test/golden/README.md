# Golden Propagation Regression Tests

Reproducible-propagation tests: each scenario runs a full `run_simulation`
and compares the resulting trajectory against a checked-in fixture within a
per-point tolerance. See `test/TEST_RESTRUCTURE_PLAN.md` for why this
directory exists and its migration history.

Run all of them:

```bash
julia --project=. test/golden/runtests.jl
```

By default only `:pr`-tier scenarios run. To also run `:nightly`-tier
scenarios (currently `gram_mars`, which is slower and depends on GRAMSuite):

```bash
SPACEAGORA_GOLDEN_TIER=all julia --project=. test/golden/runtests.jl
```

## Layout

Each scenario is a directory `test/golden/<scenario_name>/` containing:

- `config.jl` — defines a zero-argument function `build_<scenario_name>_config()`
  returning a `SimulationConfiguration`. May reference anything set up by
  `test/helpers/bootstrap.jl` (spacecraft builders, `EARTH`/`Mars`/... planets,
  `SPICE_PATH`, etc.).
- `fixture.csv` — the checked-in golden trajectory. Columns: `time`,
  `pos_atol_m`, `vel_atol_mps`, then `sc<k>_pos_<1|2|3>` / `sc<k>_vel_<1|2|3>`
  for each spacecraft `k` (1-indexed, in the order passed to
  `build_config`/`build_config_multi`).
- `metadata.toml` — `description`, `tier` (`"pr"` or `"nightly"`), and
  provenance (`git_sha`, `julia_version`, `generated_at`) written by
  `scripts/regenerate_golden.jl`.

## Adding a new scenario

1. Create `test/golden/<name>/config.jl` defining `build_<name>_config()`.
2. Write `test/golden/<name>/metadata.toml` with `description` and `tier`
   (provenance fields get filled in by the regen script).
3. Generate the fixture:
   ```bash
   julia --project=. scripts/regenerate_golden.jl <name> <spacecraft_count>
   ```
4. Add a `@testset` block for it in `test/golden/runtests.jl` calling
   `run_and_compare_golden("<name>", build_<name>_config; spacecraft_count=<n>)`.

## Regenerating a fixture

Only regenerate `fixture.csv` when the trajectory changing is an *intentional*
consequence of a model/physics change — treat a diff to a fixture as needing
extra review scrutiny, not a routine update. If a scenario's `config.jl`
didn't change and the fixture still doesn't match, that's a regression to
investigate, not a reason to regenerate.

```bash
julia --project=. scripts/regenerate_golden.jl <scenario_name> <spacecraft_count>
```
