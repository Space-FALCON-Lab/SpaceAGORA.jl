# Troubleshooting

Use this page when something is not working and you want a quick diagnostic
reference.

What to read next:

- [GRAMSuite Setup](gramsuite_setup.md)
- [Assets & Modes](../assets.md)
- [Simulation Configuration](simulation_configuration.md)

---

## First run is very slow

**Symptom:** The first `julia --project=. ...` invocation takes several minutes
before anything prints.

**Cause:** Julia compiles all package code on first load. SpaceAGORA pulls in
ODE solver, atmosphere, and geometry libraries that take time to precompile.
Subsequent runs in the same Julia session or after the compiled cache is
populated are significantly faster.

**Resolution:**
- This is expected behavior. Wait for the first run to complete.
- The compiled cache persists across sessions, so only the very first run per
  Julia install is slow.
- To warm the cache explicitly without running a full simulation:
  ```text
  julia --project=. -e "using SpaceAGORA"
  ```

---

## GRAM assets not found

**Symptom:** `julia --project=. src/cli/main.jl assets check` reports `gram_root` or
`spice_directory` as missing, or a GRAM-backed example errors on load.

**Cause:** The official NASA GRAM Suite files have not been copied into the
expected local path.

**Resolution:** Follow the full setup walkthrough in
[GRAMSuite Setup](gramsuite_setup.md). The key steps are:

1. Pull the `GRAMSuite.jl` submodule:
   ```text
   git submodule update --init --recursive --remote
   ```
2. Copy the official GRAM Suite distribution into:
   ```text
   data/GRAMSuite.jl/GRAM Suite 2.0
   ```
3. Build the native shared library (needs a C/C++ toolchain and GNU Make; no
   Fortran compiler, and CSPICE is bundled on Linux x86_64 and Windows):
   ```text
   julia --project=. scripts/ensure_gram_native.jl
   ```
4. Confirm with:
   ```text
   julia --project=. src/cli/main.jl assets check
   ```

If the submodule checkout fails with Git LFS `no space left on device`, free
disk space and see [GRAMSuite Setup](gramsuite_setup.md#git-lfs-reports-no-space-left-on-device).

---

## Native GRAM library missing or stale

**Symptom:** GRAM assets are present but `GRAMAtmosphereModel` raises a library
load error at runtime.

**Cause:** The platform-native `libGRAM` was either never built on this host, or
the GRAM tree was copied in from another machine with its `Build/lib` already
populated.

**Resolution:**

```text
julia --project=. scripts/ensure_gram_native.jl
```

If that prints `Native GRAM library already present for this host` and the load
error persists, the library on disk belongs to a different machine. The
existence check cannot distinguish it from a native build, so force a rebuild:

```text
julia --project=. scripts/ensure_gram_native.jl --clean
```

This is the common case on remote hosts, because
`scripts/remote/spaceagora-remote` mirrors the entire working tree, native
build outputs included.

A checkout that was moved or renamed on the same machine does not need
`--clean`; the vendored build helper detects the changed root from its build
manifest and rebuilds on its own.

Build prerequisites and the per-platform CSPICE story are covered in
[GRAMSuite Setup](gramsuite_setup.md#build-the-native-gram-library).

---

## Solver stalls or takes extremely long during entry

**Symptom:** The integrator slows to a crawl during an aerobraking or entry
pass, producing thousands of rejected steps.

**Cause:** The default `dt_max_atmosphere` (1 s) may be too large for a steep
entry with rapidly-changing aerodynamic forces, causing the adaptive solver to
spend time in step-rejection loops.

**Resolution:** Reduce `dt_max_atmosphere` in `IntegrationTolerances`:

```julia
SM.IntegrationTolerances(
    dt_max_atmosphere = 0.1,    # or 0.5 for moderate entries
    reltol_atmosphere = 1e-9,
    abstol_atmosphere = 1e-11
)
```

For very steep entries (ballistic coefficient > 500 kg/m², flight path angle
steeper than −15°), start with `dt_max_atmosphere = 0.05`.

---

## Output CSV is not written

**Symptom:** The simulation completes but no CSV appears in the output
directory.

**Cause:** Either `simulation_settings.results = false` or
`simulation_settings.save_csv = false`.

**Resolution:** Check the `SimulationSettings` in your configuration:

```julia
SM.SimulationSettings(
    results   = true,
    save_csv  = true,
    results_directory = "output"
)
```

If using `make_example_config`, the `results` and `results_directory`
arguments control this:

```julia
make_example_config(
    ...
    results           = true,
    results_directory = "output/my_run"
)
```

---

## Concurrent runs overwrite each other's output

**Symptom:** Parallel runs write to the same CSV or produce a file named
`simulation_results_collision_*.csv`.

**Cause:** Multiple `run_simulation` calls are writing to the same
`results_directory`.

**Resolution:** Give each run a unique `results_directory`:

```julia
for (i, config) in enumerate(configs)
    config_i = SM.SimulationConfiguration(
        config;
        simulation_settings=SM.SimulationSettings(
            config.simulation_settings;
            results_directory = "output/run_$(i)"
        )
    )
    run_simulation(config_i)
end
```

The collision file (`*_collision_*`) is written as a safety net when two
concurrent writes to the same path are detected.

---

## `isolate_state` and mutable-state aliasing

**Symptom:** Repeated calls to `run_simulation` with the same configuration
object produce different results, or the configuration appears modified after
the first run.

**Cause:** `run_simulation` deep-copies the configuration by default
(`isolate_state=true`), but some user code is passing `isolate_state=false` for
performance and then reusing a configuration object that has been mutated.

**Resolution:** Either rely on the default `isolate_state=true` (safe,
slightly slower) or ensure you never reuse a configuration instance across
overlapping or state-mutating calls when using `isolate_state=false`:

```julia
run_simulation(config)                        # safe: deep-copies config
run_simulation(config; isolate_state=false)   # fast: config must not be shared or reused
```

---

## `NRLMSISE00AtmosphereModel` network errors

**Symptom:** `NRLMSISE00AtmosphereModel(use_space_indices=true)` errors on the
first atmosphere evaluation with a network or download error.

**Cause:** `SpaceIndices` (via `SatelliteToolbox`) needs to download CelesTrak
space weather data on first use, and the download failed.

**Resolution:** Run the prewarm step manually with network access before
submitting to a compute cluster:

```julia
init_nrlmsise_space_indices!()
```

Or use fixed indices instead:

```julia
NRLMSISE00AtmosphereModel(f107a=150.0, f107=150.0, ap=4.0)
```

---

## `GRAMAtmosphereModel` not defined

**Symptom:** `GRAMAtmosphereModel` raises `UndefVarError` even though GRAM
assets are present.

**Cause:** The `SpaceAGORAGRAMSuiteExt` package extension only loads after
`GRAMSuite` is imported into the session.

**Resolution:** Call `setup_gram_example!()` (from `examples/common.jl`) or
manually import `GRAMSuite` before using `GRAMAtmosphereModel`:

```julia
# Run from the repository root:
include(joinpath(pkgdir(SpaceAGORA), "examples", "common.jl"))
setup_gram_example!()
# GRAMAtmosphereModel is now available via SM.GRAMAtmosphereModel
```

---

## Verification study threshold failures

**Symptom:** `julia --project=. src/cli/main.jl telemetry quick --enforce=1` exits non-zero and
reports threshold violations.

**Cause:** Either the simulation output diverges from the reference telemetry,
or tolerance thresholds are set tighter than the current model accuracy.

**Checklist:**

1. Confirm the run completed without solver errors (check logs for `retcode`).
2. Check `output/telemetry_cli/` for the per-metric report.
3. If you changed integration tolerances, reset them to the defaults in
   `make_example_config` and re-run.
4. If GRAM is involved, verify the asset layout with
   `julia --project=. src/cli/main.jl assets check`.
