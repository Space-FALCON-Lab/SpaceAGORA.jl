# Earth periapsis corridor example

`examples/AGORA_Earth_Periapsis_Corridor.jl` studies a high-eccentricity Earth aerobraking orbit with a fixed physical periapsis-altitude corridor of **95,000–100,000 m**. The primary objective is to remain in that corridor for **60 days (5,184,000 s)** while keeping total delivered delta-v low. The maintenance policies are one-pass-lookahead heuristics: before each apoapsis they retarget the next periapsis to the configured target altitude using a perfect-model Cartesian predictor, and they are not optimizers of total delta-v. Burn count is a secondary comparison metric.

The example has two distinct uses: an independent low-e (~700 km altitude) vacuum diagnostic for J2/J3 mean-element behavior, and a high-e Cartesian aerobraking mission that evaluates the configured horizon. The vacuum diagnostic uses its own band, separate from the aerobraking 95–100 km corridor. The Vallado–Kwok map is a coordinate conversion; current vacuum minimum-radius comparisons do not establish a universal eccentricity accuracy domain or aerobraking safety. The aerobraking state and safety prediction use Cartesian calibration and physical radial events. The configuration defaults to 60 days (5,184,000 s). Shorter diagnostic runs can help inspect behavior, but do not establish full-mission success. Vacuum results do not establish atmospheric safety or aerobraking performance.

## Run the example

Run a reduced event/output smoke with:

```sh
julia --project=. examples/AGORA_Earth_Periapsis_Corridor.jl --smoke
```

This checks that the selected case starts, produces periapsis events, and writes expected output. It does **not** evaluate the 60-day objective or support a performance claim.

Run the default 60-day favorable passive-coast case with:

```sh
julia --project=. examples/AGORA_Earth_Periapsis_Corridor.jl
```

Run all five configured cases with `--case=all`:

```sh
julia --project=. examples/AGORA_Earth_Periapsis_Corridor.jl --case=all
```

Cases are `vacuum`, `favorable`, `unfavorable`, `near_passive`, and `frequent`. The vacuum case is a separate low-e, no-atmosphere diagnostic. Favorable and unfavorable are passive aerobraking coast phases. `near_passive` and `frequent` use the same calibrated unfavorable Cartesian initial state and differ in guidance policy. The runner accepts `--case=...`, `--config=PATH`, `--horizon-s=SECONDS`, `--output=PATH`, and `--smoke`; `--horizon-s` is a diagnostic override and a shortened horizon cannot satisfy the 60-day objective. `SPACEAGORA_EXAMPLE_SMOKE=1` also enables smoke mode. Smoke mode runs a reduced event/output check, not a performance comparison.

Results are written under `output/periapsis_corridor/<case>/`, with a combined `comparison.csv` at the selected output root. Each case writes `passes.csv`, `burns.csv`, `summary.csv`, and `run_manifest.toml`. The tables distinguish physical radial-event measurements, mean-coordinate diagnostics, Cartesian policy predictions, requested actions, and delivered actions. The summary marks full-mission verification only when elapsed time reaches at least 5,184,000 s, solver termination succeeds, periapsis event coverage is complete, and no physical corridor breach is recorded. Inspect termination cause and elapsed time before interpreting metrics: a truncated run or early apoapsis-target completion does not satisfy the default 60-day requirement.

## Initial condition and fixed corridor

The aerobraking candidate begins near 97,500 m physical periapsis altitude and 60,500 km apoapsis altitude, with a 50° inclination, 27° RAAN, 90° argument of periapsis, and 180° true anomaly. These are candidate initial conditions, not measured outcomes. Calibrate the initial Cartesian state against a physical radial-event pass and use Cartesian propagation for aerobraking prediction and safety checks. The configured vehicle has 620 kg dry mass, 200 kg propellant, and a 1,000 kg/m² ballistic coefficient; remaining actuator assumptions are configuration details, not evidence of performance. The fixed corridor is measured from the configured Earth reference radius using event-based physical radius. Mean or osculating periapsis estimates never substitute for physical periapsis from radial-velocity roots.

The atmosphere is a piecewise exponential fit to the 1976 U.S. Standard Atmosphere, with one scale height per layer. The layer boundaries are the table altitudes 80, 90, 95, 100, 110, 120, 130, and 150 km. Each layer is the exponential through the table densities at its two boundaries, `H = Δh / log(ρ_lower/ρ_upper)`, so density is continuous at every boundary and decreases monotonically with altitude. The 95 km and 100 km nodes are exactly 1.3873×10⁻⁶ and 5.6044×10⁻⁷ kg/m³, so the corridor band keeps the earlier single-layer scale height of about 5,692 m. The table rows (altitude and density, 80 to 150 km), their source, the node altitudes, and the out-of-range behavior are data in `[atmosphere]` of `examples/periapsis_corridor_config.toml`; nothing is hard-coded in Julia. Altitude is geodetic height above the planet ellipsoid, the same reference as the corridor and as the altitude the engine hands to GRAM.

Validity is 80 to 150 km. Outside that range the behavior is a configured choice, not an implicit one: `below_bottom_node` and `above_top_node` are each `"extend"` (continue the nearest layer's exponential) or `"zero"` (zero density). The defaults are `"extend"`. Above 150 km the extended top layer falls off faster than the real atmosphere, so far-from-periapsis density is an underestimate and is not an atmosphere-validity claim. The drag force honors both settings; the engine's own density model (`PiecewiseExponentialAtmosphereModel`) always extends its end layers. The fit error is measured at the seven table rows that are not layer boundaries (85, 105, 115, 125, 135, 140, 145 km): the worst-case relative error is **8.29%**, at 115 km, with 3.05% at 85 km and 0.18% at 105 km. Errors above 110 km are 4–8% and lie well above the corridor; between 95 and 105 km the fit is within 0.2%. This is a simplified deterministic model, not a GRAM realization or a statement about weather variability. See the [PDAS 1976 atmosphere table](https://www.pdas.com/bigtables.html) and the [NASA report record](https://ntrs.nasa.gov/citations/19770009539) for source context. Drag uses atmosphere-relative velocity. At physical periapsis events the run reports dynamic pressure, `q = ½ρ|v_rel|²`, and the drag-to-gravity Pi ratio, `Π_drag = q/(β μ/r²)`, with the configured ballistic coefficient `β`. These values are reevaluated as the orbit shrinks.

Gravity coefficients and reference parameters come from the repository's `data/planets/gravity/EarthGGM05C.csv`; the manifest records the exact source metadata. The vacuum diagnostic isolates axisymmetric J2/J3 gravity. In the aerobraking truth propagation, gravity, drag, and commanded finite burns act together.

## Mean-element convention and limits

The low-e vacuum diagnostic uses the **Vallado–Kwok first-order J2 short-period mean-to-osculating map**, adapted from `SatelliteToolboxPropagators.jl` 1.2.0. This is a coordinate conversion, not a full Brouwer–Lyddane J2–J5 implementation. It covers J2 short-period corrections only; it omits J2 long-period terms, J3 and higher zonals, drag, and third bodies. J3 and other forces remain in Cartesian truth propagation. The inverse's `.converged` field reports coordinate-inversion convergence only; `.corridor_safety_certified` is always `false`.

In analytical full-state comparisons, position and velocity targets were missed: the high-e case reported maximum vector differences of 98.262 km and 19.056 m/s, with sampled near-periapsis differences of 6.822 km and 11.158 m/s. The full-orbit `e=0.15` comparison reported 5.327 km and 5.464 m/s. Cartesian position error includes along-track phase error and is not periapsis-radius error. Separate central+J2 vacuum minimum-radius comparisons reported differences of −0.0709 m (high-e case), −4.266 m (low-e ~700 km case), and −1.076 m (`e=0.15`). These selected cases do not establish a general eccentricity accuracy domain or performance under J3, drag, or aerobraking. For aerobraking, use Cartesian calibration and prediction; physical radial events remain authoritative.

Primary provenance is the upstream [SatelliteToolboxPropagators project](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl), version 1.2.0, pinned at commit [`743bf1785697e4855ba829d1d73cde8ffd874670`](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl/tree/743bf1785697e4855ba829d1d73cde8ffd874670). The adapted equations come from its [`osculating.jl`](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl/blob/743bf1785697e4855ba829d1d73cde8ffd874670/src/propagators/osculating.jl) (SHA-256 `82af8571a54b918a65d5d2a3f7e758ffd8fca1c2e988cd926399ab063a34b229`); the upstream [`j2osc.jl` tests](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl/blob/743bf1785697e4855ba829d1d73cde8ffd874670/test/j2osc.jl) have SHA-256 `09c3c6fec098790d45af0397f57ea60b1cfbae45ae2970469bd7773c1ea5ba89`. The algorithm is attributed upstream to Kwok in Vallado, *Fundamentals of Astrodynamics and Applications*, 4th edition, pp. 708–710. The adapted formulas do not add a package dependency. The implementation's copyright and license notice are retained with the source.

This classical-element conversion is limited to elliptic, non-equatorial, non-circular input. A successful round trip shows coordinate inversion consistency; it does not prove that physical periapsis stays inside the corridor. That requires independent Cartesian propagation, complete radial-event detection, and a full 60-day run.

## Optional GRAM truth atmosphere

The default truth (and controller-predictor) atmosphere is the piecewise 1976 U.S. Standard Atmosphere fit above. `[atmosphere] backend` in `examples/periapsis_corridor_config.toml` selects it (`"us76_piecewise"`, the default) or `"gram"`; the flag `--atmosphere=gram` (or `--atmosphere=us76_piecewise`) overrides the file. Default behavior, tests, and the fixed 95,000-100,000 m corridor and 5,184,000 s horizon are unchanged.

```sh
julia --project=. examples/AGORA_Earth_Periapsis_Corridor.jl --case=favorable --atmosphere=gram --horizon-s=350000 --output=/tmp/gram_favorable
```

With `gram`, drag density is the Earth-GRAM mean density that the engine samples through the normal environment path (`environment.atmosphere.rho_kg_m3`); `CorridorBallisticDrag` keeps drag on the planet-relative velocity and the same constant-CdA ballistic-coefficient treatment. The Cartesian predictor that picks burns is built from the same configuration, so it uses the same backend as the truth propagation (the perfect-model assumption is unchanged and explicit). Because GRAM depends on date, the predictor, which restarts at t = 0 from the current state, is run at the run epoch plus the elapsed time at the maneuver request. The run epoch is the configured 2020-01-01 00:00:00 UTC; the manifest records the epoch GRAM was built with and the epoch the engine used (`gram_core_epoch_utc`, `gram_epoch_used_by_engine_utc`).

Limits and choices to read before using a GRAM result:

- **Mean density only, deterministic.** The dynamics read GRAM mean density; `gram_seed` and perturbation scales have no effect on trajectories. A GRAM run is not a stochastic or dispersed run. GRAM winds are not applied to drag (same zero-wind assumption as the fit).
- **Serial.** GRAM calls take a global lock; the case is not parallelized.
- **No fallback.** If GRAMSuite, the `SpaceAGORAGRAMSuiteExt` extension, the native `libGRAM`, or the SPICE kernels are missing, the run fails with a message saying so; it never substitutes the US76 fit. Direct point-to-point sampling is required: the GRAM track cache, vacuum cache, and per-step freeze environment switches make the run fail rather than interpolate.
- **Entry interface is a named field.** The GRAM extension evaluates GRAM only below `EnvironmentModel.EI` (above it the planet polynomial fit is used) and hard-codes zero density above 2000 km. `[atmosphere_gram] entry_interface_km = 2000.0` sets `EI` to that ceiling, so GRAM is sampled at every altitude where it is sampled at all and there is no entry-interface shortcut to another density source. The US76 backend keeps `EI = 160 km`, which only sets step-size regimes. A consequence is that the GRAM runs use the fine atmosphere step limit (`dt_max_atmosphere_s`) below 2000 km, so the two backends do not share a step schedule.
- **Altitude definition (geodetic; decided by the user, 2026-10-14).** The corridor (95,000 to 100,000 m), the target (97,500 m), the 250 m margin, the periapsis events' altitude, the drag fit's argument, the dynamic-pressure and Pi diagnostics, and the policy's predictor are all geodetic height above the planet ellipsoid, computed by the engine's own `rtolatlong` (Bowring) from the planet-fixed position; the ellipsoid is the planet model's `Rp_e`/`Rp_p`. GRAM already received that altitude, so GRAM and the fit use one definition and no density wrapper is needed. Periapsis detection is still the radial-velocity zero crossing; only the altitude conversion at the event changed. A corridor in radius space (eccentricity band, `Pi` context) is the latitude-dependent band `corridor_radius_band` (at 50 degrees the ellipsoid surface is about 12.6 km inside `Rp_e`). `spherical_periapsis_altitude_m` (`|r| - Rp_e`), the mean-map columns and `map_*_radius_m` are spherical-radius diagnostics only. The vacuum case has no corridor; its diagnostic band stays radial because it validates the J2 mean-to-osculating map in radius. The summary records `geodetic_helper_crosscheck_max_diff_m`, a startup check of the helper against the engine routine at the equator, 50 degrees and the pole.

GRAM output is written under `output/periapsis_corridor_gram/` by default, and each case writes `atmosphere_backend` into `summary.csv`, `passes.csv`, and `run_manifest.toml` (also `density_source`, `atmosphere_fit_used_for_drag`, and an `[atmosphere_provenance]` table with the GRAM directory, library path and SHA-256, GRAMSuite.jl version and git revision, epochs, and the sampling settings). A case directory that already holds a run from the other backend is refused. The `vacuum` case has no atmosphere and is skipped under `--atmosphere=gram`. The GRAM suite exposes no version string through the repository, so the directory name and library hash identify it.

### Measured short-horizon comparison (reduced horizon, not a mission result)

All rows below are `reduced_horizon_diagnostic` runs from a local comparison output directory (`comparison.md`, `comparison_per_pass.csv`, `comparison_runs.csv`, `density_profile.csv`, `timing_warm.csv`, `gram_vs_us76.png`), with the favorable and `frequent` cases at the same horizon and the same calibrated initial state (initial positions identical, velocities within 7e-5 m/s). None of them evaluates the 60-day objective.

| Case, horizon | Passes | Periapsis US76 / GRAM, first to last (m) | Density at periapsis US76 / GRAM (kg/m3) | Dynamic pressure US76 / GRAM (Pa) | Delivered delta-v US76 / GRAM (m/s) |
|---|---|---|---|---|---|
| favorable, 350,000 s | 5 | 97,500.52 to 97,491.77 / 97,500.30 to 97,498.91 | 8.82e-7 / 8.3e-8 (pass 1) | 46.6 / 4.4 | none (passive) |
| frequent, 840,000 s | 13 | 97,499.98 to 97,465.71 / 97,499.98 to 97,465.81 | 8.82e-7 to 8.96e-7 / 9.72e-7 to 9.89e-7 | 46.5 to 47.3 / 51.2 to 52.4 | 0.007694 / 0.007749 (one burn each) |

Density at the favorable periapsis (geodetic altitude for both): GRAM is 0.966 to 1.027 times the fit over the five-pass reduced-horizon favorable run (`density_over_fit_density` in `passes.csv`; it was 0.094 under the old spherical definition).

GRAM mean density divided by the fit at the same height number (GRAM given it as geodetic height, at the favorable pass-1 longitude and time; `density_profile.csv`):

| Latitude | 95 km | 97.5 km | 100 km |
|---|---|---|---|
| 0 | 1.112 | 1.110 | 1.102 |
| 25 degrees | 1.093 | 1.070 | 1.037 |
| 50.2 degrees (favorable periapsis) | 1.007 | 1.027 | 1.036 |
| -50.2 degrees | 1.084 | 0.986 | 0.915 |

So at 95-100 km GRAM differs from the fit by roughly -8.5 to +11 percent depending on latitude, which is larger than the 0.2 percent fit error against the table in that band but small against the factor 10 altitude-definition effect. This sample is one date, one longitude, and mean density; it is not a seasonal or dispersion statement.

Warm wall-clock (second run of each case in one process, so JIT compile is excluded; `timing_warm.csv`): favorable 350,000 s took 0.48 s (US76) and 3.85 s (GRAM); `frequent` 840,000 s took 10.2 s (US76) and 26.8 s (GRAM). The first run in a fresh process is dominated by compilation (about 37 to 40 s for the favorable run and 51 to 55 s for initial-periapsis calibration). GRAM density evaluations through the engine path were 287,741 (favorable, 5 passes) and 1,426,841 (`frequent`, 13 passes; 711,039 inside the 21 predictor runs); these count engine evaluations, some of which return zero above the 2000 km ceiling without a native call, so they bound the native call count from above. Extrapolating to 60 days is an estimate only: per-pass warm cost times the 74 to 109 passes a 60-day run is expected to contain (74 if the period stays at the initial 69,900 s, 109 as recorded for the earlier full US76 favorable run) gives about 7 to 10 s (US76 favorable), 57 to 84 s (GRAM favorable), 58 to 86 s (US76 `frequent`) and 153 to 225 s (GRAM `frequent`) of warm integration, plus extra predictor runs for every additional burn decision, plus a one-time compile of about 2 minutes. The 60-day GRAM run has not been executed.

## Inclination and argument-of-periapsis sweep

The base orbit never comes near the 95 km floor, so the controller comparison on it is non-discriminating. `examples/periapsis_corridor_sweep.jl` runs the same runner over a grid of inclination and argument of periapsis, with everything else as in the base configuration (RAAN, apoapsis, vehicle, epoch, `us76_piecewise`). Each cell re-calibrates the initial periapsis to 97.5 km geodetic with the runner's own calibration path and starts at true anomaly 180 deg (apoapsis), as the base case does; a cell whose calibration fails is recorded as `calibration_failed` with the reason, never skipped.

```bash
julia --project=. examples/periapsis_corridor_sweep.jl --phase=A --output=output/periapsis_corridor_sweep
julia --project=. examples/periapsis_corridor_sweep.jl --phase=B --margin-threshold-m=500
```

Flags: `--phase=A|B|both|plot`, `--atmosphere=us76_piecewise|gram`, `--cells=`, `--workers=` (GRAM only), `--inclinations-deg=30,50`, `--omegas-deg=0,90`, `--margin-threshold-m`, `--horizon-s` (a reduced horizon is labeled `reduced_horizon_diagnostic`), `--calibration-max-iterations`, `--output`, `--config`. Defaults come from the `[sweep]` table of the config: inclinations 30, 50, 63.4, 80, 98, 116.6 deg and omega 0 to 315 deg in 45 deg steps. 63.4 and 116.6 deg are the critical inclinations (5 cos^2 i = 1, i = 63.435 and 116.565 deg) where J2 apsidal drift vanishes (derived); the other values and the 500 m threshold are assumed.

Phase A runs the passive coast (runner case `favorable`, no controller) for every cell in one Julia process and writes `passive_grid.csv` and `passive_margin_heatmap.png` (minimum margin to the 95 km floor; breaches marked LOW or HIGH, below 95 km or above 100 km). Phase B reads that grid and runs the `frequent` policy (100 m deadband, same limits as the `frequent` case) on every cell that breached or whose passive margin to the floor is below the threshold, at the cell's own argument of periapsis (the stand-alone `frequent` case forces omega = 0 through its unfavorable phase; the sweep keeps the configured omega). It writes `controlled_grid.csv` (delivered delta-v, burn count, propellant, `full_mission_verified`, controller failures, `feasible`/`infeasible` with cause) and `controlled_dv_heatmap.png`. Per-cell runner outputs are under `cells_passive/` and `cells_controlled/`. The cell (50, 90) reproduces the base `favorable` case and (50, 0) the `unfavorable` case.

GRAM backend. `--atmosphere=gram` runs the same cells with Earth-GRAM mean density (default output `output/periapsis_corridor_sweep_gram`; every CSV row carries a `backend` column, a `sweep_manifest_<phase>.toml` records the backend, and a Phase B or plot run refuses a `passive_grid.csv` from another backend). `--cells=i:omega,i:omega` restricts a run to explicit cells (Phase B with `--cells` runs exactly those cells instead of applying the margin rule). GRAM cells go through `run_monte_carlo(...; threads=:auto)` forced onto the process route (`--workers=N`, default 4): native GRAM is not thread-safe, so each worker is its own `--threads=1` OS process with its own GRAM instance, and the driver stops with an error if the result's `route` is not `:process`. A GRAM worker held 6 to 8 GB resident in Phase B, so 4 workers were killed by a 16 GB cap and Phase B was rerun with `--workers=2`; size `--workers` to memory. Example: `julia --threads=1 --project=. examples/periapsis_corridor_sweep.jl --atmosphere=gram --phase=A --cells=30:225,98:315`.

The GRAM sweep used the mean density only (no GRAM perturbations enter the dynamics), is deterministic, and uses a single epoch. On the 10 cells that breach on US76, GRAM breaches on the same 10, the `frequent` policy holds all 10 (0.222 to 0.410 m/s), and the two smallest-margin US76 passes (30, 270) and (50, 135) also pass on GRAM (margins 479 and 625 m). Comparison files: `comparison.csv` and `dv_us76_vs_gram.png` in the GRAM output directory.

Limits: Deterministic single date, mean density, one ballistic coefficient. A passive pass means the 60-day record stayed in the corridor, not that it is robust to dispersions. The corridor and horizon are never widened.

## Interpreting results

The current preliminary favorable passive-coast output records a full 5,184,000 s run with successful termination, 109 periapsis events, zero corridor breaches, minimum/maximum physical periapsis of 97,369.061941/97,500.518164 m, and zero burns / 0 m/s delivered delta-v. Its summary marks full-mission verification. For this configured initial orbit and deterministic model, it demonstrates a feasible zero-maintenance-delta-v coast, which is the lower bound on maintenance delta-v. The run predates final mean-output updates and convergence checks; insertion cost, robustness, and other initial conditions are separate questions.

The 60-day result is feasible only if every physical periapsis event remains in the fixed 95,000–100,000 m corridor and the integration completes the configured elapsed duration. If no candidate or control strategy meets both conditions, report infeasibility with the termination cause and extrema; do not widen the corridor, weaken the atmosphere, or substitute a shorter mission.

Among feasible comparisons, total delivered delta-v is the comparison metric (the policies do not minimize it). Burn count and cadence describe the strategy but are not independent success gates. Delivered delta-v is derived from actuator/mass history rather than guidance requests alone. No performance or savings claim should be made until the full-duration comparison and numerical convergence checks have completed.

In the vacuum case, `passes.csv` also records, per pass, `map_periapsis_radius_m` (the J2 short-period map evaluated at mean true anomaly 0 from that pass's mean elements), `physical_periapsis_radius_m` (the radial-event radius), and their signed difference `map_minus_physical_periapsis_radius_m`; `map_vs_physical_periapsis_radius.png` plots the difference. The J2 map is a diagnostic only and is not validated for the high-e cases.

The separate low-e, ~700 km vacuum case checks mean-element behavior and physical event diagnostics in its own diagnostic band. It is not the high-e vacuum equivalent of the aerobraking initial state and cannot be used to infer the 60-day atmospheric result.

### Comparing the `near_passive` and `frequent` policies in the sweep

`--phase=B --policy=near_passive` runs the same cells with the `near_passive` controller (burn only when the predicted next periapsis leaves the band tightened by `safety_margin_m`, then retarget to the reference altitude) and writes `*_near_passive` files, so outputs never mix with `frequent`; every CSV row has a `policy` column. On the 10 passively breaching cells it held the corridor in all 10 on both US76 and GRAM mean density with 11 to 12 burns in total, against 281 to 284 for `frequent`, with a total delivered delta-v of 3.04 m/s (US76) and 2.73 m/s (GRAM) against 3.12 and 3.07 m/s. It leaves the periapsis within 300 m of a corridor edge by construction. The result is deterministic, mean density, single epoch; the controller is a one-pass-lookahead heuristic with a perfect-model predictor, not an optimizer.
