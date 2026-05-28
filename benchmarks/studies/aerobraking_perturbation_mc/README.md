# Aerobraking Perturbation Short-Arc Sweep

This study runs a deterministic, paired short-arc sensitivity sweep for
aerobraking dynamics perturbations. The default grid spans Mars, Earth, Venus,
and Titan:

- planets: `mars`, `earth`, `venus`, `titan`
- periapsis regimes: `shallow`, `nominal`, `deep`
- apoapsis altitudes (km, per planet):
  - mars: 1 000, 1 250, 1 500, 1 750, 2 000, 2 500, 3 000, 4 000, 5 000, 7 500, 10 000, 12 000, 15 000, 20 000, 30 000
  - earth: 1 000, 1 500, 2 000, 3 000, 5 000, 7 500, 10 000, 15 000, 20 000, 30 000, 36 000, 45 000, 60 000
  - venus: 1 000, 1 500, 2 000, 3 000, 5 000, 7 500, 10 000, 15 000, 20 000, 30 000, 40 000
  - titan: 1 000, 1 500, 2 000, 3 000, 5 000, 7 500, 10 000, 15 000, 20 000, 30 000, 40 000
- dynamics cases: `two_body`, `j2`, `harmonics_low`, `srp`, `third_body_sun`,
  `gram_aero`, `full_environment`
- density cases for aero cases only: `low = 0.9`, `nominal = 1.0`,
  `high = 1.1`
- inclination sweep (deg): 15, 30, 45, 60, 75, 93, 105, 120, 135, 150, 165
- argument-of-periapsis sweep (deg): 0, 30, 45, 60, 80, 90, 120, 135, 150, 180, 210, 240, 270, 300, 330
- propagation length: `1` initial orbital period

Run the default sweep:

```bash
julia --project=. benchmarks/studies/aerobraking_perturbation_mc/main.jl
```

Run with process workers:

```bash
julia --project=. benchmarks/studies/aerobraking_perturbation_mc/main.jl --procs 4
```

Workers run `GC.gc(true)` after every case and, on Linux, make a best-effort
`malloc_trim(0)` call so freed allocator arenas can be returned to the OS.
By default, process workers are also recycled after 16 cases or whenever their
post-GC RSS is above 6 GiB, keeping slow long-run memory growth bounded while
avoiding per-case startup cost. Tune this with `--worker-max-cases N` and
`--worker-max-rss-gb GB`, or set either value to `0` to disable that trigger:

```bash
julia --project=. benchmarks/studies/aerobraking_perturbation_mc/main.jl \
  --procs 8 --worker-max-cases 8 --worker-max-rss-gb 6
```

The study also raises the outside-atmosphere solver step cap to 100 seconds by
default, matching the outside-atmosphere output cadence. Override it with
`--dt-max-orbit SECONDS` when doing accuracy/speed sweeps. The inside-atmosphere
solver step cap defaults to 5 seconds, matching the inside-atmosphere output
cadence; override it with `--dt-max-atmosphere SECONDS`.

Aero cases default to `--aero-solver auto_stiff_then_rodas5p`: try
`AutoTsit5(Rodas5P)` first, then retry only stiff-like failures such as
`MaxIters` with full `Rodas5P`. Use `--aero-solver rodas5p` to force the stiff
solver for every `gram_aero` and `full_environment` case, or
`--aero-solver auto_stiff` to disable the retry path.
The auto-stiff attempt defaults to `--aero-auto-maxiters 500000`, while the
Rodas fallback keeps `--aero-stiff-maxiters 20000000`. Lower
`--aero-auto-maxiters` to fail over sooner; tune `--aero-auto-stiff-switch-max`
to adjust OrdinaryDiffEq's internal AutoTsit5 stiff-switch threshold. Rodas5P
aero attempts use `--aero-stiff-tol-scale 10` by default, multiplying the study
tolerances by 10 to make the fallback less strict.

Aero cases also have two robustness guardrails. `--case-timeout-min 30`
terminates a single case after 30 minutes of wall-clock time; in distributed
runs the parent process recycles the timed-out worker and continues. And
`--deorbit-bailout` stops an aero case once its osculating apoapsis falls below
the atmosphere interface plus 25 km. Set `--case-timeout-min 0` or
`--no-deorbit-bailout` to disable those guards.

Run a tiny verification case:

```bash
julia --project=. benchmarks/studies/aerobraking_perturbation_mc/main.jl --smoke --norbits 1
```

Override the apoapsis grid (applied to all planets):

```bash
julia --project=. benchmarks/studies/aerobraking_perturbation_mc/main.jl \
  --apoapsis-alts 1000,2000,4500
```

Apoapsis altitudes below 1 000 km are dropped from study inputs.

Outputs are written to `output/aerobraking_perturbation_mc/<timestamp>/`.
Large tables are written as both `.csv` and `.feather`:

- `results.csv` / `results.feather`: one row per simulation
- `paired_deltas.csv` / `paired_deltas.feather`: deltas against the matching
  `two_body` baseline
- `aggregates.csv` / `aggregates.feather`: grouped summary table
- `manifest.toml`: run metadata and grid definition
- `case_*_<planet>_<regime>_<case>/orbit_chunk_manifest.feather`: one row per
  per-orbit trajectory chunk
- `case_*_<planet>_<regime>_<case>/trajectory_with_active_force_orbit_###.feather`:
  trajectory rows plus active perturbation force columns for one orbit
- `case_*_<planet>_<regime>_<case>/simulation_results.feather`,
  `active_force_history.feather`, and `trajectory_with_active_force.feather`:
  combined single-file outputs rebuilt from the orbit chunks after each case

The `harmonics_low` case uses 20x20 gravity harmonics for Mars, Earth, and
Venus. Titan uses `data/Gravity_harmonics_data/titan5.csv`, capped at degree and
order 5. Titan third-body cases include Saturn and the Sun.

Generate perturbing-force ratio summaries and plots after a run:

```bash
julia --project=. benchmarks/studies/aerobraking_perturbation_mc/plot_perturbation_force_ratios.jl \
  --rebuild-summary
```

By default the plotting pass reads the most recent timestamped dataset under
`output/aerobraking_perturbation_mc/`. Pass a run directory explicitly only when
you want to plot an older dataset.

The plotting pass writes `plots/perturbation_force_ratio_summary.feather`,
`plots/perturbation_force_ratio_summary.csv`, a filtered summary CSV, per-planet
time-history PDFs, cross-planet heatmaps/comparison plots, and ranking tables.
The force ratio is
`active_perturbation_force_mag / (sc1_mass * GM / norm(sc1_pos)^2)`. Useful
plotting options include `--metric peak|p95|p50|max_in_atmosphere`,
`--plot-set all|time|heatmaps|rankings`, and
`--density-case nominal|low|high|all` (default: `nominal`). Time-history plots include analytical
basic/detailed perturbation-parameter overlays by default; use
`--analytical-overlays none|basic|detailed|both` to control them. The plotting
pass also writes per-planet `analytical_basic_vs_detailed_*.pdf` comparison
plots for the selected time-history cuts, with nominal-density aero only and separate J2, higher-degree
harmonics, SRP, third-body, and aerodynamic analytical curves. Summary generation uses multiple threads
by default: if the plotting script is launched from a single-threaded Julia
process, it restarts itself with `--threads=auto`. Cap summary worker tasks with
`--summary-threads N` or `SPACEAGORA_AERO_PERTURB_SUMMARY_THREADS=N` when memory
bandwidth is tight. Set `SPACEAGORA_AERO_PERTURB_NO_THREAD_REEXEC=1` to disable
the automatic restart, or `SPACEAGORA_AERO_PERTURB_AUTO_THREADS=N` to choose a
different restart thread count.

Per-case trajectory output is streamed by orbit to limit memory consumption.
The study saves approximately one point every 100 seconds outside the
atmosphere and every 5 seconds through the atmosphere passage, then flushes the
current orbit chunk to Feather before continuing. After the solve, chunks are
combined one case at a time into the familiar single-file Feather outputs for
post-processing. The normal `simulation_results.csv` path is disabled for this
study so 8 process workers can run on a 64 GB machine without retaining full
multi-orbit trajectories during propagation.

GRAM density is evaluated through the vacuum-predicted trajectory cache for
this study. The driver turns on `SPACEAGORA_VACUUM_GRAM_CACHE`, sets a
24-point, 900-second look-ahead, uses a 5 km absolute inertial position
deviation threshold for cache rebuilds, and disables the native GRAM track cache
before launching worker processes. Override these study defaults with the
matching `SPACEAGORA_AERO_PERTURB_*` environment variable names when comparing
against direct point sampling.

GRAM density columns that require runtime-only solver state are reported as
`NaN` when they cannot be reconstructed after the solve; orbital-element deltas
and final-state metrics are still computed from the returned trajectory.
