# Aerobraking Perturbation Short-Arc Sweep

This study runs a deterministic, paired short-arc sensitivity sweep for
aerobraking dynamics perturbations. The default grid spans Mars, Earth, and
Venus:

- planets: `mars`, `earth`, `venus`
- periapsis regimes: `shallow`, `nominal`, `deep`
- apoapsis altitudes (km, per planet):
  - mars: 300, 500, 750, 1 000, 1 500, 2 000, 3 000, 4 500, 7 000, 12 000
  - earth: 600, 1 000, 2 000, 5 000, 10 000, 20 000, 36 000, 60 000
  - venus: 600, 1 000, 2 000, 5 000, 10 000, 20 000, 40 000
- dynamics cases: `two_body`, `j2`, `harmonics_low`, `srp`, `third_body_sun`,
  `gram_aero`, `full_environment`
- density cases for aero cases only: `low = 0.75`, `nominal = 1.0`,
  `high = 1.25`
- propagation length: `30` initial orbital periods

Run the default sweep:

```bash
julia --project=. benchmarks/studies/aerobraking_perturbation_mc/main.jl
```

Run with process workers:

```bash
julia --project=. benchmarks/studies/aerobraking_perturbation_mc/main.jl --procs 4
```

Run a tiny verification case:

```bash
julia --project=. benchmarks/studies/aerobraking_perturbation_mc/main.jl --smoke --norbits 1
```

Override the apoapsis grid (applied to all planets):

```bash
julia --project=. benchmarks/studies/aerobraking_perturbation_mc/main.jl \
  --apoapsis-alts 500,1000,2000,4500
```

Outputs are written to `output/aerobraking_perturbation_mc/<timestamp>/`.
Large tables are written as both `.csv` and `.feather`:

- `results.csv` / `results.feather`: one row per simulation
- `paired_deltas.csv` / `paired_deltas.feather`: deltas against the matching
  `two_body` baseline
- `aggregates.csv` / `aggregates.feather`: grouped summary table
- `manifest.toml`: run metadata and grid definition
- `case_*_<planet>_<regime>_<case>/simulation_results.csv` /
  `simulation_results.feather`: regular trajectory data for that simulation
- `case_*_<planet>_<regime>_<case>/active_force_history.csv` /
  `active_force_history.feather`: time history of the active perturbation force
  vector and magnitude
- `case_*_<planet>_<regime>_<case>/trajectory_with_active_force.csv` /
  `trajectory_with_active_force.feather`: regular trajectory data plus the
  active perturbation force columns

Regular trajectory CSV saving is enabled by default because the aerodynamic
force history reuses the saved drag, lift, and crosswind channels. Use
`--no-save-simulation-csv` only for focused non-aero debugging.

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
