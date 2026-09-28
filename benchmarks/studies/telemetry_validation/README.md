# Telemetry validation: Odyssey and Venus Express aerobraking reconstructions

This study is the runnable form of the internal technical record
[`docs/spaceagora_aerobraking_reconstruction_record.md`](../../../docs/spaceagora_aerobraking_reconstruction_record.md)
(also `docs/spaceagora_aerobraking_reconstruction_record.pdf`). It executes
the Mars Odyssey and Venus Express (VEX) aerobraking reconstructions exactly
as documented there, through the `SpaceAGORA.TelemetryVerification` engine,
with one frozen manifest per run under `manifests/`.

The manifests are copies-of-record: `odyssey_tolson` and `vex_venusgram`
mirror the `odyssey`/`vex` scenarios in `test/telemetry_benchmark_manifest.toml`
(the nightly/PR sentinel path); the MarsGRAM comparator and the ±1σ envelope
exist only here. If the record configuration ever changes, both places must
change together.

## Runs

| Run | Atmosphere | GRAM needed | Full-profile runtime* | Role |
|---|---|---|---|---|
| `odyssey_tolson` | Tolson accelerometer-derived density (tabulated flight) | no | ~53 min | **Benchmark of record** — the 0.667% nMAE result |
| `odyssey_marsgram` | MarsGRAM TES Mapping Year 2 (MY25), deterministic | yes | ~4.3–5.8 h | Climatology comparator (context, not a twin gate) |
| `odyssey_tolson_sigma_minus` | Tolson density at −1σ | no | ~53 min | Density-uncertainty envelope (informational) |
| `odyssey_tolson_sigma_plus` | Tolson density at +1σ | no | <53 min (impacts early) | Density-uncertainty envelope (informational) |
| `vex_venusgram` | VenusGRAM, deterministic | yes | ~18 min | VEX 2014 reconstruction |

\* Recorded on the record's Apple arm64 system; expect different wall times here.

Recorded full-profile values to compare against:

| Quantity | Odyssey: Tolson | Odyssey: MarsGRAM | VEX: VenusGRAM |
|---|---|---|---|
| Apoapsis MAE / RMSE / max (km) | 162.60 / 179.70 / 320.95 | 1399.52 / 1757.47 / 3445.60 | 254.34 / 397.93 / 948.24 |
| Apoapsis normalized MAE | 0.00667 | 0.05745 | 0.0701 |
| Periapsis MAE / RMSE / max (km) | 4.16 / 5.77 / 13.50 | 2.53 / 3.17 / 9.78 | 0.422 / 0.617 / 2.965 |
| Drag-decay ratio (sim/flight) | 0.984 | 0.840 | 1.262 median / 1.283 span-weighted |
| Event alignments (apo / peri, km) | −15.422 / −9.851 | −153.010 / −9.857 | −12.585 / +0.253 |
| Fitted Cr | pinned 1.3 | pinned 1.3 | 1.15 (grid {1.15, 1.30, 1.45}) |

The apoapsis worst-point metric is platform-sensitive (record §6: max_abs
320.95 km on arm64 vs 938.57 km on x86_64 with ~13% mission-duration change);
MAE/nMAE are the stable comparison signals.

## Commands

```bash
# Record scenarios (Tolson needs no GRAM; VEX needs the native GRAM build):
julia --project=. benchmarks/studies/telemetry_validation/run_validation.jl full --runs=core

# Everything in the record (adds the MarsGRAM comparator and the ±1σ envelope):
julia --project=. benchmarks/studies/telemetry_validation/run_validation.jl full --runs=all --enforce=true

# Fast pipeline check (first orbits only; does NOT reproduce record values):
julia --project=. benchmarks/studies/telemetry_validation/run_validation.jl quick --runs=odyssey_tolson

# Figure 1 input: matched accelerometer-vs-MarsGRAM density comparison
# (needs native GRAM + the m01_ab_v2.bsp kernel under the GRAM Suite SPICE tree):
julia --project=. benchmarks/studies/telemetry_validation/density_compare.jl

# Figures 1-5 of the record's visual section, from whatever outputs exist:
julia --project=. benchmarks/studies/telemetry_validation/plot_results.jl
```

`--enforce=true` applies the manifest gate tolerances to `odyssey_tolson`,
`odyssey_marsgram`, and `vex_venusgram` only; the ±1σ envelope runs are never
enforced (the +1σ case impacts before the end of the campaign, so its
full-span metrics are not comparable — record §6).

Outputs land in `results/<hostname>/` (gitignored):

```
results/<hostname>/
  <run>/telemetry_orbit_accuracy_summary.csv   # per-event metrics, calibration, gates
  <run>/telemetry_orbit_accuracy_errors.csv    # per-orbit truth/sim/error rows
  density_compare/odyssey_accel_vs_marsgram_density.csv
  figures/figure{1..5}_*.png
```

## What each run does (and does not) show

These are **calibrated mission reconstructions**, not blind predictions and
not component-by-component validation (record §1, §7, §8):

- No pass or time interval is withheld. The constant per-event altitude
  alignments (`b = -median(e_1..e_10)`, ±500 km cap — a candidate at the cap
  is a model-mismatch alarm and is reported unapplied) are estimated from the
  same event series that is scored.
- `odyssey_tolson` largely removes atmospheric-climatology error, so its
  0.667% nMAE evaluates the aggregate propagation + free-molecular aero +
  maneuver replay + apsis-detection chain. The Tolson density itself embeds a
  flight-side aerodynamic model, so the run is only weakly independent of
  aerodynamic assumptions. Cr is pinned at 1.3 and Cd scale at 1.0 by
  deliberate decision of record (an earlier Cr search selected 1.45 and
  absorbed atmospheric over-drag).
- `odyssey_marsgram` shows what restoring the selected climatology does to
  the same chain (84.0% vs 98.4% recovered decay). It is not blind MarsGRAM
  validation: the TES family was selected using Odyssey density data from
  this same campaign.
- `vex_venusgram` exercises Venus-specific geometry, element-frame handling,
  maneuver replay, and geodesy against plot-digitized truth. It fits Cr on
  the scored campaign and uses the legacy comparison axis (no
  `epoch_orbit_offset`), which adds ~0.5–1 orbit registration uncertainty to
  the apoapsis budget. Its ~26–28% apoapsis over-decay is an aggregate
  drag-corridor discrepancy; VEX alone cannot split atmosphere from
  aerodynamics.

## Relationship to the nightly sentinel

`benchmarks/studies/telemetry_orbit_accuracy_study.jl` +
`test/telemetry_benchmark_manifest.toml` remain the CI/nightly path (the
Tolson digital-twin regression sentinel). This directory exists to reproduce
the *complete* record — including the MarsGRAM comparator and the ±1σ
envelope, which are deliberately not part of the nightly — without touching
the sentinel configuration.

## Data prerequisites

- `data/telemetry/Odyssey/*.feather`, `data/telemetry/VEx/*.feather` — in-tree.
- `data/Gravity_harmonics_data/{Mars50c,MGNP180U}.csv` — in-tree.
- Native GRAM build under `data/GRAMSuite.jl/GRAM Suite 2.0` (see
  `scripts/ensure_gram_native.jl`) — required by `odyssey_marsgram`,
  `vex_venusgram`, and `density_compare.jl` only.
- `data/GRAMSuite.jl/GRAM Suite 2.0/SPICE/spk/missions/m01_ab_v2.bsp` —
  required by `density_compare.jl`.
