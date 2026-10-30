# Periapsis border switchability study

Purpose: test whether a small maneuver can switch an Earth aerobraking orbit between periapsis-raising,
static and lowering regimes under Sun and Moon third-body gravity, and whether that is cheaper than stepping
the periapsis directly. The metric is the per-pass peak free-stream heat rate `Q = 0.5 rho V_rel^3` held in
the band `[0.14, 0.32] W/cm2` (hard limit 0.54 W/cm2), with `V_rel` the atmosphere-relative speed and `rho`
the piecewise US76 fit of `examples/periapsis_corridor_config.toml` (valid 80-150 km; above the top node the
fit extends its last layer, and the band, about 134-147 km, lies inside the table).

Physics is the corridor example's (`examples/AGORA_Earth_Periapsis_Corridor.jl` is included, nothing is
copied); this study adds only `NBodyGravityModel` (Sun, Moon; SPICE), the heat-rate diagnostic, calibration to
a target peak Q, and plane angles. No GRAM, no `src/` change.

## Files

| File | Role |
|---|---|
| `study_config.toml` | Every named value, labeled SOURCED / DERIVED / ASSUMED |
| `common.jl` | Truth runs, calibration to peak Q, perturber-plane angles, process-route helper |
| `run_grid.jl` | Part B: passive 60-day grid, with and without Sun+Moon (`--pilot` runs one cell serially) |
| `select_geometries.py` | Picks the Part C geometries from the grid by a stated rule |
| `analyze_grid.py`, `plot_geometries.py` | Part B statistics, heat maps, geometry figure |
| `run_switch.jl` | Part C: `--phase=states\|kicks\|verify\|steps` (finite differences, extrapolation check, direct step) |
| `analyze_switch.py` | Part C/D tables and the switchability bar chart |
| `write_manifest.py` | Manifest: config, kernel hashes, git summary, versions, wall times |
| `results_root.py` | Shared helper: resolves the results root for the Python scripts |

## Commands (one Julia job at a time, memory capped)

```bash
cd <repo>
run() { systemd-run --user --scope -p MemoryMax=16G julia --project=. "$@"; }
run benchmarks/studies/periapsis_border_switchability/run_grid.jl --pilot          # one cell, serial
run benchmarks/studies/periapsis_border_switchability/run_grid.jl --workers=2      # 48 cells, process route
python3 benchmarks/studies/periapsis_border_switchability/analyze_grid.py
python3 benchmarks/studies/periapsis_border_switchability/select_geometries.py
for p in states kicks; do run benchmarks/studies/periapsis_border_switchability/run_switch.jl --phase=$p --workers=2; done
python3 benchmarks/studies/periapsis_border_switchability/analyze_switch.py        # writes switchability.csv
for p in verify steps; do run benchmarks/studies/periapsis_border_switchability/run_switch.jl --phase=$p --workers=2; done
python3 benchmarks/studies/periapsis_border_switchability/analyze_switch.py
python3 benchmarks/studies/periapsis_border_switchability/plot_geometries.py
python3 benchmarks/studies/periapsis_border_switchability/write_manifest.py
```

### Results root and reproduction

All study outputs go under one results root. Its value is `results_root` in `study_config.toml`
(default `output/periapsis_border_switchability`, relative to the repo root; `output/` is gitignored).
Override order:

- Julia drivers: env `SPACEAGORA_BORDER_RESULTS`, else `results_root`. `--out=DIR` overrides the output
  directory of one run (`run_static.jl` keeps its own default, `output/periapsis_static_region`).
- Python scripts: `--results-root DIR`, else env `SPACEAGORA_BORDER_RESULTS`, else `results_root`.
  A positional `OUT_DIR` argument, where a script takes one, overrides both.

The earlier border and border2 grids (the 50-degree grid and the 23.4/116.6/63.4-degree grids) were written to a
scratch directory outside the repository, and that data was left where it is. Those are the seeds of the
attribution drivers, which read `<results_root>/border` and `<results_root>/border2`. To reproduce or to
re-analyze from that data, point the root at the scratch directory:

```bash
SPACEAGORA_BORDER_RESULTS=<scratch-root> julia --project=. benchmarks/studies/periapsis_border_switchability/run_attribution.jl --workers=4 --out=$O
python3 benchmarks/studies/periapsis_border_switchability/analyze_attribution.py $O part1 --results-root <scratch-root>
```

### Which command writes which file

Run every command from the repo root, one Julia job at a time, memory capped. `O` is the output directory.

| Results file(s) | Command | Inputs |
|---|---|---|
| `grid/<cell>/passes_3b.csv`, `passes_ctl.csv`, `calibration.csv`, `grid_summary.csv`, `grid_mc.toml`, `pilot_row.csv` (pilot) | `run_grid.jl [--pilot] [--workers=N] [--out=O]` | `study_config.toml`, corridor example |
| `grid_summary_with_derived.csv`, `fig_heatmap_drp60.png`, `fig_hand_vs_measured.png` | `analyze_grid.py O` | `grid_summary.csv` |
| `geometries.csv`, `selection_scores.csv` | `select_geometries.py O` | `grid_summary.csv` |
| `fig_geometries_alt_q.png` | `plot_geometries.py O` | `grid/<cell>/passes_3b.csv`, `geometries.csv` |
| `kick_states.csv` | `run_switch.jl --phase=states --workers=2 --out=O` | `geometries.csv`, `grid/<cell>/calibration.csv` |
| `kicks/<tag>/<start>/<location>/`, `switch_fd.csv` | `run_switch.jl --phase=kicks --workers=2 --out=O` | `kick_states.csv` |
| `switch_derivatives.csv`, `switchability.csv`, `fig_switchability_bars.png` | `analyze_switch.py O` | `switch_fd.csv`, `geometries.csv` |
| `switch_verify.csv` | `run_switch.jl --phase=verify --workers=2 --out=O` | `kick_states.csv`, `switch_fd.csv` |
| `steps/<tag>/<start>/`, `direct_step.csv` | `run_switch.jl --phase=steps --workers=2 --out=O` | `geometries.csv` |
| `manifest.toml` | `write_manifest.py O` | every file above |

Static-region sweep (`run_static.jl`, `run_static_refband.jl`, `analyze_static.py`) writes `static_run_summary.csv`,
`control_run_summary.csv`, `race_run_summary.csv`, `*_mc.toml`, `static_cells.csv`, `control_cells.csv`,
`race_cells.csv`, `race_summary.csv`, `correction_cost.csv`, `ref50/<tag>/passes_3b.csv` and the `fig_*.png` files.
`analyze_static.py` reads the 50-degree reference from env `REF_DIR`, else `reference_dir` in `static_config.toml`
(a repo-relative path, not changed by the results-root setting).

Attribution (`run_attribution.jl`, `analyze_attribution.py`, `estimator_check.jl`, `analyze_estimator.py`) writes
`runs/<tag>/passes_<variant>.csv`, `attribution_run_summary.csv`, `attribution_mc.toml`, `attribution_table.csv`,
`attribution_contributions.csv`, `spectral_peaks.csv`, `quadrupole_amplitude_check.csv`, `fig_overlay_*.png`,
`fig_periodograms.png`, `estimator_series.csv`, `estimator_convention_check.toml`,
`estimator_quadrature_convergence.csv`, `estimator_stats.csv` and `fig_estimator_vs_truth.png`.

## Method summary

- **Run route.** Grid and kick campaigns go through `run_monte_carlo` on the process route (SPICE is not
  thread-safe). The router feature `density_family="gram_point"` is used only because it restricts the
  candidate routes to none/process (the corridor sweep's device); the run asserts `route == :process`.
- **Ephemerides.** The engine refuses N-body gravity under `SimpleEphemeridesModel`, so every run (with and
  without third bodies) uses `SpiceEphemeridesModel` (ITRF93 Earth orientation). The corridor example's
  analytic GMST frame is therefore replaced; geodetic altitude does not depend on longitude.
- **Vehicle.** Constant-beta drag of the corridor example. `beta = wet_mass / (Cd * A)`, `A` from the
  `examples/AGORA_Earth.jl` three-link reference areas (5.74 + 2 x 2.85 = 11.44 m2), `Cd = 2.2` ASSUMED
  (planner hand figure, not computed from the free-molecular tables), wet mass 820 kg (corridor config):
  beta = 32.58 kg/m2. With the 640 kg dry mass of `AGORA_Earth.jl` the same formula gives 25.4 kg/m2, the
  planner's hand figure.
- **Calibration.** The seed periapsis altitude is iterated until the first-pass peak Q equals
  `sqrt(0.14 * 0.32) = 0.2117 W/cm2` (tolerance 1e-3 relative, at most 6 iterations; 2-3 used). Initial
  periapsis 140.70 km geodetic, `V_rel` 10.26 km/s; Q at `V_rel` = 10.26 km/s gives h(0.54) = 127.96 km,
  h(0.32) = 134.65 km, h(0.14) = 146.74 km (band 12.09 km). Both runs of a cell start from the same state.
- **Plane angles.** `n = r_b x v_b` (geocentric J2000, SPICE `spkezr`), `k = n x h`, `omega' = atan2(h . (k x
  e_hat), k . e_hat)`, `i' = acos(h . n)`; a perturber in the equatorial plane gives the usual omega, i.
  The hand prediction is the double-averaged `dr_p/dt = -a (15/8) (mu_b / (n r_b^3)) e sqrt(1-e^2)
  sin(2 omega') sin^2(i')`, summed over Sun and Moon and integrated over the passes.
- **Finite differences.** Impulsive kick (inertial tangential = along velocity, radial, orbit-normal) on the
  unkicked state at three leg locations (apoapsis; inbound crossing of 2000 km geodetic altitude; outbound
  true anomaly 135 deg), at two starts (campaign start and the first periapsis at or after day 30), +/- each of
  0.05, 0.2, 1.0 m/s, 15 days of the full truth model. Drift-rate change = slope of (kicked - unkicked)
  periapsis altitude against time; slopes and step are least squares through the origin over the six probes.
  Both runs restart from the same extracted state.
- **Step compensation.** Every kick also steps the next periapsis. The compensated lever is the net rate change
  per total m/s when the kick's step is cancelled by a tangential apoapsis kick (superposition, linear
  response assumed). For the tangential-at-apoapsis kick itself the net is exactly zero.
- **switch_dv** = target rate change / compensated lever of the best (direction, location); linear
  extrapolation beyond 1 m/s is flagged (`extrapolation_factor`) and checked by the `verify` phase.
  `switch_dv_max` (ASSUMED proposal) = 0.5 x the measured direct 10 km step cost.

## Geometry selection rule (fixed before any kick result)

Crossing: the Sun+Moon periapsis history has an interior turning point between day 10 and 50 with both legs
over 2 km, ranked by in-band time, top 2. Static: not V-shaped, run completed, longest total in-band time,
top 2. Worst case: largest Sun+Moon effect among the remaining completed cells (scale only). The first
draft ranked the difference to the no-third-body run and picked window-edge extrema; it was replaced after
the grid was seen and before any kick was run.

## Known limits

- Passive runs whose periapsis plunged end early (engine `Terminated`, orbit decayed to impact); their
  "60 d" change is the value at termination (marked `*` in the heat map).
- Each exit time uses a constant baseline drift rate fitted over 15 days; the near-static geometries show a
  13-15 day periapsis oscillation (peak-to-peak up to ~10 km) that a constant rate does not capture.
- The Moon enters through SPICE ephemeris time from `SpiceEphemeridesModel`; the double-averaged hand
  formula uses instantaneous r x v planes of the Moon, so it is a check, not a model.

## Static-region sweep (extension)

Question: do critical-inclination (F1: 63.4 / 116.6 deg) or near-ecliptic (F2: 23.4 deg) orbits hold an extended static periapsis?
Files: `static_config.toml` (named values), `run_static.jl` (`--pilot`, `--phase=grid|control|race`), `run_static_refband.jl`
(heat-limited altitude columns for the previous 50-deg grid), `analyze_static.py OUT metrics|race|plots`; `write_manifest.py` appends the sweep config and route/workers.
`common.jl` is unchanged. Sequence (one Julia job at a time, memory capped, process route via `run_monte_carlo`):

```bash
run() { systemd-run --user --scope -p MemoryMax=32G julia --project=. "$@"; }
D=benchmarks/studies/periapsis_border_switchability; O=/path/to/out
run $D/run_static.jl --pilot --out=$O
run $D/run_static.jl --phase=grid --workers=4 --out=$O                 # 288 cells, 60 d
run $D/run_static_refband.jl /path/to/previous_50deg_out $O
REF_DIR=/path/to/previous_50deg_out python3 $D/analyze_static.py $O metrics                                 # static_cells.csv, control_cells.csv, race_cells.csv
run $D/run_static.jl --phase=control --workers=4 --out=$O
run $D/run_static.jl --phase=race --workers=4 --out=$O                  # 180 d + finite-difference lever
python3 $D/analyze_static.py $O race; python3 $D/analyze_static.py $O plots; python3 $D/write_manifest.py $O
```

Result: no family shows a static window beyond 26 d (50-deg reference: 30 d); the no-third-body control stays static for about 59 d, so the drift is the Sun+Moon term.
The static window is the longest run of passes with peak Q in the band and |alt - alt_0| <= X (X = 6 km, `static_config.toml`). The correction-cost estimate replays the passive
history with a deadband and the Keplerian tangential-at-apoapsis lever (checked by finite difference at the end of the 180-day runs). Details: `tasks/2026-10-15_periapsis-region-border-targeting/static_region_sweep_results.md`.

## Attribution and estimator check (extension)

Question: what drives the 41-110 km periapsis motion of the near-ecliptic (F2) cells and the 13-15 day oscillation, and does a single-averaged third-body model reproduce them?
Files: `attribution_config.toml` (cells, variants), `run_attribution.jl` (`--pilot`, `--workers=N`, `--out=DIR`; variants = Sun/Moon/neither, J3 off, drag off, Moon or Sun frozen at the start epoch, built by editing the effector tuple, no `src/` change),
`analyze_attribution.py OUT part1|part2|amp` (tables, overlays, periodograms, quadrupole amplitude check), `estimator_check.jl OUT` (offline single-averaged estimator: Gauss variational equations averaged over the Kepler ellipse, J2 secular rates, per-pass body positions) and `analyze_estimator.py OUT`; `write_manifest.py OUT` records the sweep config and route/workers.
The drivers read the calibrated seeds from `<results_root>/border` and `<results_root>/border2` (see Results root above).

```bash
run() { systemd-run --user --scope -p MemoryMax=32G julia --project=. "$@"; }
D=benchmarks/studies/periapsis_border_switchability; O=/path/to/out
run $D/run_attribution.jl --workers=4 --out=$O        # 98 runs, process route
python3 $D/analyze_attribution.py $O part1; python3 $D/analyze_attribution.py $O part2; python3 $D/analyze_attribution.py $O amp
run $D/estimator_check.jl $O 256; python3 $D/analyze_estimator.py $O; python3 $D/write_manifest.py $O
```

Result: the F2 null cells are driven by the Sun (its longitude moves only 59 deg in 60 d, so the double-averaged sin^2 i' formula does not apply), the 13-15 day oscillation is the Moon (half the sidereal month shifted by apsidal motion; amplitude matches the averaged quadrupole expression), J3 is below 2.4 km, drag matters only near decay, and the single-averaged model follows the drag-off truth to 0.25-1.4 km RMS in 12 cells. Details: `tasks/2026-10-15_periapsis-region-border-targeting/attribution_results.md`.

## Estimator v1 and selection skill (task A2)

Question: can a cheap per-pass mean-element estimator rank geometries by their 60-day maintenance dV and flag decays, offline, on the 336 truth cells?
Files: `estimator_v1.jl` (library: single-averaged exact Sun and Moon tide, secular J2 with the first-order J2 mean/osculating map, a pass-integrated drag term from the Gauss equations, band and decay), `run_selection_skill.jl` (driver; phases `selftest`, `cells`, `attrib`, `convergence`, `timing`, `scan`), `selection_skill_config.toml` (resolutions, thresholds, folds, scan grid), `analyze_selection_skill.py` (A1 policy on the predicted histories, baselines B1 double-averaged hand model, B2 quiet-set rule, B3 random, metrics, thresholds, figures).
The Julia half is serial and single process (SPICE is not thread-safe); the estimator has no fitted parameter. Inputs: `<results_root>/border`, `border2`, `border3` (truth) and `<results_root>/oracle/oracle_cell_labels.csv` (A1 labels).

```bash
run() { systemd-run --user --scope -p MemoryMax=16G julia --project=. "$@"; }
D=benchmarks/studies/periapsis_border_switchability; R=<results-root>     # holds border/, border2/, border3/, oracle/
for p in selftest attrib cells convergence timing scan; do run $D/run_selection_skill.jl --phase=$p --results-root=$R; done   # outputs under $R/selection
python3 -I $D/analyze_selection_skill.py --results-root $R                                                                       # tables and figures under $R/selection
```

| Results file(s) | Command | Inputs |
|---|---|---|
| `pred/<tag>.csv`, `pred_selfcal/<tag>.csv`, `hand/<tag>.csv`, `pred_cells.csv` | `run_selection_skill.jl --phase=cells` | truth `calibration.csv`, labels |
| `attrib_series.csv`, `j2_investigation_raw.csv` | `--phase=attrib` | `border3/runs/<tag>/passes_v{0,3,5}.csv` |
| `convergence_raw.csv`, `timing_candidates.csv`, `timing_summary.toml` | `--phase=convergence`, `--phase=timing` | as above |
| `scan_summary.csv`, `scan_hist.csv.gz`, `scan_timing.toml` | `--phase=scan` | `selection_skill_config.toml [scan]` |
| `selection_cell_table.csv`, `selection_metrics.csv`, `threshold_table.csv`, `confusion_*.csv`, `lookahead_error*.csv`, `attribution_cells_error.csv`, `dv_error_table.csv`, `decay_time_error.csv`, `safety_margin_cv.csv`, `fig_pred_vs_true_dv.png`, `fig_error_vs_lookahead.png`, `fig_confusion.png`, `fig_metrics_by_split.png`, `fig_seed_calibration.png` | `analyze_selection_skill.py --phase=cells` | the files above |
| `scan_candidates.csv`, `scan_top5.csv`, `scan_verifiable_top5.csv`, `fig_scan_heatmaps.png` | `analyze_selection_skill.py --phase=scan` | `scan_*` |

Result and limits: `tasks/2026-10-15_periapsis-region-border-targeting/selection_skill_results.md`. In short, on the 336 cells Spearman 0.997, top-10% recall 0.94-1.00, regret 0-0.065 m/s, all 133 decays predicted, but the false-safe threshold (0) fails (0.022-0.053 per fold; all cases near-boundary passive breaches); 3 ms per 60-day candidate.

## Blind verification of the estimator picks

Rule fixed in `blind_config.toml` before the runs (45 cells: per epoch the 5 predicted-cheapest, 5 random predicted-safe, 5 random predicted decay/breach, none already in the 336). `R` is the results root.

```bash
python3 -I $D/select_blind.py --results-root $R                    # writes $R/blind/blind_cells.csv (needs $R/selection/scan_candidates.csv)
run $D/run_blind.jl --pilot --results-root=$R --only=<tag>         # one cell, serial
run $D/run_blind.jl --results-root=$R --workers=4                  # process route, 45 cells
python3 -I $D/analyze_blind.py --results-root $R                   # labels (A1 pipeline), metrics, false-safe definitions, figures
python3 -I $D/part_manifest.py --results-root $R
```

## Superposition check

Rule in `superposition_config.toml` (8 completed cells: per inclination the 10th and 90th percentile of the A1 DP dV). Schedules come from the A1 library; Julia replays them as time-event tangential kicks (`--phase=replay`) or through the engine's thruster path (`--phase=finite`).

```bash
python3 -I $D/analyze_superposition.py select --results-root $R
run $D/run_superposition.jl --phase=null --results-root=$R         # zero-kick check + finite-difference lever probes
for i in 0 1 2; do
  python3 -I $D/analyze_superposition.py schedule --iter $i --results-root $R   # i=0 from the passive history, i>=1 repair increments
  run $D/run_superposition.jl --phase=replay --iter=$i --results-root=$R
done
run $D/run_superposition.jl --phase=finite --iter=0 --sched=lazy_s0 --only=<tag>,<tag> --results-root=$R   # at least 2 cells (a single sample is routed :none)
python3 -I $D/analyze_superposition.py analyze --results-root $R
```

## Estimator v2 (short-period terms predicted from the Pi metrics)

`estimator_v2.jl` adds to v1 the first-order Sun/Moon short-period (orbital-frequency) terms: the initial mean elements are the osculating ones minus the J2 and third-body short periods, and each pass's osculating periapsis state is mean plus both. `shortperiod_analytic.jl` holds the closed-form (Pi form) quadrupole and octupole expressions for the periapsis-radius rate, short period, and the oscillation amplitude, period and phase; `shortperiod_derivation.py` derives them with sympy and checks them against a direct integration. Config: `estimator_v2_config.toml`. Only `--results-root` carries a path; `R` as above.

```bash
python3 -I $D/shortperiod_derivation.py --results-root $R
for p in selftest cells attrib osc planets blind; do run $D/run_estimator_v2.jl --phase=$p --results-root=$R; done
run $D/run_estimator_v2.jl --phase=cells --variant=quadrupole --results-root=$R; run $D/run_estimator_v2.jl --phase=attrib --variant=quadrupole --results-root=$R
run $D/run_estimator_v2.jl --phase=timing --results-root=$R       # only with the machine otherwise idle
python3 -I $D/analyze_estimator_v2.py --results-root $R --phase all
```

Outputs under `$R/estimator_v2/` (v1-vs-v2 tables and figures prefixed `fig_`). Results and limits: `tasks/2026-10-15_periapsis-region-border-targeting/estimator_v2_results.md`.

## Stage E: baseline comparison

Question: does choosing the geometry with estimator v2 save propellant versus plain step-and-stay maintenance on a baseline geometry? `R` is the results root; named values in `baseline_closedloop_config.toml` (Part 2).

```bash
python3 -I $D/compare_baselines.py --results-root $R          # Part 1: no truth run; writes $R/compare_baselines/ (set_policy_stats.csv, ratios.csv, drift_floor*.csv, set_membership.csv, assumptions.txt, fig_*.png)
python3 -I $D/select_closedloop_cells.py --results-root $R    # Part 2: stratified 30-cell list, closedloop_cells.csv, written BEFORE any truth run
run $D/run_baseline_closedloop.jl --verify-tide=I0500_R270_W090_E3 --horizon-days=6 --results-root=$R   # predictor check, serial: no-tide / Arm A / Arm B next-periapsis shift vs truth
for arm in A B; do run $D/run_baseline_closedloop.jl --arm=$arm --smoke --results-root=$R; done   # 2 cells, 5 d each
for arm in A B; do run $D/run_baseline_closedloop.jl --arm=$arm --results-root=$R --workers=4; done   # full 60-day closed-loop runs, sequential, one Julia job, machine otherwise idle
python3 -I $D/analyze_closedloop.py --results-root $R          # per stratum and arm, vs open loop, failures
```

Part 1 uses `oracle/oracle_cell_labels.csv`, `estimator_v2/selection_cell_table.csv` and `blind/blind_cell_table.csv` under `R`. Baseline sets (B-ref50, B-grid, B-median, B-quiet) are defined by membership in the designed 336-cell sweep, which is an assumption; every dV is open loop (superposition, perfect forecast, unbounded burns, 60 d). Decayed cells report only the cost to termination, are kept in separate rows, and are never pooled with completed ones. Part 2 flies the corridor example's controller (`:near_passive`, its own predictor) in the engine with Sun and Moon, the corridor being the heat-rate band at the first-pass relative speed; US76 piecewise atmosphere, no GRAM. Results and limits: `tasks/2026-10-15_periapsis-region-border-targeting/baseline_comparison_results.md`.

Closed-loop arms (`run_baseline_closedloop.jl`, `analyze_closedloop.py`): both arms fly `select_corridor_maneuver` (`:near_passive`, same corridor, 1000 m margin, finite-thrust burn plans) in the truth engine with Sun and Moon. Arm A (`--arm=A`) predicts with the example's finite-thrust engine predictor rebuilt study-locally with SPICE ephemerides and `NBodyGravityModel` (the example's own predictor has no third bodies); Arm B (`--arm=B`) predicts with `estimator_v2.jl` re-initialized from the measured osculating state (impulsive dv). The example's `CorridorGuidancePolicy` hard-wires its predictor, so a study-local `PredictorPolicy` wrapper replaces only the decision hook (the example's callbacks and policy state are reused); no `src/` or `examples/` change. Two small study-local edits: `estimator_v1.jl` skips re-including `common.jl` when it is loaded, and `predict_history2` (`estimator_v2.jl`) takes an optional `init_override` (default unchanged). Results: `tasks/2026-10-15_periapsis-region-border-targeting/baseline_comparison_results.md`, section "Closed loop".

Follow-up closed-loop test (`select_followup_cells.py`, `analyze_followup.py`; config `[followup]`): new strata `quiet_not_top5`, `estimator_top_not_quiet`, `quiet_in_top` (disjoint from the first 30; written to `closedloop_cells_followup.csv` and, merged, `closedloop_cells_all.csv` before any run). All 24 quiet-set cells rank inside the estimator's top 10%, so the requested `quiet_not_top` is empty; `quiet_not_top5` (quiet cells outside the top 5%) is the named fallback. Arm C (`--arm=C`) uses the example's own predictor, unmodified (no Sun, no Moon).

```bash
python3 -I $D/select_followup_cells.py --results-root $R
run $D/run_baseline_closedloop.jl --arm=B --cells-file=closedloop_cells_followup.csv --results-root=$R --workers=4   # 18 new cells -> closedloop_B_followup/
run $D/run_baseline_closedloop.jl --arm=C --cells-file=closedloop_cells_all.csv --results-root=$R --workers=4        # 48 cells -> closedloop_C_all/
python3 -I $D/analyze_followup.py --results-root $R
```
Results: `tasks/2026-10-15_periapsis-region-border-targeting/baseline_comparison_results.md`, section "Follow-up closed loop".

## Stage C: planet generality (Mars)

Question: does the Earth-built estimator (`estimator_v2.jl`) predict periapsis histories at Mars without refitting, and do the Earth findings (geometry matters, the critical-inclination "quiet set") carry over? Plan and pre-registered thresholds: `tasks/2026-10-15_periapsis-region-border-targeting/mars_generality_plan.md` (written before any Mars truth run; addenda record every later addition). Results: `.../mars_generality_results.md`. The estimator is frozen: the only change to it is parametrization (planet record), proven bit-identical on Earth.

Files (all new unless marked):

| File | Role |
|---|---|
| `mars_generality_config.toml` | Every named Mars value, labeled SOURCED / DERIVED / ASSUMED |
| `mars_common.jl` | Mars truth runs (SPICE planet, Mars50c degree 3 order 0, Sun, single-layer exponential atmosphere, constant-beta drag), body-equator element frame, calibration to the band-center peak Q, periapsis table |
| `run_mars_generality.jl` | Driver, `--phase=cells\|pitable\|pi\|convergence\|pilot\|truth\|predict\|diag\|predict_diag\|predict_j3` |
| `analyze_mars_generality.py` | A1 labels on the Mars truth (`oracle_lib.py` unchanged, Mars GM set from the caller), metrics, pre-registered verdict table (`--variant j3` scores the unfrozen diagnostic) |
| `analyze_mars_diag.py`, `select_mars_diag_cells.py` | Which term fails: J3 exposure and the truth-variant (j2only / nosun / nodrag) comparison on 12 diagnostic cells |
| `estimator_v2_j3_diag.jl` | UNFROZEN DIAGNOSTIC: estimator v2 plus the orbit-averaged J3 term (derived, no fitted coefficient); never enters a verdict |
| `select_mars_closedloop_cells.py`, `run_mars_closedloop.jl`, `analyze_mars_closedloop.py` | Closed-loop arm B at Mars (example controller, Mars port of the predictor), cell rule, criteria, quiet-set shortcut test |
| `analyze_pi_collapse.py` | Pi tables and the Earth-plus-Mars collapse view |
| `c0_earth_regression.jl` | Earth regression test of the planet-record parametrization (frozen sources against the edited ones) |
| `estimator_v1.jl`, `common.jl`, `shortperiod_analytic.jl`, `run_estimator_v2.jl` (edited) | Parametrization only: `EstContext` gains `center` and `rot` (Earth: "EARTH", `nothing`), `PlanetRecord` and `earth_record`, `perturber_angles` and `bodystate` take `center` and `rot` keywords with the Earth defaults |

Commands (repo root, one Julia job at a time, memory capped; `R` is the results root, outputs under `$R/mars`):

```bash
run() { systemd-run --user --scope -p MemoryMax=16G julia --project=. "$@"; }
D=benchmarks/studies/periapsis_border_switchability; R=<results-root>
run $D/run_mars_generality.jl --phase=cells --results-root=$R           # cells.csv (128 cells), written before any run
run $D/run_mars_generality.jl --phase=pitable --results-root=$R         # constants, perturber distances, pole drift, Phobos/Deimos guard
run $D/run_mars_generality.jl --phase=convergence --results-root=$R     # numerical resolution of the estimator at Mars (no truth)
run $D/run_mars_generality.jl --phase=truth --workers=4 --results-root=$R   # process route, route asserted
run $D/run_mars_generality.jl --phase=predict --results-root=$R         # frozen estimator, all cells
python3 -I $D/analyze_mars_generality.py --results-root $R              # labels, metrics, verdict_table.csv, figures
python3 -I $D/select_mars_diag_cells.py --results-root $R; run $D/run_mars_generality.jl --phase=diag --workers=4 --results-root=$R
run $D/run_mars_generality.jl --phase=predict_diag --results-root=$R; run $D/run_mars_generality.jl --phase=predict_j3 --results-root=$R
python3 -I $D/analyze_mars_generality.py --results-root $R --variant j3; python3 -I $D/analyze_mars_diag.py --results-root $R
python3 -I $D/select_mars_closedloop_cells.py --results-root $R         # closedloop_cells.csv, written before any closed-loop run
run $D/run_mars_closedloop.jl --verify=<tag> --horizon-days=6 --results-root=$R      # predictor check
run $D/run_mars_closedloop.jl --workers=4 --results-root=$R             # arm B, frozen predictor (add --predictor=j3diag for the diagnostic arm)
run $D/run_mars_generality.jl --phase=pi --results-root=$R; python3 -I $D/analyze_pi_collapse.py --results-root $R
python3 -I $D/analyze_mars_closedloop.py --results-root $R
```

Earth regression test of the parametrization: copy `common.jl`, `estimator_v1.jl`, `estimator_v2.jl`, `shortperiod_analytic.jl`, `estimator_v2_config.toml`, `study_config.toml`, `attribution_config.toml` as they were before the edit into a directory at the same depth as this one (`common.jl` includes `../../../examples`), then `run $D/c0_earth_regression.jl --src=<frozen dir> --out=<out>/frozen --results-root=$R` and the same with `--src=$D --out=<out>/new`, and `cmp` the two trees. The test reads the Earth truth under `$R/border*` and `$R/oracle`.

Headline (details and limits in the results note): the frozen estimator passes the rank criteria pooled (Spearman 0.787, top-10% recall 0.846) and the pass-1 offset (22 m RMS), and fails the lookahead-error, false-safe, dV-error and classification criteria; the failure is the omitted J3 term (the frozen estimator matches a truth with J3 removed to 0.023 km median RMS, and the unfrozen J3 diagnostic matches the full truth to 0.01-0.09 km).

## Higher-harmonic parameter (Kaula)

Question: which gravity-harmonic terms must an estimator carry, and what does truncating the series cost? Pi_n = J_n (R/p)^n with an explicit orbit-averaged geometry factor Gamma_n(e, i, omega); Pi_K(N) is the worst-case periapsis-radius effect (km over the horizon) of all degrees above N, with the tail beyond the gravity file taken from Kaula's rule A_n ~ K/n^2 (K fitted per planet from the repo file, with its scatter). Derivation: `project_specs/kaula_harmonic_parameter.md`; pre-registered rules and pass criteria (written before any Mars cell was evaluated): `tasks/2026-10-15_periapsis-region-border-targeting/kaula_preregistration.md`; results: `.../kaula_parameter_results.md`. Python only (numpy, pandas, scipy, matplotlib; no Julia, no `src/` change).

| File | Role |
|---|---|
| `kaula_lib.py` | Gravity-file reader (finds `data/Gravity_harmonics_data` also from a worktree; `SPACEAGORA_GRAVITY_DIR` overrides), degree amplitudes, K fit, zonal secular geometry (`zonal_rates`, `gamma_max`) |
| `kaula_calibrate.py` | A_n, K and scatter per planet (Earth GGM05C, Mars50c, Venus MGNP180U), zonal roughness ratios, `fig_kaula_fit.png` |
| `kaula_cells.py` | Per-cell Pi_n, rho_n, D_n, Pi_K for the 336 Earth or 128 Mars cells (`--stage analytic`), and the conservative zonal Cowell check (`--stage cowell`, process pool) |
| `kaula_cowell.py` | Point mass + zonals propagator (DOP853, periapsis events); no drag, no third bodies |
| `kaula_selftest.py` | Textbook J2 rates, J3 and J4 closed forms, Kozai coefficient 15/8, analytic against Cowell J3 drift |
| `analyze_kaula.py` | Pre-registered criteria table, confusion tables, rank statistics, J4 sizing, figure |

Commands (repo root; `R` is the results root, outputs under `$R/kaula`; needs the Stage C Mars outputs and the Earth truth and v2 predictions under `R`):

```bash
D=benchmarks/studies/periapsis_border_switchability; R=<results-root>
python3 -I $D/kaula_selftest.py
python3 -I $D/kaula_calibrate.py --results-root $R
for pl in Earth Mars; do python3 -I $D/kaula_cells.py --planet $pl --stage analytic --results-root $R; done
for pl in Earth Mars; do python3 -I $D/kaula_cells.py --planet $pl --stage cowell --workers 4 --results-root $R; done   # about 2 min Earth, 1 min Mars; not a timing result
python3 -I $D/analyze_kaula.py --results-root $R
```

Headline (details and limits in the results note): measured K is 9.07e-6 (Earth), 1.02e-4 (Mars) and 9.90e-6 (Venus), scatter x1.30, x1.40 and x1.18; Mars J3 is on the rule (roughness 1.05), Mars needs J3 because K is 11 times larger and R/r_p = 0.969. The rule "carry degree n when rho_n > 4.63e-3" (threshold from Earth only) flags exactly the 64 Mars cells at 25.2 and 93.522 deg that needed J3 (one decay-endpoint cell missed, none of the 64 critical-inclination cells flagged), but also 32 borderline Earth cells; the along-history J3 prediction matches the Stage C frozen error (Spearman 0.92, median ratio 0.94). J4 alone is not needed at Earth and marginal at Mars; J5 is as large as J4 there, so the tail does not close with one more degree. Pre-registered criterion C8a (the worst-case envelope as a bound for degree 3) failed on the short-period offset at the critical inclinations. Tesseral and resonance effects are not covered.

## Estimator v3 and >=20x20 truth

Directive: every truth study uses gravity harmonics of at least degree 20 and order 20, and the estimator gains conditional higher-degree zonal terms driven by the Kaula parameter. Plan and pre-registered thresholds: `tasks/2026-10-15_periapsis-region-border-targeting/estimator_v3_plan.md` (addenda A to C record every later change and when); blind pre-registration: `estimator_v3_blind_preregistration.md`; results: `estimator_v3_results.md`.

Truth audit (read from the code): the Earth truth (336 cells) and the Stage C Mars truth (128 cells) used `GravitationalHarmonicsModel(3, 0, ...)` (central, J2, J3; hard-coded in the corridor example for Earth, `[gravity]` degree 3 order 0 for Mars). The new truth is degree 20 order 20 (EarthGGM05C file 99x99, Mars50c 50x50); Venus has a gravity file but no truth setup here and was not run.

| File | Role |
|---|---|
| `estimator_v3.jl` | Estimator v3: `predict_history3` = v2 plus the zonal degrees 3 to L chosen by the cumulative-tail form of the Kaula criterion (`ZonalOpts`: budget 1 km, `rule = :tail`, or `:per_degree`); exact zonal acceleration averaged over the Kepler ellipse by the same Gauss quadrature as the tide; no degree carried = v2 arithmetic (bit-identical) |
| `run_estimator_v3.jl` | Driver: `--phase=selftest` (averaged rates vs the disturbing-function rates, vs Vern9 Cowell integration, critical inclinations, node convergence), `mars`, `earth`, `mars_attrib` (gravity-term attribution), `timing` |
| `c0_v3_regression.jl` | Earth regression (the `c0_earth_regression.jl` pattern): v2 against v3 with an empty carried set (must be byte-identical) and against v3 by criterion |
| `run_earth_g20.jl` | Earth truth at degree/order from the command line (default 20x20); writes `study_config_g20.toml` (= `study_config.toml` + a `[gravity]` table, which `common.jl` `build_args` honors) |
| `mars_generality_config_g20.toml` | The Mars config at 20x20 (`mars_common.jl` takes degree and order from the config; new truth variants `zonal` = 20x0, `deg3` = 3x0 and `j2only`) |
| `earth_g20_prepare.py`, `earth_blind_labels.py` | Lay the new Earth truth out as the old analysis scripts expect (`border`, `border2`, `static_cells.csv`); A1 labels for the blind Earth cells |
| `analyze_earth_g20.py` | Earth scoring of any prediction directory with the A2 machinery, and the per-cell old-vs-new truth difference |
| `analyze_attrib_gravity.py` | Mars gravity-term attribution: truth against truth (tesseral, zonal 4-20, J3) and predictor against truth variant |
| `select_blind_g20.py`, `blind_g20_config.toml` | The blind cells (64 Mars, 40 Earth; seeds and rules in the config, drawn before any blind run) |
| `analyze_closedloop_g20.py` | Closed-loop check summary (Arm B with the v3 and frozen v2 predictors) |
| `common.jl`, `mars_common.jl`, `run_mars_generality.jl`, `run_mars_closedloop.jl`, `analyze_mars_generality.py`, `select_mars_closedloop_cells.py` (edited) | Optional `[gravity]` override and `--config`; truth variants; `--variant v3`; `--predictor=v3`; folds for off-grid RAAN |

Commands (repo root; `R` is a NEW results root for the 20x20 truth, `OLD` the root of the 3x0 truth; one Julia job at a time, memory capped):

```bash
run() { systemd-run --user --scope -p MemoryMax=16G julia --project=. "$@"; }
D=benchmarks/studies/periapsis_border_switchability; G=$PWD/$D/mars_generality_config_g20.toml
run $D/run_estimator_v3.jl --phase=selftest --results-root=$R
# truth at >= 20x20 (process route via run_monte_carlo, route asserted)
run $D/run_mars_generality.jl --phase=cells --config=$G --results-root=$R; run $D/run_mars_generality.jl --phase=truth --workers=4 --config=$G --results-root=$R
run $D/run_earth_g20.jl --phase=cells --results-root=$R; run $D/run_earth_g20.jl --phase=grid --workers=4 --results-root=$R
python3 -I $D/earth_g20_prepare.py --results-root $R --old-root $OLD; python3 -I $D/oracle_run.py --data-root $R --out-dir $R/oracle
python3 -I $D/analyze_earth_g20.py --results-root $R --mode truthdiff --old-root $OLD
# baseline estimators on the new truth (copy mars/estimator_numerics.toml and planet_constants.csv from the old Mars root first)
run $D/run_mars_generality.jl --phase=predict --config=$G --results-root=$R; run $D/run_mars_generality.jl --phase=predict_j3 --config=$G --results-root=$R
run $D/run_selection_skill.jl --phase=cells --results-root=$R; run $D/run_estimator_v2.jl --phase=cells --results-root=$R
# attribution: truth variants on all 128 cells, then the predictors against each
cp $R/mars/cells.csv $R/mars/diag_cells.csv; run $D/run_mars_generality.jl --phase=diag --variants=zonal,deg3,j2only --workers=4 --config=$G --results-root=$R
# estimator v3 (development sets), scoring
run $D/run_estimator_v3.jl --phase=mars --config=$G --results-root=$R; run $D/run_estimator_v3.jl --phase=earth --results-root=$R
run $D/run_estimator_v3.jl --phase=mars_attrib --config=$G --results-root=$R; run $D/c0_v3_regression.jl --out=$R/c0_v3 --results-root=$R
python3 -I $D/analyze_mars_generality.py --results-root $R [--variant j3|v3|v3a]; python3 -I $D/analyze_attrib_gravity.py --results-root $R
for d in selection estimator_v2 estimator_v3; do python3 -I $D/analyze_earth_g20.py --results-root $R --mode score --pred-dir $R/$d; done
# blind sets (root $B): select_blind_g20.py writes the cells and the epoch-shifted Mars config; then the same runs with --config=$B/mars_blind_config.toml, run_earth_g20.jl --cells= --epochs=
python3 -I $D/select_blind_g20.py --results-root $B
```

Headline (details, tables and limits in the results note). The Earth truth hardly changes at 20x20 (same status in all 336 cells, A1 labels move by at most 0.049 m/s, rank correlation 1.0000), so the Earth selection-skill and cost findings stand; the Mars truth changes a great deal (34 of 128 cells change completion status, median altitude-history change 3.1 km). The averaged zonal theory is right (v3 with every zonal degree matches a 20x0 truth to 0.03 km median, selftests agree with Cowell integration to 1.1% of the effect), but at Mars the omitted non-zonal terms dominate: v3 matches the 20x20 truth to about 2.4 to 3 km median RMS, the tesseral/sectorial share of the truth change is about 3 km against 0.4 to 0.8 km for zonal degrees 4 to 20. On the blind Mars set v3 fails the Stage C criteria (as pre-registered) but is better than v2 and the J3 diagnostic; at Earth v3 is not worse than v2 and reduces the 50-60 day altitude error from 0.77 to 0.53 km on the blind cells.

## Tesseral and sectorial terms (derivation)

Question: why do the non-zonal terms move the Mars periapsis by kilometers, what is the closed form of their per-pass effect, which dimensionless parameter says when to carry them, and can an estimator carry them? Derivation and results: `tasks/2026-10-15_periapsis-region-border-targeting/tesseral_derivation.md`; pre-registered parameter test (written before any per-cell value): `.../tesseral_preregistration.md`. Python only (numpy, scipy, pandas; no Julia, no `src/` change).

| File | Role |
|---|---|
| `tesseral_lib.py` | Gravity coefficients (unnormalized), vectorized Cunningham acceleration and potential, Kepler arc, Gauss rates, per-pass kick by Gauss quadrature (`gauss_kick`), per-term unit kicks for all (n, m) at once (`term_kicks`), complex inclination function, Hansen coefficients by FFT, Kaula's F_nmp, the Lagrange-series per-pass kick with the commensurability filter (`series_kick`) |
| `tesseral_selftest.py` | Acceleration vs potential gradient, zonal part vs `kaula_cowell`, Kaula F vs Fourier definition, series reconstruction of the potential, series kick vs quadrature, per-term kicks vs the full field |
| `tesseral_cowell.py` | Vectorized Cowell (all variants in one state vector): `--set mech` (10 orbits x 41 field variants, no drag, 40 passes), `--set drag` (6 orbits, 20x0 and 20x20, Stage C exponential atmosphere, 60 d) |
| `tesseral_analyze.py` | `--set mech` (map, slow-term and truncated-series predictions against Cowell), `--set loc` (where along the arc the r_p kick is collected), `--set drag` (open-loop map along the reference) |
| `tesseral_param.py` | Pi_T, D_T, D_T_hist for the 128 Mars and 336 Earth cells, the deterministic map along the Mars truth with the IAU rotation angle, the pre-registered criteria; `--tail` the degree 21-50 and truncation statement |
| `tesseral_propagator.py` | Minimal per-pass stepper (drag + field by quadrature, timing from the mean-anomaly equation) with and without the tesseral field against the drag Cowell set |

Commands (repo root; `R` the 20x20 results root, `OLD` the 3x0 root for the Earth `pi_cells.csv`):

```bash
D=benchmarks/studies/periapsis_border_switchability
python3 -I $D/tesseral_selftest.py
python3 -I $D/tesseral_cowell.py --set mech --results-root $R; python3 -I $D/tesseral_cowell.py --set drag --results-root $R   # about 2 and 4 min, single process
python3 -I $D/tesseral_analyze.py --set mech --results-root $R    # about 15 min (series truncation study)
python3 -I $D/tesseral_analyze.py --set loc --results-root $R; python3 -I $D/tesseral_analyze.py --set drag --results-root $R
python3 -I $D/tesseral_param.py --results-root $R --pi-cells $OLD/mars/pi_cells.csv; python3 -I $D/tesseral_param.py --results-root $R --tail
python3 -I $D/tesseral_propagator.py --results-root $R
```

Headline (details and limits in the note): the all-tesseral per-pass periapsis-radius kick at the Odyssey-like Mars orbit is 210 to 350 m (Earth cells 25 m); it survives averaging because m nu (nu = w_rot T / 2 pi = 0.69) is close to an integer for the orders with power (m = 3, 10), and drag sweeps nu to about 0.14 through the commensurabilities 2/3, 1/2, 1/3, where each crossing adds about 10 per-pass kicks. Degrees 2 and 3 collect their kick over the whole orbit, not at periapsis. The Kaula-series per-pass kick (with the filter 2 sin(pi(j - m nu))/(j - m nu)) equals the quadrature to 2e-9 and the first-order map matches Cowell to 2.7% median. Pi_T = sigma_n (R/r_p)^n with the accumulation gain gives D_T: pre-registered test 6 of 7 pass (magnitude 1.34, Earth all below 1 km, contrast 34; cell ranking at Mars fails, Spearman 0.26). A per-pass stepper that carries the 20x20 kick tracks the drag Cowell truth at its own zonal floor (0.02 to 1.1 km per bin) against 0.5 to 6 km without it; as an open-loop correction along a fixed reference it holds only about 10 days.

## Numerical predictor at Mars (option 3)

Question: if a short numerical propagation with the full gravity field replaces the closed-form estimator at Mars, what does it buy and what does it cost? Plan and pass/fail thresholds (written before the runs): `tasks/2026-10-15_periapsis-region-border-targeting/numerical_predictor_plan.md`; results: `.../numerical_predictor_results.md`.

N(field) is the Mars truth engine itself (`mars_common.jl`: SPICE Mars, Sun, exponential atmosphere with constant-beta drag, truth tolerances), passive, restarted from the initial or measured osculating state, with Mars50c harmonics of a chosen degree and order; no engine code is re-implemented and nothing under `src/` or `examples/` changes. It writes the estimator-v3 output schema (`pred_<name>/`, `pred_<name>_selfcal/`, `pred_cells_<name>.csv`), so `analyze_mars_generality.py --variant <name>` scores it unchanged (its only edit: the variant prefix `num` is accepted). Names: `num<D>x<O>` (truth tolerances T1), `...L` (ASSUMED looser tolerances T2: reltol 1e-6, abstol 1e-8, dt_max 120 s / 5 s), `...T` (tighter T0: reltol 1e-11, abstol 1e-13, integrator-error reference), `numhyb<F>_<base>` (hybrid: v3 for all cells, N for the F-percent cells with the lowest v3 predicted dV).

| File | Role |
|---|---|
| `run_numerical_predictor.jl` | Driver, `--phase=probe\|ladder\|timing\|verify\|closedloop`; includes `run_mars_closedloop.jl` (policy wrapper, corridor config); the closed-loop predictor is the Arm A finite-thrust rollout with the Mars setup and the numerical field |
| `analyze_numerical_predictor.py` | `--phase=hybrid` (builds the mixed v3/N prediction folders), `summary` (consistency, pass counts against v3, per-bin tables, cost table), `closedloop` (per-stratum table beside frozen v2 and v3) |

```bash
run() { systemd-run --user --scope -p MemoryMax=16G julia --project=. "$@"; }
D=benchmarks/studies/periapsis_border_switchability; R=<blind 20x20 results root>; C=$R/mars_blind_config.toml   # or the 128-cell root with mars_generality_config_g20.toml
run $D/run_numerical_predictor.jl --phase=probe --fields=20x20,4x4,2x0 --tols=T1,T2 --only=<tag> --results-root=$R --config=$C       # scale of the campaign (not a table number)
run $D/run_numerical_predictor.jl --phase=ladder --tol=T1 --workers=4 --results-root=$R --config=$C            # 8 fields x 64 cells, process route, route asserted
run $D/run_numerical_predictor.jl --phase=ladder --fields=16x16 --tol=T1 --workers=4 --results-root=$R --config=$C
run $D/run_numerical_predictor.jl --phase=ladder --fields=2x0,4x4,8x8,12x12,16x16,20x20 --tol=T2 --workers=4 --results-root=$R --config=$C
for f in 2x0 3x0 8x0 20x0 4x4 8x8 12x12 16x16 20x20; do python3 -I $D/analyze_mars_generality.py --results-root $R --variant num$f; done   # and num<f>L for the T2 fields
python3 -I $D/analyze_numerical_predictor.py --results-root $R --phase hybrid --base num20x20L --fractions 0.10,0.20
for f in 10 20; do python3 -I $D/analyze_mars_generality.py --results-root $R --variant numhyb${f}_num20x20L; done
run $D/run_numerical_predictor.jl --phase=timing --fields=2x0,3x0,8x0,20x0,4x4,8x8,12x12,16x16,20x20 --tols=T1,T2 --ncells=8 --reps=2 --results-root=$R --config=$C   # ONE job, machine otherwise idle
run $D/run_numerical_predictor.jl --phase=verify --tag=<tag> --field=20x20 --tol=T1 --horizon-days=6 --results-root=$R --config=$C                  # restart consistency from the measured state
run $D/run_numerical_predictor.jl --phase=ladder --fields=20x20 --tol=T0 --only=<8 tags> --workers=4 --results-root=$R --config=$C                    # integrator error of the truth settings
for v in "T1 fixed" "T1 perpass" "T2 fixed"; do set -- $v; run $D/run_numerical_predictor.jl --phase=closedloop --field=20x20 --tol=$1 --band=$2 --workers=4 --results-root=$R --config=$C; done
python3 -I $D/analyze_numerical_predictor.py --results-root $R --phase summary; python3 -I $D/analyze_numerical_predictor.py --results-root $R --phase closedloop
```

Headline (details, tables and limits in the results note): the full 20x20 field removes the Stage C failures only when it is the full field (measured against a truth that shares the same gravity file, atmosphere, drag and Sun, so this is a perfect-model upper bound). With truth tolerances the engine reproduces the truth to 0 m and restarts from a measured state to 7e-5 m; no truncated field meets the 1 km criterion at every lookahead (16x16: first four bins only); at 20x20 and the looser tolerance 55 of 56 blind sub-criteria pass (v3: 17) at about 2.1 s per 60-day candidate (8.9 s at truth tolerances), 230 to 1000 times the v3 cost and not field dependent. Scanning 6,336 candidates takes 3.6 core-hours (T2) or 15.6 (T1); a hybrid refining the top 20% fits one core-hour only at T2 and recovers 15 of the 39 v3 failures.

## Estimator v4 (tesseral)

Directive: build and test v3 plus the per-pass tesseral/sectorial kick of the derivation. Plan, pre-registered thresholds and addenda (written before the runs they govern): `tasks/2026-10-15_periapsis-region-border-targeting/estimator_v4_plan.md`; results: `.../estimator_v4_results.md`.

| File | Role |
|---|---|
| `estimator_v4.jl` | `tesseral_pass_kick` (Gauss quadrature of the non-zonal field along the Kepler arc, inbound/outbound halves, timing terms), `RotationModel` (engine frame routine), carry rule D_T (Pi_T in Julia), `predict_history4`, `calibrate_seed4`; with the rule false it returns the v3 prediction unchanged |
| `run_estimator_v4.jl` | Driver: `--phase=selftest\|mars\|earth\|timing`, `--lt=N` (reduced field), `--name=`, `--oracle-time` (UNFROZEN attribution diagnostic) |
| `v4_reference_kicks.py` | Python reference values for selftests a and b |
| `c0_v4_regression.jl` | Earth regression (v3 against v4, and a forced-kick diagnostic) |
| `select_blind_v4.py`, `blind_v4_config.toml` | The second blind Mars set (64 new cells, disjoint from the development and v3-blind sets) |
| `analyze_estimator_v4.py` | Stage C pass counts, gains and losses against v3, T5b, per-cell RMS, pass-time error at commensurability crossings, beside v3 and the numerical predictor |
| `analyze_mars_generality.py`, `run_mars_closedloop.jl` (edited) | variant prefix `v4`; `--predictor=v4` (kick forced on) |

```bash
run() { systemd-run --user --scope -p MemoryMax=16G julia --project=. "$@"; }
D=benchmarks/studies/periapsis_border_switchability; G=$PWD/$D/mars_generality_config_g20.toml   # R: 20x20 truth root
python3 -I $D/v4_reference_kicks.py --results-root $R
run $D/run_estimator_v4.jl --phase=selftest --which=a,c,d,f,n --config=$G --results-root=$R
run $D/run_estimator_v4.jl --phase=mars --config=$G --results-root=$R; for L in 6 8 12; do run $D/run_estimator_v4.jl --phase=mars --lt=$L --name=v4_L$L --config=$G --results-root=$R; done
python3 -I $D/analyze_mars_generality.py --results-root $R --variant v4; python3 -I $D/analyze_estimator_v4.py --results-root $R --variants v3,v4,num20x20L --tag dev
run $D/run_estimator_v4.jl --phase=earth --results-root=$R [--study=$R/study_config_g20.toml for the blind Earth root]   # then diff -rq against estimator_v3/{pred,pred_selfcal}
python3 -I $D/select_blind_v4.py --results-root $B2 --exclude-root $B      # then truth: run_mars_generality.jl --phase=truth --workers=4 --config=$B2/mars_blind_config.toml
run $D/run_estimator_v4.jl --phase=timing --config=$G --results-root=$R      # ONE job, machine otherwise idle
run $D/run_mars_closedloop.jl --predictor=v4 --workers=4 --config=$B/mars_blind_config.toml --results-root=$B
```

Headline (details and limits in the results note). The kick is correct (Python reference 1.5e-9, Cowell 0.1%), removes the pass-1 offset (0.30 to 0.02 km) and cuts the 0-20 d altitude error to 0.09-0.85 km, at 25 ms per 60-day prediction (scan of 6,336 candidates: 0.045 core-hours); Earth is byte-identical to v3 (376 cells). The frozen v4 does not meet the pre-registered criteria on the new blind set (33 of 56 Stage C sub-criteria against v3 21 and the numerical predictor 55; T5b fails beyond 20 d): the stepper's own pass times drift (the drag-driven change of a is 1% off), which sets the wrong rotation phase for later kicks. With the truth pass times fed to the kick only (unfrozen diagnostic) v4 reaches 55 of 56. Two timing corrections found along the way (outbound leg on the post-drag ellipse; base period without the transient half kick) apply to the carried path only.

## Estimator v5 (pass timing)

Directive: find and fix the pass-timing error of the closed-form periapsis estimator under drag-driven orbit decay (the stepper shared by v1 to v4). Plan, pre-registered thresholds and the addendum (written before the blind runs): `tasks/2026-10-15_periapsis-region-border-targeting/estimator_v5_plan.md`; results and the error budget: `.../estimator_v5_results.md`.

| File | Role |
|---|---|
| `estimator_v5.jl` | `TimingOpts`, `drag_window` (numerical integration of the real perturbed drag window, invariants accumulated along the path), `predict_history5` (v4 plus the three timing corrections; works with and without the tesseral kick, so it also corrects Earth), `calibrate_seed5`. With `window = :kepler`, `post_drag_period = false`, `kick_post_drag = false`, `zonal_sp = false` it is v3 (rule false) or v4 (rule true) bit for bit |
| `run_estimator_v5.jl` | Driver: `--phase=selftest\|mars\|earth\|timing`, `--name=`, `--window=`, `--nwindow=`, `--kickpd=`, `--zsp=`, `--oracle-time=1` (UNFROZEN attribution diagnostic) |
| `diag_v5_lib.jl` | Cowell propagator of one orbit leg with switchable forces (point mass, J2, zonal, non-zonal 20x20, Sun and Moon, drag): the reference of the diagnostics (reproduces the Mars engine truth to 2 ms and 0.3 m per orbit) |
| `diag_v5_timing.jl`, `diag_v5_earth.jl` | Pass-timing error budget by truth-state injection (`E_half`, `E_full`, radius errors) and by 40-pass (Earth 30) chains, per force subset, Mars and Earth |
| `diag_v5_window_check.jl` | R2: the drag increment of one pass against the Cowell reference, Kepler quadrature against the integrated window at several resolutions |
| `c0_v5_regression.jl` | Earth regression (24 cells): legacy options against v3, no-drag identity, drag runs |
| `analyze_estimator_v5.py`, `analyze_v5_criteria.py`, `analyze_v5_budget.py`, `analyze_earth_v5.py`, `analyze_v5_closedloop.py` | Stage C counts and pass-time tables beside v3/v4/N; pass/fail of the plan's criteria; budget tables; Earth altitude and pass-time error by lookahead bin; closed-loop table |
| `analyze_mars_generality.py`, `run_mars_closedloop.jl` (edited) | variant prefix `v5`; `--predictor=v5` (kick forced on) |

```bash
run() { systemd-run --user --scope -p MemoryMax=16G julia --project=. "$@"; }
D=benchmarks/studies/periapsis_border_switchability; G=$PWD/$D/mars_generality_config_g20.toml   # R: 20x20 truth root (development), B2/B: blind roots
run $D/run_estimator_v5.jl --phase=selftest --config=$G --results-root=$R                      # identities R1
run $D/run_estimator_v5.jl --phase=mars --config=$G --results-root=$R; run $D/run_estimator_v5.jl --phase=earth --results-root=$R
python3 -I $D/analyze_mars_generality.py --results-root $R --variant v5; python3 -I $D/analyze_estimator_v5.py --results-root $R --variants v3,v4,v5,num20x20L --tag dev
python3 -I $D/analyze_v5_criteria.py --results-root $R --tag dev --set dev; python3 -I $D/analyze_earth_v5.py --results-root $R --variants estimator_v4,estimator_v5
run $D/run_estimator_v5.jl --phase=mars --config=$B2/mars_blind_config.toml --results-root=$B2     # blind-2 (once, frozen file); same for $B (old blind, flagged seen)
run $D/run_estimator_v5.jl --phase=earth --study=$B/study_config_g20.toml --results-root=$B        # Earth blind g20
run $D/diag_v5_timing.jl --root=$R --cells=<tags> --variants=full,nodrag,j2drag,drag2b --stepper=v4|v5 [--mode=chain --npass=40]   # budget; diag_v5_earth.jl for Earth
run $D/c0_v5_regression.jl --out=<dir> --results-root=$R; run $D/diag_v5_window_check.jl --root=$R
run $D/run_estimator_v5.jl --phase=timing --config=$G --results-root=$R                          # ONE job, machine otherwise idle
run $D/run_mars_closedloop.jl --predictor=v5 --workers=4 --config=$B/mars_blind_config.toml --results-root=$B
```

Headline (details and limits in the results note). The pass-time error is the accumulated error of the drag-driven change of a, not of the period formula (exact to 0.01 s per orbit without drag). Three corrections, each found by measurement: the drag increments are integrated along the real perturbed path instead of the Kepler ellipse (the first-order quadrature overestimated the change of a by 0.84% at Mars, which a fresh period turns into 6 to 8 s per orbit), the outbound half of the tesseral kick is evaluated on the post-drag ellipse (27 m of periapsis radius per orbit at Mars), and an uncarried pass takes half the post-drag period (Earth: 72 to 180 s per orbit, half of the drag-induced period change was missing). Pooled pass-time RMS on the Mars development cells falls 238 / 1,184 / 4,925 / 21,512 / 71,185 / 165,322 / 294,091 s (v4) to 33 / 231 / 1,004 / 3,317 / 10,596 / 33,631 / 84,016 s by 5-day lookahead bin 0-5 .. 50-60 d; on blind-2 105 / 830 / 4,677 s (v4) to 23 / 150 / 965 s. Earth completed cells: 207 / 699 / 1,267 / 2,104 / 3,236 / 5,014 / 8,290 s to 12 / 53 / 181 / 558 / 1,315 / 2,649 / 5,886 s on the 336 development cells. The Stage C selection score does not reach the pre-registered levels: 28 of 46 (development; v4 18), 33 of 56 (blind-2; v4 33, N(20x20) 55), 36 of 56 (old blind; v4 30). The rotation phase of the kicks is still set by pass times that are still 3,300 s (development) to 6,400 s (blind-2) off at 20 to 30 d; fed the truth pass times the v5 kick reaches 55 of 56. What now limits the timing is the periapsis radius (5 to 15 m level errors from J2 x zonal, Sun and tesseral cross terms of the osculating-to-mean map; at Earth 80 to 140 m per orbit from the Moon) amplified by the density (Mars scale height 6.31 km: 0.16% per 10 m).

## Estimator status (closed 2026-10-19)

The closed-form estimator line (v1 to v5) is closed. No further estimator version is planned. Numbers below are from the local notes of each stage (Stage C, v3, v4, v5, numerical predictor); nothing here is a new measurement.

**Earth: validated.** v2 against the Earth truth: pooled RMS altitude error 0.12 / 0.13 / 0.18 / 0.31 / 0.42 / 0.52 / 0.68 km at 0-5 / 5-10 / 10-20 / 20-30 / 30-40 / 40-50 / 50-60 d, pass-1 offset 135 m (Stage C note, section 4.2). v5 improves it on the 172 completed development cells: pooled altitude RMS 0.119 / 0.106 / 0.156 / 0.222 / 0.294 / 0.367 / 0.466 km and pooled pass-time RMS 12 / 53 / 181 / 558 / 1,315 / 2,649 / 5,886 s (v5 note, section 7).

**Mars: validated for ranking and for prediction up to about 10 days; open-loop 60-day prediction is NOT validated.** Ranking inside the critical-inclination set: Spearman 0.989 over the 64 critical cells (Stage C note, section 6). Short lookahead, v5 on blind-2: pooled pass-time RMS 23 / 150 / 965 s at 0-5 / 5-10 / 10-20 d and pooled altitude RMS 0.020 / 0.033 / 0.232 km (v5 note, sections 1.3 and 5). The pre-registered Stage C criteria are not met at 60 days. Stage C sub-criteria passed of 56 on blind-2: v3 21, v4 33, v5 33, numerical predictor N(20x20, looser tolerances) 55; v5 with the truth pass times fed to the kick (unfrozen diagnostic) 55 of 56 (v5 note, section 5; numerical predictor note, section 3).

**Cause.** The score is limited by pass timing alone, and the timing is limited by the periapsis radius: level errors of 5 to 15 m (second-order J2 x zonal, Sun and tesseral cross terms of the osculating-to-mean map) are amplified by the density (0.16% per 10 m of periapsis radius at Mars) into the drag-driven change of the semi-major axis, and the pass-time error is the double integral of that (v5 note, section 8). The remedies are a second-order J2 osculating-to-mean map or numerical refinement of the periapsis state. Both were deliberately not pursued.

**Use guidance.** Where Mars needs prediction deeper than about 10 days, use the 20x20 numerical predictor (`run_numerical_predictor.jl`): 7.2-8.9 s per 60-day candidate at truth tolerances and 1.6-2.1 s at the looser tolerances, 230 to 1000 times the closed-form v3 cost; it is a refiner, not a scan tool (numerical predictor note, sections 1.4 and 4). It is a perfect-model upper bound: it shares the truth's gravity file, atmosphere, drag and Sun. The closed-loop controller does not need long open-loop prediction, because it re-anchors every apoapsis from the measured state; closed-loop results are identical for v4 and v5 (a one-pass predictor cannot see the timing) and the predictor choice did not change them materially (v5 note, section 10).

**Why it is closed.** The project goal is the Pi parameters as a planet-general feature set for a multi-planet meta-learning extension of an aerobraking reinforcement-learning policy. What that needs are one-pass-ahead features and a fast surrogate, not 60-day open-loop accuracy, and the numerical predictor covers the deeper cases.

## Planet generality at Venus (development set)

Question: do the frozen estimators v3 and v5, the apoapsis-decay check and the Pi bookkeeping carry to a third planet, and what are the Venus values of the Pi groups? Plan and pre-registered expectations (written before any Venus truth run): `tasks/2026-10-15_periapsis-region-border-targeting/venus_generality_plan.md`; results: `.../venus_generality_results.md`. New files only; the estimators and the Mars truth machinery are reused unchanged (`venus_common.jl` swaps the planet).

| File | Role |
|---|---|
| `venus_generality_config.toml` | Every named Venus value, labeled SOURCED / DERIVED / ASSUMED (atmosphere: single exponential derived from the Venus Express table; MGNP180U 20x20; VEX vehicle; 32-cell grid) |
| `venus_common.jl` | Venus setup (`MarsSetup` with the Venus planet and derived atmosphere), the cell task |
| `run_venus_generality.jl` | Driver, `--phase=cells\|pitable\|convergence\|convergence_refined\|pilot\|truth\|predict` (process route, route asserted) |
| `venus_analyze.py`, `venus_apoapsis_pi_check.py`, `venus_pi_cells.py` | Estimator errors by lookahead bin and the V1-V6 verdicts; the D1V/D2V apoapsis checks; the per-cell Pi table, group ranges and the hand check |

Commands (repo root, one Julia job at a time, memory capped; `R` is the results root, outputs under `$R/venus`):

```bash
D=benchmarks/studies/periapsis_border_switchability; R=<results root>; run="systemd-run --user --scope -p MemoryMax=16G julia --project=."
for ph in cells pitable convergence convergence_refined; do $run $D/run_venus_generality.jl --phase=$ph --results-root=$R; done
$run $D/run_venus_generality.jl --phase=truth --workers=4 --results-root=$R     # 32 cells, process route
$run $D/run_venus_generality.jl --phase=predict --results-root=$R               # frozen v3 and v5
python3 -I $D/venus_analyze.py --results-root $R; python3 -I $D/venus_apoapsis_pi_check.py --results-root $R; python3 -I $D/venus_pi_cells.py --results-root $R
```

## Venus with GRAM

The Venus generality test again, with the native Venus-GRAM MEAN density as truth (the exponential section above used a single layer derived from the Venus Express table). Same 32-cell grid, vehicle, 20x20 MGNP180U gravity, Sun third body and thresholds; the frozen estimators v3 and v5 are unchanged. Plan and pre-registered expectations (written before any GRAM truth run): `tasks/2026-10-15_periapsis-region-border-targeting/venus_gram_plan.md`; results: `.../venus_gram_results.md`.

Needs the native GRAM build and the Venus-GRAM and SPICE data of `data/GRAMSuite.jl` (see `docs/src/user/gramsuite_setup.md`). Truth drag is the corridor example's `CorridorBallisticDrag(..., :sampled)` on the engine's GRAM density with EI = 2000 km (below EI the GRAM extension is used, above it a planet polynomial), direct sampling, no density perturbation (`SPACEAGORA_GRAM_DENSITY_PERTURBATION` unset). The truth sweep goes through `run_monte_carlo` with `density_family="gram"`; the route is asserted to be `:process` (native GRAM is not thread-safe).

| File | Role |
|---|---|
| `venus_gram_config.toml` | Overrides only (GRAM settings, Pi-rule and fit windows, calibration iterations); every other table is inherited from `venus_generality_config.toml`. Labels SOURCED / DERIVED / ASSUMED |
| `venus_gram_common.jl` | Pi-based initial periapsis (root of 0.5 rho_GRAM V^3 = Q_target at the two-body first-periapsis site), per-cell fits to GRAM mean density (piecewise primary, single exponential secondary), GRAM truth configuration, pass table (adds longitude and local scale height), cell task, process-route launcher |
| `venus_gram_generality.jl` | Driver, `--phase=guess\|convergence\|convergence_refined\|pilot\|truth\|predict`, `--fit=piecewise\|single` |
| `venus_gram_apoapsis_check.py`, `venus_gram_pi_cells.py` | D1V/D2V and the Pi table with the GRAM local scale height (companion: the start-fit H) |
| `venus_gram_diagnostics.py` | Criteria G1-G5, fit-induced density error per pass, Pi-rule site error, side-by-side with the exponential run |

`venus_analyze.py` is reused unchanged (the results subdirectory stays `venus`, so use a separate results root). Commands (repo root, one Julia job at a time, memory capped; `R` is the GRAM results root, `S` a second root for `--fit=single` whose `truth`, `cells.csv`, `fit_nodes.csv`, `truth_run_summary.csv`, `planet_constants.csv`, `body_distances.csv`, `estimator_numerics.toml` and `truth_mc.toml` are symlinks to `$R/venus`):

```bash
D=benchmarks/studies/periapsis_border_switchability; R=<results root>; run="systemd-run --user --scope -p MemoryMax=16G julia --project=."
for ph in guess convergence convergence_refined; do $run $D/venus_gram_generality.jl --phase=$ph --results-root=$R; done   # no truth
$run $D/venus_gram_generality.jl --phase=truth --workers=4 --results-root=$R     # 32 cells, process route
$run $D/venus_gram_generality.jl --phase=predict --results-root=$R               # frozen v3 and v5, piecewise fit
$run $D/venus_gram_generality.jl --phase=predict --fit=single --results-root=$S  # secondary single-exponential fit
python3 -I $D/venus_analyze.py --results-root $R; python3 -I $D/venus_gram_apoapsis_check.py --results-root $R; python3 -I $D/venus_gram_pi_cells.py --results-root $R
python3 -I $D/venus_gram_diagnostics.py --results-root $R --single-root $S --exp-root <exponential Venus root>
```

Headline (details and limits in the results note): the same 17 of 32 cells decay as with the exponential atmosphere, altitude predictions and pass-1 offsets are unchanged (V1, V2 pass), but the peak-Q error of the estimators rises from 1.4-24% to 5-540% over the bins (V3 fails from 5 d on) because an altitude-only density model cannot represent GRAM's longitude and local-time structure; the apoapsis-decay check keeps its median with the local scale height (D1V 1.017) but not its 5% tightness (D2V 0.80); the Pi-rule initial periapsis lands the first-pass Q in the band in 16 of 32 cells (the Sun's half-orbit tide shifts the periapsis by up to 3.6 km); the local scale height varies 2.1 to 2.7 times across one periapsis band and 3 to 44 km over a campaign.

## Density buffer (analysis)

Can a k-sigma buffer from previously measured densities replace a fixed periapsis-altitude margin (Q = 0.5 rho V^3, so a density error maps to an altitude buffer of about H_local k sigma_lnrho)? Python-only analysis of the Venus-GRAM mean-density truth passes (`density_buffer_analysis.py`; pre-registration and results notes `tasks/2026-10-15_periapsis-region-border-targeting/density_buffer_preregistration.md` and `density_buffer_results.md`). The quantity is y_n = ln(rho_n / rho_pw,n) against the start-fit model; predictors are the start model alone, persistence, EWMA (N = 3, 5, 10) and persistence plus a longitude-slope term; sigma is leave-one-cell-out.

```bash
python3 -I benchmarks/studies/periapsis_border_switchability/density_buffer_analysis.py --results-root <GRAM results root> --selftest
```

Headline: persistence gives sigma 0.054 ln (0.24 km), a 99% buffer of 0.65 km against the fixed 1000 m (decay cost 12% against 20% per pass); EWMA is worse at every N; coverage is over-conservative at k = 1 and near nominal at k = 1.64 and 2.33; the error is autocorrelated (lag-1 0.875). Limits: GRAM MEAN density has no random dispersion (sigma is a lower bound), one epoch, 15 completed cells flown uncontrolled, only 6 with long in-window runs.

## Stage F

Does a margin derived from previously measured densities replace the fixed 1000 m periapsis margin, with native GRAM density as truth, at Venus, Mars and Earth, with GRAM mean and GRAM perturbed density on the same cells? Plan (pre-registered, with addenda): `tasks/2026-10-15_periapsis-region-border-targeting/stagef_plan.md`; results: `.../stagef_results.md`. No exponential or US76 fallback is used. Needs native GRAM and SPICE as in 'Venus with GRAM' (Earth GRAM data is 11 GB).

| File | Role |
|---|---|
| `stagef_config.toml` | Every named value, labeled SOURCED / DERIVED / ASSUMED; per-planet tables merged with `[common.*]`; analysis thresholds (identical across planets) |
| `stagef_cells.py` | Seeded cell draw: 3 epochs x 4 inclinations x 4 cells per planet, 24 development and 24 held out |
| `stagef_common.jl` | Reuses `venus_gram_common.jl` and the frozen estimators; adds the per-planet setup, the perturbed-density arm (`SPACEAGORA_GRAM_DENSITY_PERTURBATION=pass`, reseed per pass), the in-run density measurement, the fixed / adaptive margin policy wrapper, and the self-contained case functions |
| `stagef_run.jl` | Driver, `--phase=feasibility\|cal\|f1\|f2\|selftest`: one `run_monte_carlo` campaign per phase over the planets and groups given, process route asserted |
| `stagef_analysis.py` | Tests R0 to R7, development-cell rule parameters (`f2params`), F2 comparison (`f2`, `f2track`), cross-host compare (`xhost`), `selftest` |

```bash
D=benchmarks/studies/periapsis_border_switchability; R=<results root>; run="systemd-run --user --scope -p MemoryMax=16G julia --project=."
python3 -I $D/stagef_cells.py --results-root $R                       # cell lists, before any run
$run $D/stagef_run.jl --phase=cal --results-root=$R --workers=4       # Pi-rule altitude, start fit, mean-density calibration (all feasible planets, one campaign)
$run $D/stagef_run.jl --phase=f1  --results-root=$R --workers=4       # open and fixed-margin keeper x mean and perturbed, per cell
python3 -I $D/stagef_analysis.py --phase=f1 --results-root $R         # R0-R7
python3 -I $D/stagef_analysis.py --phase=f2params --results-root $R   # predictor and error quantiles from development cells only
$run $D/stagef_run.jl --phase=f2  --results-root=$R --workers=4       # adaptive margin on held-out cells
python3 -I $D/stagef_analysis.py --phase=f2 --results-root $R
```

`--planets=` and `--groups=planet:epoch,...` split a phase across hosts (whole planet-epoch groups, so mean/perturbed pairs stay together); each worker is about 2.7 GB resident, so a 16 GB scope allows about 5 workers.

Headline (details and limits in the results note): the GRAM perturbation is white across passes at all three planets (lag-1 of ln r about 0) with std 0.125 (Venus), 0.32 (Mars), 0.07 (Earth), so EWMA beats persistence by 14 to 25% (the Stage E ranking inverts) and the buffer is planet-specific (95% quantile 0.24, 0.52, 0.13 scale heights, ratio 4.1). On held-out cells the adaptive corridor lowers soft-limit violations from 29.6% to 2.4% (Venus) and from 0.9% to 0.2% (Earth) at equal dV, and loses at Mars, where the random dispersion exceeds the corridor. At Venus most of the gain comes from tracking the systematic density drift (UNFROZEN diagnostic), at Earth from the quantile buffer.

## Stage F at 20x20

Stage F used truth gravity degree 3, order 0 at Earth and Mars (copied from earlier configs; Venus was already 20x20). This rerun uses degree 20, order 20 (Mars50c, GGM05C) for Earth and Mars, with the same cells, seeds, thresholds and rules, and the controller predictor v5 with tesseral terms over the same field as the truth. Plan addendum (pre-registered): `tasks/2026-10-15_periapsis-region-border-targeting/stagef20_plan.md`; results: `stagef20_results.md` next to it.

| File | Role |
|---|---|
| `stagef20_config.toml` | `stagef_config.toml` with only `results_subdir` and the Mars and Earth degree and order changed |
| `stagef20_cells.py` | Copy of `stagef_cells.py` reading the new config (same seed; the cell lists are byte-identical to Stage F) |
| `stagef20_run.jl` | Copy of `stagef_run.jl` whose config reads go to `stagef20_config.toml`; same phases and options |
| `stagef20_analysis.py` | Copy of `stagef_analysis.py` reading the new config; same tests, plus the same `f2params`, `f2`, `f2track`, `xhost` phases |

```bash
D=benchmarks/studies/periapsis_border_switchability; R=<results root>; run="systemd-run --user --scope -p MemoryMax=16G julia --project=."
python3 -I $D/stagef20_cells.py --results-root $R
$run $D/stagef20_run.jl --phase=cal --results-root=$R --planets=mars,earth --workers=4
$run $D/stagef20_run.jl --phase=f1  --results-root=$R --planets=mars,earth --workers=4
python3 -I $D/stagef20_analysis.py --phase=f1 --results-root $R
python3 -I $D/stagef20_analysis.py --phase=f2params --results-root $R
$run $D/stagef20_run.jl --phase=f2 --results-root=$R --planets=mars,earth --workers=4
python3 -I $D/stagef20_analysis.py --phase=f2 --results-root $R
```

Venus is not rerun: copy its Stage F cells, calibrations, runs, adaptive parameters and summary rows into `<results root>/stagef20/venus` so the three-planet tables and R7 are computed by the same script.

Headline (details and limits in the results note): no pre-registered verdict flips. The Mars adaptive-margin LOSS survives (4 hard-limit passes against 1 for the fixed margin, dV 6.4 times higher, corridor collapsed at every decision), the Venus and Earth WINS survive, and R7 planet generality still fails (q0.95 buffer 0.24, 0.51, 0.13 scale heights). The Mars horizon stays 30 days; a 60-day Mars rerun was not run.

## Stage G

Stage F's adaptive rule buffers both corridor edges by one-sided error quantiles; at Mars (20x20) the error spread exceeds the band and the corridor collapses at every decision. Stage G tests an asymmetric band: the thermal-protecting (high-heat) edge keeps its quantile buffer e_hi, the low-heat edge has no buffer (nominal predicted Q_low), and a band narrower than the 200 m minimum is anchored at the thermal edge instead of collapsing to its center. It runs on the held-out cells of the Stage F at 20x20 tree and is compared with the existing fixed 1000 m margin and Stage F adaptive runs (F1 is not rerun). Plan (pre-registered): `tasks/2026-10-15_periapsis-region-border-targeting/stageg_plan.md`.

| File | Role |
|---|---|
| `stageg_common.jl` | Includes `stagef_common.jl` unchanged; adds the asymmetric edge rule (`stageg_edges`), its policy wrapper and the run case function |
| `stageg_run.jl` | Driver, `--phase=f2g\|selftest`; reads `stagef20_config.toml` and the Stage F at 20x20 results tree, writes `runs/adaptive_asym_pert_<tag>/` |
| `stageg_analysis.py` | Pre-registered criteria G1 to G6 against fixed and Stage F adaptive (reuses `stagef20_analysis.py`) |

```bash
D=benchmarks/studies/periapsis_border_switchability; R=<results root>; run="systemd-run --user --scope -p MemoryMax=16G julia --project=."
$run $D/stageg_run.jl --phase=selftest --results-root=$R
$run $D/stageg_run.jl --phase=f2g --results-root=$R --workers=4
python3 -I $D/stageg_analysis.py --results-root $R
```

Headline (details in the results note): on the 24 held-out cells per planet the asymmetric band never collapses at Mars (0 of its decisions, against 1461 for Stage F adaptive) and removes the hard-limit passes (0, against 1 for the fixed margin and 4 for Stage F adaptive) and almost all soft violations (1.5%, against 9.8% and 6.9%). It fails the pre-registered Mars dV criterion: mean dV 1.42 m/s is 3.3 times the fixed margin's 0.43 m/s (limit 2 times), at half of Stage F adaptive's 2.76 m/s. The accepted cost is large: 31.6% of Mars passes below the low-heat limit (fixed 9.4%) and a 22% slower median apoapsis decay. Earth and Venus do not lose (soft violations 0.15% and 1.0%, dV within 1.3% of fixed), with a 6% (Earth) and 16% (Venus) slower apoapsis decay than fixed.

## Stage H

Stage G's Mars dV gate failed because the thermal-buffered band is only about 0.30 ln wide and the controller burns almost every pass. Stage H relaxes the low-heat edge by a fraction of the development-cell lower-side error quantile: hi = alt(2 Q_low / V^3 exp(-(yhat + g))) with g = e_lo_frac * |e_lo| >= 0, which moves the edge up and WIDENS the band (the opposite sign of Stage F's negative e_lo, which narrows it). The thermal edge, near_passive policy (which already skips burns inside the band), 200 m anchoring and all other settings are Stage G's; g = 0 is Stage G. e_lo_frac is chosen on the Mars development cells only (grid 0.25, 0.5, 0.75, 1.0, ascending, first value meeting D1 to D4 = G1, G2, G3, G5 against the fixed margin on those cells), then evaluated once on the held-out cells of all three planets with the Stage G criteria. Plan (pre-registered): `tasks/2026-10-15_periapsis-region-border-targeting/stageh_plan.md`.

| File | Role |
|---|---|
| `stageh_common.jl` | Includes `stagef_common.jl` only; adds the relaxed edge rule (`stageh_edges`), its policy wrapper and the run case function |
| `stageh_run.jl` | Driver, `--phase=dev\|held\|selftest --frac=F`; writes `runs/adaptive_h<FFF>_pert_<tag>/` |
| `stageh_analysis.py` | `--mode dev\|held --frac F`: selection criteria D1 to D4, held-out G1 to G6 against fixed, Stage F adaptive and Stage G |

```bash
D=benchmarks/studies/periapsis_border_switchability; R=<results root>; run="systemd-run --user --scope -p MemoryMax=16G julia --project=."
$run $D/stageh_run.jl --phase=selftest --results-root=$R
$run $D/stageh_run.jl --phase=dev --frac=0.25 --planets=mars --workers=4 --results-root=$R
python3 -I $D/stageh_analysis.py --results-root $R --mode dev --frac 0.25 --planets mars
$run $D/stageh_run.jl --phase=held --frac=0.25 --workers=4 --results-root=$R
python3 -I $D/stageh_analysis.py --results-root $R --mode held --frac 0.25
```

Headline (details in the results note): e_lo_frac 0.25 was selected on the Mars development cells (dV 1.84x fixed, limit 2x). On the held-out cells it does not carry over: Mars mean dV 0.981 m/s is 2.28x the fixed margin's 0.430 (limit 2x, so G3 fails; Stage G 1.417, Stage F adaptive 2.764), with 0 hard-limit passes (fixed 1), soft violations 0.9% (fixed 9.8%), 6.1 burns per run (fixed 2.6) and 41.0% of passes below the low-heat limit (fixed 9.4%, Stage G 31.6%), median apoapsis decay 26% slower than fixed. Earth (WIN) and Venus (SAFETY_TRADE) do not lose and stay within 1.10x of fixed dV, with 6.6% (Earth) and 19.4% (Venus) slower decay.

## Stage I

Stages F to H flew a corridor of the study's own band (Q 0.14 to 0.32 W/cm2, ln width 0.827), which the Mars density error (spread 1.085 ln) exceeds, so the two-sided adaptive rule collapsed there. Stage I flies new Mars flights against the corridor of Falcone and Putnam 2023 (IEEE TAES 59(3), PR-DRL): 0.05 to 0.25 W/cm2, target 0.15 (ln width 1.609, spread/band 0.67 at Mars), with their four violation tiers (low below 0.05, soft 0.25 to 0.30, medium 0.30 to 0.45, hard above 0.45). The band values are SOURCED (paper); that their heat rate (on the Odyssey solar panels) and this study's vehicle Q are comparable in absolute terms is ASSUMED, not established. Everything else is Stage F at 20x20: the same 48 cells per planet (24 development, 24 held out, byte-identical lists), GRAM mean plus perturbation truth, 20x20 gravity, process route. Redone for the new band: first-pass calibration to Q = 0.15, the fixed 1000 m comparator, and the predictor, e_hi and e_lo (refit from re-flown development keeper flights). Plan (pre-registered): `tasks/2026-10-15_periapsis-region-border-targeting/stagei_plan.md`.

| File | Role |
|---|---|
| `stagei_config.toml` | `stagef20_config.toml` with only `results_subdir`, the `[common.band]` values and the `[stagei]` table changed |
| `stagei_common.jl` | Includes the Stage F, G and H common files unchanged; adds the calibration case (target Q = 0.15) and a dispatcher over the fixed, Stage F adaptive, Stage G and Stage H arms |
| `stagei_run.jl` | Driver, `--phase=cal\|dev\|held\|selftest`, process route asserted |
| `stagei_analysis.py` | `--phase integrity\|params\|held\|selftest`: calibration and epoch checks, development refit, paper-tier counts and the pre-registered criteria I1 to I5 (Mars) and J1, J2 (Earth, Venus) |

```bash
D=benchmarks/studies/periapsis_border_switchability; R=<results root>; run="systemd-run --user --scope -p MemoryMax=16G julia --project=."
python3 -I $D/stagef20_cells.py --results-root $R   # then copy <R>/stagef20/<planet>/cells_base.csv to <R>/stagei/<planet>/ (the lists are the Stage F20 cells)
python3 -I $D/stagei_analysis.py --phase selftest; $run $D/stagei_run.jl --phase=selftest --results-root=$R
$run $D/stagei_run.jl --phase=cal  --results-root=$R --workers=4
$run $D/stagei_run.jl --phase=dev  --planets=mars --results-root=$R --workers=4
python3 -I $D/stagei_analysis.py --phase params --planets mars --results-root $R
$run $D/stagei_run.jl --phase=held --planets=mars --results-root=$R --workers=4
python3 -I $D/stagei_analysis.py --phase held --planets mars --results-root $R   # and the same for earth,venus
```

Use a distinct `--host=` stamp for every `dev` and `held` call: the summary and `mc_*.toml` file names carry only the phase and the host stamp, so two calls with the same stamp overwrite each other's rows (the pass and decision tables are not affected).

Headline (details and limits in the project log; 24 held-out cells per planet): at the paper width the two-sided adaptive rule does not collapse at Mars (0 of 1223 decisions; development spread/band 0.69 against 1.31 at the old band), cuts the over-0.25 fraction from 3.9% (fixed 1000 m) to 0.9% with no hard-tier pass, but its mean dV is 0.91 m/s, 2.37 times the fixed arm's 0.385 m/s, so it fails the pre-registered 2x dV limit. The asymmetric arms G and H pass every Mars criterion (dV 1.45 times fixed, over-0.25 fraction 0.5%, no hard-tier pass) at the cost of 10% and 15% of passes below 0.05 W/cm2 (fixed 3.2%) and a 35% and 41% slower median apoapsis decay. Pre-registered Mars verdict ASYMMETRIC_ONLY(G, H). Earth and Venus do not lose with any of the three arms (dV within 1.4% of fixed); at Venus the fixed 1000 m margin itself flies 8.5% of its passes above 0.45 W/cm2 at this band and all adaptive arms remove them. The Mars horizon stays 30 days, and the paper's heat-rate convention is not shown equal to this study's.

## Stage L

Stage I gave four isolated arms at the paper corridor. Stage L replaces them by a continuous family, so a sweep traces the (dV, burns, apoapsis decay, violation fraction) frontier per planet: thermal-edge buffer e_hi(p) = the p-quantile (p in 0.80, 0.90, 0.95, 0.99) of the Stage I development one-pass-ahead error, low-heat offset g = s |e_lo| with signed s (s = -1 is Stage F's low-edge buffer, s = 0 Stage G, s > 0 the Stage H direction; s in -1, 0, +0.5), plus a fixed-margin baseline curve (500, 1000, 2000, 4000 m): 16 arms. The 16 arms are flown on 12 development half-cells per planet, a pre-registered rule selects at most 4 members (non-dominated at over-0.25 fraction <= 1%), and those go to the 24 held-out cells against Stage I's fixed 1000 m arm. Calibrations, cells and the development refit are read from the Stage I tree.

| File | Role |
|---|---|
| `stagel_config.toml` | `stagei_config.toml` with only `results_subdir` and a `[stagel]` table (grids, ceiling, thresholds) changed |
| `stagel_common.jl` | Includes the Stage F and H common files unchanged; adds the arm code and one case function (fixed margin m, or a family member through the unchanged Stage H edge rule with a signed offset) |
| `stagel_run.jl` | Driver, `--phase=repro\|dev\|held\|heldall\|selftest`, `--dry-run`, process route asserted |
| `stagel_analysis.py` | `--phase params\|repro\|select\|held\|diag\|integrity\|selftest`: family quantiles, bit-identity check against Stage I, development selection, criteria L1 and L2, L3 trends, per-pass decay versus heat rate |

```bash
D=benchmarks/studies/periapsis_border_switchability; R=<results root containing stagei/>; run="systemd-run --user --scope -p MemoryMax=16G julia --project=."
python3 -I $D/stagel_analysis.py --phase selftest; $run $D/stagel_run.jl --phase=selftest --results-root=$R
python3 -I $D/stagel_analysis.py --phase params --results-root $R
$run $D/stagel_run.jl --phase=repro --results-root=$R --host=lrepro; python3 -I $D/stagel_analysis.py --phase repro --results-root $R
$run $D/stagel_run.jl --phase=dev  --planets=mars --results-root=$R --host=ldev_mars --workers=4
python3 -I $D/stagel_analysis.py --phase select --planets mars --results-root $R
$run $D/stagel_run.jl --phase=held --planets=mars --results-root=$R --host=lheld_mars --workers=4
python3 -I $D/stagel_analysis.py --phase held --results-root $R   # then diag and integrity; select with all planets at once rewrites analysis/l_dev.json for those planets only
```

Use a distinct `--host=` stamp for every call. Headline (held-out cells, details and limits in the project log): at the paper width no family member that the development selection picks keeps the fixed arm's campaign speed while meeting the violation ceiling at Mars or Venus (L1 and L2-Venus fail on the 90% decay condition: selected Mars member 72% of fixed decay, Venus 77%), whereas at Earth the selected member passes with dV 1.003x and decay 0.977x; every low-heat-tolerant member trades decay for fewer over-0.25 passes, and decay per pass rises with the pass heat rate.

## Stage N

Two parts. (1) Offline: the irreducible versus recoverable corridor buffer and a chance-constrained edge rule replayed on the flown Stage I and F20 fixed-arm pass tables (`stagen_analysis.py`, frozen; `stagen_unfrozen.py`, written after the frozen results). Finding: the white GRAM perturbation is 81 to 93% of the buffer at the three planets, no causal estimator beats the controller's predictor by more than a few percent, and an empirical-quantile edge at 1% is calibrated at Mars in the perfect-execution replay (13 of 1272 passes over 0.25 W/cm2). (2) Flight (Mars): the same rule flown closed loop on the 24 held-out cells of Stages I and L.

| File | Role |
|---|---|
| `stagen_analysis.py`, `stagen_unfrozen.py` | Offline decomposition, ceiling, Kalman filter and replay; the supplement (not frozen) |
| `stagen_flight_config.toml` | `stagel_config.toml` with only `results_subdir` and a `[stagen]` table (epsilon, cap grid, thresholds) changed |
| `stagen_flight_shim.jl` | Compatibility with current main: the old guidance tick and ephemerides hook names were removed, the frozen Stage F-L files still use them; forwarders, loaded first |
| `stagen_flight_common.jl` | Online Kalman filter (same recursion as `stagen_analysis.kalman`), the edge rule (thermal edge `lo = alt(2 Q_high / V^3 exp(-(m + t s)))`, low-heat edge of Stage G capped to `lo + w`), a Stage H hook copy, the case function |
| `stagen_flight_run.jl` | Driver, `--phase=smoke\|dev\|held\|repro\|selftest`, `--dry-run`, process route asserted, runs on disk skipped |
| `stagen_flight_xcheck.jl` | Re-flies one held-out cell of the Stage I fixed, G and F arms with the frozen Stage L case function (bitwise comparison with the stored tables) |
| `stagen_flight_analysis.py` | `--phase params\|select\|held\|repro\|integrity\|selftest`: filter constants (asserted equal to the offline record), cap selection on the dev half-cells, criteria C1 to C5, closed-loop innovations versus the replay |

```bash
D=benchmarks/studies/periapsis_border_switchability; R=<results root containing stagei/ and stagel/>; run="systemd-run --user --scope -p MemoryMax=16G julia --project=."
python3 -I $D/stagen_flight_analysis.py --phase selftest; $run $D/stagen_flight_run.jl --phase=selftest --results-root=$R
python3 -I $D/stagen_flight_analysis.py --phase params --results-root $R --n-results <stagen_data/n_results.json>
$run $D/stagen_flight_run.jl --phase=smoke --results-root=$R --workers=2 --host=smoke
$run $D/stagen_flight_run.jl --phase=dev   --results-root=$R --workers=4 --host=ndev;  python3 -I $D/stagen_flight_analysis.py --phase select --results-root $R
$run $D/stagen_flight_run.jl --phase=held  --results-root=$R --workers=4 --host=nheld; python3 -I $D/stagen_flight_analysis.py --phase held --results-root $R --n-results <stagen_data/n_results.json>
```

Use a distinct `--host=` stamp for every call. First run on a fresh checkout of current main: `julia examples/setup.jl` (the studies load through the examples environment). Headline (Mars, 24 held-out cells; details, sources and deviations in the project log): the selected arm (cap 4000 m) has 3 of 1184 passes over 0.25 W/cm2 (fixed 0.0393, conservative against the 1% target), no hard pass, dV 1.54x fixed, but median apoapsis decay 0.704x fixed (mean Q 0.73x): C1 fails on calibration (too conservative), C4 fails as in Stage L. The estimator itself stays calibrated in closed loop (innovation exceedance 13 of 1160 against the replay's 13 of 1272); the replay's cost-free mean Q came from aiming at the edge, which a banded `near_passive` controller does not do. The comparators are statistically, not bitwise, reproduced on the current code.

## Stage N2

Stage N-flight found the 1% edge calibrated in closed loop but the realized violation too low (0.25%) and decay 0.704x, because `near_passive` aims at the band center. Stage N2 changes only the aim point (target at Q_high exp(-(c + g_a)), g_a in ln Q below the edge, clamped to the band) and the edge multiplier t, and calibrates both so that the predicted overall P(Q > 0.25) is 1% (single tier). Calibration uses the Stage N-flight dev pass tables: the controller's predicted position ln(Q_pred/Q_edge) is a noise-driven random walk, and a Monte Carlo of it (validated at the Stage N operating point) shows that with t unchanged no aim reaches 1%, so t is relaxed per aim offset.

| File | Role |
|---|---|
| `stagen2_config.toml` | `stagen_flight_config.toml` with only `results_subdir` and a `[stagen2]` table (grid, thresholds, Monte Carlo settings) changed |
| `stagen2_common.jl` | Includes `stagen_flight_common.jl` (shim, Stage H files, online filter); adds the edge-side aim target, the policy copy and the case function |
| `stagen2_run.jl` | Driver, `--phase=smoke\|dev\|held\|repro\|selftest`, `--dry-run`, process route asserted |
| `stagen2_analysis.py` | `--phase calibrate\|select\|held\|repro\|integrity\|selftest`: Monte Carlo calibration, dev selection, criteria C1 to C5, frontier comparison |

```bash
D=benchmarks/studies/periapsis_border_switchability; R=<results root containing stagei/, stagel/, stagen_flight/>; run="systemd-run --user --scope -p MemoryMax=16G julia --project=."
python3 -I $D/stagen2_analysis.py --phase selftest; $run $D/stagen2_run.jl --phase=selftest --results-root=$R
python3 -I $D/stagen2_analysis.py --phase calibrate --results-root $R
$run $D/stagen2_run.jl --phase=smoke --results-root=$R --workers=2 --host=n2smoke
$run $D/stagen2_run.jl --phase=dev   --results-root=$R --workers=4 --host=n2dev;  python3 -I $D/stagen2_analysis.py --phase select --results-root $R
$run $D/stagen2_run.jl --phase=held  --results-root=$R --workers=4 --host=n2held; $run $D/stagen2_run.jl --phase=repro --results-root=$R --workers=2 --host=n2repro
python3 -I $D/stagen2_analysis.py --phase held --results-root $R; python3 -I $D/stagen2_analysis.py --phase repro --results-root $R
```

Use a distinct `--host=` stamp for every call. Headline (Mars, 24 held-out cells; details in the project log): selected arm g_a 0.10, t 1.52: over-0.25 8 of 1237 (0.65%, CP90 [0.32%, 1.16%] contains 1%; fixed 3.9%), no hard pass, dV 1.30x fixed, decay 0.884x (C4 needs 0.90), mean Q 0.913x, 4.4 burns per run. The (dV, decay) frontier moved: the arm dominates Stages F, G, H, the Stage L fm04000 and F-equivalent members and the Stage N arm at an over fraction below 1%; only `fam_p090_sp050` (lower dV, 0.719x decay) stays on the frontier beside it.
