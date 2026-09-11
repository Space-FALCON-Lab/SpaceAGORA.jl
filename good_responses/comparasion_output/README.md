# Comparison Output — 30 orbits, `gve_sma`, `gve_ecc`, and `gve_inc` schedules

Re-run of the [comparison_test_report.md](../comparison_test_report.md) scenario with
**30 target orbits** (instead of 10) and the **`gve_sma`**, **`gve_ecc`**, and
**`gve_inc`** helper-selection schedulers (instead of `naive_next_entering`) for both
codes. Same geometry otherwise: 10 helpers @ 1050 km, target @ 1000 km, 200 km laser
range, B=100, P=10 kW. Target inclination is 0° for `gve_sma`/`gve_ecc`, and 1° for
`gve_inc` (helper inclination stays 0° throughout — a nonzero relative inclination is
needed for `gve_inc` to have anything to actually optimize).

Outputs are organized by scenario key (schedule, plus `_it<N>deg` suffix when target
inclination != 0°): `gve_sma/`, `gve_ecc/`, `gve_inc_it1.0deg/`, each containing:

- `prototype/` — Kuang's prototype code (`1_Kuang's Prototype Code`), run via
  `scripts/run_prototype_analysis.jl <schedule> [target_inc_deg]`. ~1700-2000 adaptive
  Vern9 steps, runtime ~6 s.
- `spaceagora/` — SpaceAGORA integration (`2_SpaceAGORA.jl/ORACLE/run_case2_laser_links.jl
  --schedule <schedule> --target-inclination-deg <deg> --orbits 30 --timeseries-points 20000`),
  post-processed by `scripts/run_spaceagora_analysis.jl <schedule> [target_inc_deg]` from
  the feather output in
  `spaceagora_raw/`. 18923 saved steps, runtime ~26 s.

[scripts/encounter_analysis.jl](scripts/encounter_analysis.jl) holds the shared,
code-agnostic plotting/analysis functions used by both adapter scripts. Both scripts
take the schedule name (`gve_sma`, `gve_ecc`, or `gve_inc`) as their first command-line
argument (default `gve_sma`) and the target inclination in degrees as their second
(default `0.0`).

Both codes and all three schedules found the same 4 range-feasible encounters (helpers
1–4 only; helpers 5–10 never come within 200 km over 30 orbits — geometric feasibility
is governed by orbit geometry, not the scheduler). Laser-on duty cycles differ
meaningfully between schedules (e.g. helper02: ~47% under `gve_sma`, ~35% under
`gve_ecc`, ~64% under `gve_inc` at 1° inclination), confirming the scheduler choice is
actually being exercised — see `<scenario>/prototype/plot3_contact_windows/contact_windows.md`
vs `<scenario>/spaceagora/plot3_contact_windows/contact_windows.md`.

## Deliverables (per pair/encounter), in each of `prototype/` and `spaceagora/`

1. **`plot1_relative_speed/`** — one PNG per range-feasible encounter:
   |v_target − v_helper| (relative-speed magnitude) vs. time since encounter start, plus
   `combined_relspeed_vs_time_since_encounter.png` overlaying every encounter on one plot.
2. **`plot2_range_and_rangerate/`** — one PNG per pair (full 30-orbit timeline): range
   (km, left axis) with horizontal reference lines at 0 and +200 km (max range — the
   within-range band), overlaid with signed range-rate (m/s, right axis; negative =
   closing, positive = opening); plus `combined_range_rangerate_vs_time_since_encounter.png`
   (two stacked panels — range on top, range-rate on bottom — all encounters overlaid,
   x-axis reset to time-since-encounter-start).
3. **`plot3_contact_windows/`** — `contact_windows.csv` / `.md` (per-encounter geometric
   contact-window duration vs. real laser-on duration and duty cycle) plus
   `contact_windows_bar.png` (grouped bar chart, x-axis labels in short `h02_vs_t01` form).

All titles use a slightly reduced font size (`titlefontsize=10`), 3 lines: main title,
code identity ("ORACLE prototype code" / "ORACLE with SpaceAGORA engine"), and the
scenario property list ("Helpers @ 1050 km, target @ 1000 km, 0° inclination, 200 km
laser range, B=100, P=10 kW, gve_sma" or "..., gve_ecc"), all centered at the same font
size. Every legend is a single horizontal row placed below its plot (`legend=:outerbottom,
legend_column=-1`; the twin-axis and 2-panel figures use a phantom/dedicated legend
series so both axes' entries share one row).

## Key numbers

### `gve_sma`

| Pair | Prototype window (s) | Prototype laser-on (s) | SpaceAGORA window (s) | SpaceAGORA laser-on (s) |
|---|---|---|---|---|
| helper01_vs_target | 2591.2 | 2591.2 (100%) | 2590.0 | 2580.0 (99.6%) |
| helper02_vs_target | 5253.3 | 2493.1 (47.5%) | 5280.0 | 2480.0 (47.0%) |
| helper03_vs_target | 5319.6 | 2914.5 (54.8%) | 5360.0 | 2900.0 (54.1%) |
| helper04_vs_target | 4399.2 | 2497.5 (56.8%) | 4483.6 | 2493.6 (55.6%) |

### `gve_ecc`

| Pair | Prototype window (s) | Prototype laser-on (s) | SpaceAGORA window (s) | SpaceAGORA laser-on (s) |
|---|---|---|---|---|
| helper01_vs_target | 2591.2 | 2591.2 (100%) | 2590.0 | 2580.0 (99.6%) |
| helper02_vs_target | 5253.6 | 1883.7 (35.9%) | 5280.0 | 1870.0 (35.4%) |
| helper03_vs_target | 5313.6 | 1443.8 (27.2%) | 5370.0 | 1430.0 (26.6%) |
| helper04_vs_target | 4484.9 | 2166.0 (48.3%) | 4523.6 | 2153.6 (47.6%) |

### `gve_inc_it1.0deg` (1° target inclination)

| Pair | Prototype window (s) | Prototype laser-on (s) | SpaceAGORA window (s) | SpaceAGORA laser-on (s) |
|---|---|---|---|---|
| helper01_vs_target | 2054.9 | 1573.6 (76.6%) | 2130.0 | 1570.0 (73.7%) |
| helper02_vs_target | 4636.7 | 3012.4 (65.0%) | 4700.0 | 2990.0 (63.6%) |
| helper03_vs_target | 4773.3 | 2718.0 (56.9%) | 4780.0 | 2700.0 (56.5%) |
| helper04_vs_target | 4549.8 | 2251.2 (49.5%) | 4553.6 | 2243.6 (49.3%) |

Note: the geometric windows here are shorter than the 0°-inclination cases above — the
target now has an out-of-plane component of motion relative to the (0°-inclination)
helpers, so conjunctions are less head-on and the 200 km range window is crossed faster.

## Notes / fixes required to run this scenario

- Prototype: `test16_options.jl` itself was **not modified** — its two dead-wiring quirks
  (`T_seconds` hardcoded to 63071s regardless of `opts.orbits`; `gve_schedule` hardcoded
  to `:none` regardless of `opts.schedule`) were worked around directly in
  `run_prototype_analysis.jl` by computing `T_seconds` from `ORBITS=30.0` and passing
  `gve_schedule=<schedule>` explicitly.
- SpaceAGORA: fixed a real bug in
  [`ORACLE/functions/3_Dynamics.jl`](../../2_SpaceAGORA.jl/ORACLE/functions/3_Dynamics.jl) —
  `num_steps_to_save`/`data_rate` were hardcoded to 1000 points regardless of the
  documented `--timeseries-points` CLI option. Now wired to `opts.timeseries_points`
  (mirrors the pattern already used in `Paper_Numerical_Verification_Code/run_verification.jl`),
  which was needed here to resolve the ~2.5–5 minute laser-on windows at adequate
  resolution (used `--timeseries-points 20000` ≈ 9.5 s/sample vs. the previous fixed
  ≈189 s/sample for a 30-orbit run).
- "Encounter" = a maximal continuous time interval where range ≤ 200 km (laser range),
  independent of whether the laser was actually firing during that interval.

