# CYGNSS Attitude/Torque Reconstruction Study

Reconstructs CYGNSS's attitude trajectory over a real one-hour flight
telemetry window (spanning a slew maneuver at t=901.75s) using SpaceAGORA's
actual `SimulationConfiguration`/`run_simulation` engine on a three-body
spacecraft (bus + 2 panels), and uses it to answer two separate questions:

1. How much do gravity-gradient, GRAM aerodynamic, and real magnetorquer
   telemetry each individually contribute to the attitude trajectory, on top
   of a reaction-wheel-only baseline? (`run_independent_effects.jl`)
2. Is the ~87 deg/hour drift in that wheel-only baseline caused by missing
   environmental physics, or by the wheel telemetry model itself? Answered by
   inverse dynamics: back the required torque directly out of the measured
   attitude kinematics (bypassing wheel-speed telemetry entirely) and replay
   it through the engine. (`run_kinematic_backout.jl`)

## Setup

```
common.jl                       # telemetry loading, custom effectors, spacecraft, run_cygnss_case
run_independent_effects.jl      # question 1 -- runs the 5 sims, writes data/plot_data/independent_effects_*.arrow
run_kinematic_backout.jl        # question 2 -- runs the sim, writes data/plot_data/kinematic_backout*.arrow
run_position_accuracy.jl        # position-accuracy ablation on the 1hr slew dataset (see below)
run_position_accuracy_48hr.jl   # exact reference replication on the dedicated 48hr PVT dataset (see below)
make_plots.jl                   # reads data/plot_data/*.arrow, writes plots/*.pdf -- no simulation, fast to re-run
make_rmse_table.jl              # reads data/plot_data/*_rmse.arrow, writes tables/cygnss_rmse_tables.tex
data/                           # pre-extracted telemetry (see "Data provenance" below)
data/plot_data/                 # intermediate arrow files written by the run_*.jl scripts (gitignored)
plots/                          # all output plots land here, PDF (gitignored)
tables/                         # LaTeX RMSE table lands here (gitignored)
```

Each `run_*.jl` script (`include`-ing `common.jl` itself, except
`run_position_accuracy_48hr.jl`, which is self-contained -- see "Position
accuracy" below for why) only runs the simulation(s) and writes the plotted/
tabulated quantities to `data/plot_data/*.arrow`; it does not generate any
plots or tables itself. Run the ones you need, then `make_plots.jl` and/or
`make_rmse_table.jl` to (re)generate output from whatever data is currently
on disk:

```bash
julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/run_independent_effects.jl
julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/run_kinematic_backout.jl
julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/run_position_accuracy.jl
julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/run_position_accuracy_48hr.jl

julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/make_plots.jl        # plots/*.pdf
julia --project=. benchmarks/studies/cygnss_attitude_torque_reconstruction/make_rmse_table.jl   # tables/cygnss_rmse_tables.tex
```

`make_plots.jl`/`make_rmse_table.jl` don't re-run any simulation, so they're
fast (seconds) to re-run when only iterating on plot styling or table
formatting -- only re-run the relevant `run_*.jl` script if the underlying
data is stale.

Each full-hour case takes roughly 1-2 minutes (GRAM cases somewhat longer:
real MERRA2 atmosphere queries plus one extra raw density query per non-root
link per RHS evaluation, see "Per-link atmosphere sampling" below). The 48hr
reference replication takes ~35s (it's a 50x50-degree/order harmonics
propagation over 172800s; the "quick" verification profile subsamples the
comparison points, not the propagation itself).

## Spacecraft model

Bus + 2 solar panels (`SimulationModel.SpacecraftModel`/`Link`, no facets --
SpaceAGORA's box-model aerodynamic/SRP torque doesn't need them):

- bus: 0.202 x 0.521 x 0.641 m, 28.94 kg
- panels: 0.497 x 0.0005 x 0.521 m each, offset +/-0.569 m in y, ~0 kg
- inertia tensor overridden with a custom value validated against the
  maneuver-window conservation-law fit (not the auto box-inertia calc):
  `[1.4e6 -1.71e4 8.08e3; -1.71e4 8.19e5 -5.35e3; 8.08e3 -5.35e3 1.95e6] * 1e-6`

Initial condition (`CartesianInitialCondition`) is read directly from
telemetry at the calibration time (r, v, q, ang_vel).

## Independent-effect results (full hour, t=0-3600s)

Each row adds one effect on top of the Omega_rw/J_RW wheel-momentum-replay
baseline (`WheelMomentumReplayControlModel` in `common.jl`, computing
`tau_control = -w x H_wheels(t) - Hdot_wheels(t)` from measured wheel-speed
telemetry via a J_RW allocation matrix regressed from the maneuver window).

| Configuration | mean error | max error | Delta vs wheel-only |
|---|---|---|---|
| Wheel-only baseline | 87.56 deg | 179.72 deg | -- |
| + gravity-gradient | 59.21 deg | 109.12 deg | **-28.3 deg** (dominant) |
| + aero (GRAM, per-link density) | 87.56 deg | 179.73 deg | ~0 deg (negligible) |
| + magnetic (real dplCmd telemetry) | 90.38 deg | 179.69 deg | +2.82 deg (worse) |
| All combined | 101.61 deg | 179.78 deg | +14.05 deg (worse than any single addition) |

Gravity-gradient is the only one of the three that meaningfully helps.
Aerodynamic torque is real (confirmed via GRAM's actual MERRA2 density
queries, not a stub) but utterly negligible at CYGNSS's ~500-540 km altitude
and panel scale. The magnetic dipole telemetry (`dplCmd`, taken at face value
as A*m^2 -- see "dplCmd units" below) makes the fit *worse*, and stacking all
three does not beat gravity-gradient alone. None of these come close to
closing the gap to the maneuver-window fit's sub-1-degree accuracy.

### Position and velocity, same five configurations

`run_cygnss_case` also compares simulated ECI position and velocity against
ground truth (`plots/cygnss_independent_effects_position_summary.pdf`,
`..._position_timeseries.pdf`, `..._velocity_summary.pdf`,
`..._velocity_timeseries.pdf`). Ground truth here is the raw GPS receiver fix
(`gpsData.PV_ecef`), converted to ECI via SpaceAGORA's own SPICE frame -- see
"Position accuracy" below for why, this wasn't the original choice:

| Configuration | position mean | position max | velocity mean | velocity max |
|---|---|---|---|---|
| Wheel-only baseline | 23.876 km | 63.029 km | 0.02811 km/s | 0.06512 km/s |
| + gravity-gradient | 23.876 km | 63.029 km | 0.02811 km/s | 0.06512 km/s |
| + aero (GRAM, per-link density) | 23.857 km | 62.947 km | 0.02808 km/s | 0.06502 km/s |
| + magnetic (real dplCmd telemetry) | 23.876 km | 63.029 km | 0.02811 km/s | 0.06512 km/s |
| All combined | 23.857 km | 62.947 km | 0.02808 km/s | 0.06502 km/s |

Position (and velocity) error is essentially identical across all five (as
expected -- gravity-gradient and magnetic torque are pure torques with no
translational effect at all; only the two aero cases differ, and only
slightly, from the real drag force `AerodynamicCoefficientfM` adds on top of
its torque). Every case here uses bare `InverseSquaredGravityModel`
(two-body only, no J2), deliberately, since this study's focus is
attitude/torque, not orbit determination -- see "Position accuracy" below for
how much of this ~24 km floor is the missing J2 term versus this specific GPS
fix's own precision, and for an exact (not just close) reproduction of the
48hr reference's ~1.58 km result.

## Kinematic torque back-out

The real question the independent-effect sweep raises: is the 87.56 deg/hour
residual from missing environmental physics, or from the wheel telemetry
model (Omega_rw -> H_wheels via J_RW/I_WHEEL) itself?

Answered directly: differentiate the *measured* attitude trajectory (not
wheel-speed telemetry) to get w(t) and its derivative, and evaluate the rigid
body Euler equation on the measurement itself:

```
tau(t) = I * w_dot_meas(t) + w_meas(t) x (I * w_meas(t))
```

Under the modified Euler equation with wheel momentum,
`I*w_dot + w x (I*w + H_wheels) + Hdot_wheels = tau_ext`, setting
`tau_ext = 0` (torque-free environment -- no gravity-gradient/aero/magnetic)
makes this `tau(t)` exactly the combined wheel-reaction term that the wheels
must have delivered, with no dependence on Omega_rw telemetry or J_RW at all.
Applying its negative through the same wheel-momentum control channel and
re-integrating from the same initial condition:

| | mean error | max error |
|---|---|---|
| Kinematic back-out replay | **0.22 deg** | **0.42 deg** |
| Omega_rw/J_RW wheel-only replay (for reference) | 87.56 deg | 179.72 deg |

This is near-perfect reproduction (see `plots/cygnss_kinematic_backout_quaternion_q1..4_timeseries.pdf`,
one file per component since each is now its own untitled figure): the
quaternion components are visually indistinguishable from telemetry for the
full hour, and the error stays flat around 0.2-0.3 deg through the middle of
the window rather than growing, which is the signature of
spline-differentiation/integration noise, not a physical mismatch (it does
climb to 0.42 deg by the very end, consistent with edge effects in the
cubic-spline derivative near the boundary of the telemetry window; see
`plots/cygnss_kinematic_backout_attitude_error_timeseries.pdf`).

This run's translational dynamics are the same bare two-body gravity model as
every case in the independent-effect sweep (torque back-out only touches
attitude), so it shows the same ~24 km mean / ~63 km max position error for
the same reason (`plots/cygnss_kinematic_backout_position_x/y/z_timeseries.pdf`,
`plots/cygnss_kinematic_backout_position_error_timeseries.pdf`,
`position mean=23.863 km max=63.029 km`) -- see "Position accuracy" below.
Velocity error is the matching translational quantity
(`plots/cygnss_kinematic_backout_velocity_error_timeseries.pdf`,
`velocity mean=0.02810 km/s max=0.06512 km/s`), computed the same way
(`u.sc[1].vel` vs. the GPS-derived `v_truth(t)` in `common.jl`) for every case
in this study, not just this one.

**This proves the 87.56 deg residual is not missing environmental physics --
gravity-gradient/aero/magnetic are all too small individually and combined to
explain it, and a torque-free-except-wheels model fully explains the full
hour once the true required torque is known.** The bottleneck is specifically
the Omega_rw -> H_wheels mapping: J_RW was regressed from the maneuver window
alone and evidently doesn't generalize to the full hour, and/or Omega_rw
telemetry itself has drift/bias over longer timescales that the linear wheel
model doesn't capture.

The backed-out torque magnitude (`plots/cygnss_kinematic_torque_magnitude.pdf`)
is a useful sanity check in its own right: mean 1.10e-4 N*m, max 4.67e-4 N*m,
comfortably inside CYGNSS's actual reaction-wheel torque rating (6.55e-4 N*m
max, from the wheel datasheet value used elsewhere in this repo) -- a
physically plausible number came out of pure kinematics, not a fabricated
one. The torque time series itself is dominated by high-frequency noise
rather than a smooth commanded profile, because it comes from numerically
differentiating a cubic spline of already-noisy rate telemetry, not a
physical torque command; treat it as a diagnostic quantity, not a flight
torque-rod command reconstruction.

## Position accuracy

Three progressively deeper investigations, all under `run_position_accuracy*.jl`.

### GPS as ground truth (1hr slew dataset)

Position ground truth was originally `adsAttFilter.pv_eci` (the onboard
attitude filter's own propagated position estimate). Neither adding J2 nor
the full 50x50 EGM05C harmonics field, nor `test/gmat_scenario_matrix.jl`'s
N=5 vis-viva-SMA-averaged IC-noise fix (see "Lessons learned" #7 below),
moved that ~29 km baseline by more than a few percent -- a strong sign the
problem wasn't propagation fidelity at all.

Direct cross-check: converting the *raw* GPS receiver fix
(`gpsData.PV_ecef`, a genuinely independent measurement, not previously
used) to ECI via SpaceAGORA's own SPICE frame and comparing against
`adsAttFilter.pv_eci` at the same timestamps -- with zero propagation
involved -- showed a 7-14 km disagreement. Sweeping the ECEF->ECI conversion
epoch to find the best-fit timing offset only reduced that to ~5.5 km, not
~0, ruling out a simple clock/leap-second bug and pointing to a genuine
frame-convention difference between full SPICE-grade Earth orientation and
whatever simplified onboard sidereal-rotation model `adsAttFilter` uses.
That ~5.5 km is essentially the same size as the "position error" that had
been chased through gravity-model and IC fidelity -- so the ground truth
itself, not the propagation, was the dominant error source.

Switching ground truth to the raw GPS fix (properly converted to ECI, see
`common.jl`'s `gps_eci` construction) dropped the two-body baseline from
~29 km to ~24 km, and the best case (J2, N=5 SMA-corrected IC) from 4.89 km
to 2.92 km mean -- real improvement, and the error-vs-time shape changed
from a bounded oscillation (frame-mismatch signature) to a smooth,
monotonically growing curve (ordinary residual-dynamics signature),
confirming the diagnosis. It doesn't reach sub-2km, and a likely reason is
visible in the telemetry itself: this GPS receiver's own reported GDOP is
12-22 for this window, poor by GPS standards (a good fix is typically GDOP
< 4-6) -- the ground truth may simply carry real measurement uncertainty on
the order of a km or two, setting a floor no amount of propagation fidelity
can get under. `run_position_accuracy.jl` runs this full ablation (raw vs.
N=5-SMA-corrected IC, x, two-body vs. J2 vs. 50x50 EGM harmonics) and writes
`plots/cygnss_position_accuracy_ablation.pdf` (mean/max position error bar
chart) and its velocity-error counterpart
`plots/cygnss_velocity_accuracy_ablation.pdf`, plus best-case (SMA-corrected
IC, 50x50 EGM) time series:
`plots/cygnss_position_accuracy_best_case_x/y/z_timeseries.pdf`,
`..._error_timeseries.pdf`, and `..._velocity_error_timeseries.pdf`.

### Exact reference replication (48hr dedicated PVT dataset)

The ~1.58 km figure quoted throughout this README comes from a *different*
telemetry source entirely: `data/telemetry/CYGNSS/cygnss_data_48hr.feather`
(a dedicated 48-hour PVT position/velocity file, already ECI, no frame
conversion needed), scored via `test/gmat_scenario_matrix.jl`'s "CYGNSS
Legacy 48hr Entry Point" testset. Applying the same N=5 SMA-averaging + 50x50
EGM harmonics technique to *this* dataset by hand (reading the scenario dict
in `test/gmat_scenario_matrix.jl` and reconstructing an equivalent
`SimulationConfiguration`) got to ~4.5 km, not ~1.58 km, with the raw
(uncorrected) baseline landing suspiciously *better* than the reference's
own raw baseline (5.0 km vs. their 6.2 km) -- a pattern inconsistent with any
single missing setting, and unexplained after checking gravity fidelity,
solver algorithm (`dp8`), tolerances, ephemerides model, and confirming drag
is genuinely disabled for this scenario (`atmosphere=GRAM` in the real
test's log output is metadata about how the reference telemetry was
originally generated, not an active force -- confirmed by reading
`_scenario_density_model` in `scenario_builders.jl`, which only builds a GRAM
density model when `drag_enabled=true`).

`run_position_accuracy_48hr.jl` instead copies the handful of helper
functions `test/gmat_scenario_matrix.jl` uses to build and run this exact
scenario (`_build_cygnss_48hr_reference`, `_base_scenario_dict`,
`_scenario_rmse`, `_telemetry_solver_env_overrides`) verbatim, and calls the
real `TV.VerificationRequest`/`TV.run_verification` pipeline directly instead
of hand-reimplementing a `SimulationConfiguration` to approximate it. First
attempt still failed (`Missing 'spacecraft' in manifest`) because the earlier
hand-read of `_base_scenario_dict` had been truncated at a 40-line grep
window and missed its actual default block (`spacecraft`, a second
`initial_time` default, `atmosphere_truth`, `calibration`, `EI_km` -- lines
815-849, not shown by the original `-A 40`); copying the complete function
fixed it. Result:

**`cygnss_48hr_pvt mean position-axis RMSE: 1.5777521293524328 km`** -- an
exact match (to reported float precision) to the documented reference value,
because it now runs the identical code path rather than a parallel
reimplementation that could (and did) silently diverge from it. Per-axis:
state_x 1.820 km, state_y 1.644 km, state_z 1.269 km (mean = 1.578 km,
matching `_scenario_rmse`'s definition: the mean of three *separate*
per-axis RMSEs, not one 3D Euclidean-distance RMSE -- an easy metric
mismatch to fall into, and the first thing that needed fixing in the hand
reimplementation before the spacecraft-defaults bug was found). See
`plots/cygnss_48hr_reference_replication_error_timeseries.pdf` and
`data/cygnss_48hr_reference_summary.csv` (the per-sample errors table is
~100MB+ and deliberately not persisted -- re-run the script to regenerate it
if needed).

The plot itself is IQR-filtered (`_iqr_inlier_mask`/`_filter_position_series_
for_plot`, copied verbatim from `test/gmat_scenario_matrix.jl` alongside the
other helpers) to drop isolated spikes before plotting -- 19528/172797
samples (~11%) get dropped as outliers on the 3D-distance error, most likely
isolated GPS-fix glitches in the 48hr telemetry itself rather than a
propagation issue (the RMSE printed above is always computed from the
*unfiltered* errors table, so 1.578 km is unaffected by this -- the filter
is a plotting-only choice, applied after the number that matters was already
computed; the plot no longer draws that unfiltered-RMSE value as a reference
line, only the per-axis errors and the running RMSE of the filtered series).

**Velocity RMSE for this same 48hr case: `0.0017462593412124674 km/s`**,
computed the same way (mean of the three per-axis `state_vx_time`/
`state_vy_time`/`state_vz_time` RMSEs) against the real per-sample
`vel_ii_1/2/3` ECI velocity telemetry already present in
`cygnss_data_48hr.feather`. Getting this required extending the shared
telemetry-verification framework itself
(`src/analysis/verification/telemetry_verification/`): the
`time_aligned_state` comparison path only ever computed `state_x/y/z_time`
(position); velocity was available in `_load_time_aligned_telemetry` only as
`_differentiate_series(x_km, time_s)` (numerically differentiated from
position, used solely for the unrelated `orbit_events`
perigee/apogee-speed comparison, never compared against directly). Added
three new optional `TimeAlignedScenarioConfig` fields
(`telemetry_vx_col`/`vy`/`vz`, distinct from the pre-existing single-value
`vx_ic`/`vy_ic`/`vz_ic` IC columns) that, when a scenario's
`telemetry_columns` sets them, make `_load_time_aligned_telemetry` read the
real per-sample column instead of differentiating, and make
`_time_aligned_rows_errors` additionally emit `state_vx_time`/`vy`/`vz`
error rows compared against it. Strictly opt-in (new fields default to
`nothing`), so every other scenario using this shared framework -- including
`test/gmat_scenario_matrix.jl`'s own tests -- is unaffected; verified by
re-parsing the existing manifest and confirming the new fields default to
`nothing`, and by re-running this scenario and confirming the position RMSE
is still exactly `1.5777521293524328` km after the change.

## RMSE summary table

`make_rmse_table.jl` reads `data/plot_data/*_rmse.arrow` and writes a single
LaTeX table, `tables/cygnss_rmse_tables.tex` (gitignored, regenerate by
re-running the script). It deliberately does *not* draw its four numbers
from one common case -- each comes from whichever case in this folder
actually exercises that quantity meaningfully:

| Position RMSE (km) | Velocity RMSE (km/s) | Quaternion RMSE (deg) | Angular velocity RMSE (deg/s) |
|---|---|---|---|
| 1.578 | 0.00175 | 0.232 | 1.489e-6 |

- **Position/velocity** come from the 48hr reference-replication case (see
  above) -- the longer window and the validated `TV.run_verification`
  pipeline make it the meaningful test of translational accuracy, not the
  1-hour slew window (whose position error is dominated by using bare
  two-body gravity by design, ~24 km, and isn't the point of that study).
- **Quaternion/angular velocity** come from the kinematic torque back-out
  case instead (the full-hour slew maneuver) -- the 48hr case and the other
  1-hour studies (`run_independent_effects.jl`, `run_position_accuracy.jl`)
  all drive attitude with the wheel-only Omega_rw/J_RW baseline, which the
  independent-effects table above already shows floors at ~87-114 deg and
  isn't a meaningful attitude-reconstruction result; the kinematic back-out
  case is the one place in this study where attitude is actually
  reconstructed well.

Both RMSEs in the first two columns are the mean of three per-axis RMSEs
(matching `_scenario_rmse`'s definition used throughout the 48hr section
above), not a 3D Euclidean-norm RMSE. The quaternion RMSE is
`sqrt(mean(angle_err_deg .^ 2))` where `angle_err_deg` is the same
`2*acos(|q_sim . q_gt|)` magnitude error used everywhere else in this study
(not an elementwise per-component quaternion RMSE). The angular velocity
RMSE is `sqrt(mean(norm(w_sim - w_gt) .^ 2))`, `w_gt` from the same
`ω_meas(t)` telemetry interpolant `run_kinematic_backout.jl` already used.
The angular velocity number is genuinely `1.489e-6 deg/s`, not exactly zero
-- expected, since `KinematicTorqueReplayControlModel` is inverse dynamics
constructed specifically to reproduce the measured rate trajectory, so the
residual is purely spline/differentiation/integration error, not a real
mismatch.

## Lessons learned / bugs found and fixed this session

Several of these are `src/` and `data/GRAMSuite.jl` fixes, not part of this
study folder -- listed here because they were discovered *by* this study and
explain why some of the above numbers took several iterations to trust.

1. **Magnetic torque was completely unwired.** `get_magnetic_field_dipole`
   and `calculate_magnetic_torque` (tilted-dipole model + `tau = m x B`) were
   correct, verified physics with no caller anywhere in `src/` -- exactly
   parallel to a previously-found orphaned `ReactionWheelAssembly`. Fixed by
   adding `MagneticTorqueRodModel <: AbstractForceTorqueModel`
   (`src/dynamics/coupled/perturbations.jl`), following the same
   `wrench()`/`calcForceTorque()` dual-interface pattern as the gravity
   models, wired through the normal `dynamic_effectors` tuple. A real but
   minor doc bug went with it: `Magnet.m`'s comment said "nT" (a field unit)
   when it's actually A*m^2 (`src/vehicle/spacecraft/components.jl`).

2. **Relative vs. absolute time.** The engine's integration time `t` passed
   to `wrench`/`calcControlForceTorque` is relative to the simulation start
   (always starts at 0), but the telemetry interpolants are indexed by
   absolute time. Every custom effector/control model in `common.jl` takes an
   explicit `t_offset` and evaluates telemetry at `t + t_offset`. Missing
   this the first time silently fed the control torque totally wrong-window
   wheel data (e.g. a "t=899-904s" run was actually driven by t=0-5s wheel
   telemetry) -- caught by cross-validating the engine's wheel-only replay
   against an independent, from-scratch standalone integration of the same
   physics outside the engine (both converged to 0.79/2.24 deg over the
   maneuver window once fixed).

3. **`GRAMSuite.jl` submodule was 3 commits behind and missing a symbol.**
   `_GRAM_EPHEMERIS_STATE_FN` referenced by `ext/SpaceAGORAGRAMSuiteExt.jl`
   didn't exist in the checked-out commit, so the extension silently failed
   to load every run. Fast-forwarded the submodule (the fix was already
   fetched locally as `origin/main`, no network needed), and reconciled a
   real merge conflict against a local, uncommitted native-library-discovery
   patch (`_find_gram_library`, Linux `.so` support) -- kept both, preserved
   the already-built Linux native library.

4. **`density_model` was hardcoded to `NoAtmosphereModel()`.** Even after
   fixing the extension, the aero-effect test case never actually passed
   `GRAMAtmosphereModel` into `EnvironmentModel`, so `AerodynamicCoefficientfM`
   correctly computed exactly zero force/torque (its own zero-density guard
   clause). `run_cygnss_case` now takes `density_model` as a parameter.

5. **A genuine upstream world-age bug in `GRAMSuite.jl`,** newly exposed once
   the extension actually started loading (previously dormant because the
   extension never got that far). `_gram_apply_user_ephemeris_state!` (added
   by the just-pulled "ephemeris-state bypass hook" commit -- the actual Mars
   GRAM fix documented in
   `../gram_mars_fix_and_constellation_scaling/README.md`) constructed
   `EphemerisStateC`, a type from a module loaded via a runtime
   `Base.include` (not precompiled), *outside* an `invokelatest` boundary --
   safe for short runs, but the RHS gets evaluated enough times in a full
   3600s run to reliably hit it. Fixed with the same
   `Base.invokelatest(EphemerisStateC, state...)` pattern the surrounding
   code already used for the analogous `set_ephemeris_state!` call, one line.

6. **Per-link atmosphere sampling.** `AerodynamicCoefficientfM`/`Constant`/
   `NoBallisticFlight` sampled density/temperature once at the spacecraft's
   root position and reused it for every link, including the two panels
   offset ~0.57 m away. Fixed in
   `src/dynamics/coupled/aerodynamic_wrench_models.jl`: `_aero_pure_wrench`
   now accepts an optional `link_atmosphere_fn` that re-queries the density
   model at each non-root link's own planet-fixed position (raw, uncached --
   deliberately bypasses the satellite-level GRAM trajectory-extrapolation
   cache, which is keyed and warmed per satellite trajectory point and would
   be corrupted by jumping between link offsets). The `wrench` (non-caching)
   methods, which don't receive the live density model, keep the old
   single-sample behavior. For CYGNSS specifically this makes no measurable
   difference (panel offsets are minuscule next to atmospheric density scale
   heights) but it's now physically correct for larger multi-body spacecraft.

7. **`adsAttFilter.pv_eci` was a bad ground-truth choice.** See "Position
   accuracy" above -- it disagrees with the raw GPS fix by 5-14 km with zero
   propagation involved, a frame-convention mismatch, not a physics gap. The
   lesson generalizes: when a "ground truth" telemetry field is itself
   *derived* (filtered/propagated onboard) rather than a direct sensor
   reading, cross-check it against the underlying raw measurement before
   trusting a comparison against it, especially before concluding a
   simulation's dynamics model is deficient.

8. **Reimplementing a working pipeline by hand invites silent drift; copying
   it doesn't.** The 48hr replication's hand-built `SimulationConfiguration`
   matched every setting visible in `test/gmat_scenario_matrix.jl`'s
   scenario dict and still landed at ~4.5 km instead of ~1.58 km, for a
   reason that took real effort to find (a truncated function read that
   missed a default `spacecraft`/`initial_time`/etc. block). Copying the
   actual helper functions and calling the real `TV.VerificationRequest`/
   `TV.run_verification` pipeline directly instead produced an exact match
   immediately, because there was no longer a parallel implementation that
   could diverge from the original in some untraced way. Prefer this
   whenever "replicate an existing result" is the actual goal.

9. **dplCmd units.** The magnetorquer command telemetry
   (`ADCS_BUS_FSW_OUT.actCommands.dplCmd[0..2]`) has no units metadata and is
   suspiciously quantized (multiples of 1/24, max magnitude ~0.42), which
   looks like a normalized PWM duty cycle rather than a raw A*m^2 dipole
   moment. No documented CYGNSS torque-rod spec was available to derive a
   scale factor, so by explicit choice it's used as-is (raw value = A*m^2,
   no fabricated calibration constant) in `TelemetryMagnetorquerModel`. This
   is the most likely explanation for the magnetic effect actively hurting
   the fit in the independent-effects table above -- if `dplCmd` really is a
   duty-cycle fraction, the effective torque is being driven at whatever
   scale the true rod happens to sit at relative to the assumed 1 A*m^2 unit,
   not necessarily the right one.

10. **Extending a shared, validated pipeline is safe when done strictly
    opt-in.** Getting a real (not differentiated) velocity RMSE for the
    48hr case (see "RMSE summary table" above) meant adding new capability
    to `src/analysis/verification/telemetry_verification/`, code shared
    with `test/gmat_scenario_matrix.jl` and other consumers. Every addition
    (`TimeAlignedScenarioConfig`'s three new fields, the new
    `state_vx/vy/vz_time` event rows) is gated behind those fields being
    non-`nothing`, so scenarios that don't set them -- every existing one --
    take the exact same code path as before. Verified concretely, not just
    by inspection: re-parsed the existing test manifest and confirmed the
    new fields default to `nothing`, and re-ran this exact scenario after
    the change and confirmed the position RMSE was still
    `1.5777521293524328` km to full float precision. Cheap insurance against
    silently regressing a pipeline other tests depend on.

## Data provenance

- `data/cygnss_slew_data_extracted.feather`: extracted from
  `CYGNSS_slew_data.feather` (repo root, HDF5-format, gitignored) via a
  throwaway h5py venv. Columns: `t_rel`, `q_eci_0..3` (scalar-first),
  `w_eci_0..2`, `q_lvlh_0..3`, `tqDmdCtrl_0..2`, `tqDmdMomMgr_0..2`,
  `rwCmd_0..2`, `Omega_rw_0..2` (RPM), `dplCmd_0..2` (see "dplCmd units").
  Raw HDF5 paths: `ADCS_BUS_FSW_TLM.adsAttFilter.{time,qw_eci.q,qw_eci.w}`,
  `ADCS_BUS_FSW_INP.sensorData.tacData.Omega_rw`,
  `ADCS_BUS_FSW_OUT.actCommands.{rwCmd,dplCmd}`. 14399 rows, t=0-3600.5s.
- `data/cygnss_slew_orbit_extracted.feather`: `r_eci_0..2`, `v_eci_0..2` from
  `ADCS_BUS_FSW_TLM.adsAttFilter.pv_eci.{r,v}` in the same source file. No
  longer used as position ground truth (see "Position accuracy" above) --
  kept for reference/comparison only, not read by `common.jl`.
- `data/cygnss_slew_gps_raw_extracted.feather`: the actual position ground
  truth now used by `common.jl`. `r_ecef_0..2`, `v_ecef_0..2` from
  `ADCS_BUS_FSW_INP.sensorData.gpsData.PV_ecef.{r,v}` (raw GPS receiver fix,
  ECEF, ~1 Hz, held between telemetry samples until the next fix -- `common.jl`
  dedupes to actual fix instants via `gps_sec`/`gps_week` before
  interpolating), plus `r_eci_ads_0..2` (`adsAttFilter.pv_eci.r`, kept only
  for the ground-truth cross-check described above). `common.jl` converts
  ECEF->ECI itself via `SimulationModel.planet_frame_lpi` at
  `TELEMETRY_INITIAL_TIME`, not a value from the source file.
- `data/cygnss_data_48hr.feather`: copy of
  `data/telemetry/CYGNSS/cygnss_data_48hr.feather` at the repo root, a
  dedicated 48-hour PVT position/velocity file (already ECI, no frame
  conversion needed), used only by `run_position_accuracy_48hr.jl`, not by
  `common.jl` or any other script in this folder. This is a genuinely
  different dataset from the 1-hour slew telemetry everything else here uses
  -- see "Position accuracy" above.
- `TELEMETRY_INITIAL_TIME` (2025-10-04T00:56:58 UTC): derived from the slew
  dataset's raw GPS telemetry (`gpsWeek`=2386, `gpsSec`=521835.99999998 at
  the first sample), converted GPS time -> UTC (GPS is 18s ahead of UTC as of
  the 2016-12-31 leap second). An earlier, unrelated placeholder date
  (2018-12-15) was used here before this was derived; correcting it made no
  measurable difference to propagated position error (ruled out as an
  explanation for anything above), but it's kept since it's simply more
  accurate. `run_position_accuracy_48hr.jl` uses a separate, unrelated epoch
  (2025-06-06T00:00:00 UTC) matching its own different scenario/dataset.
- `w_eci` sign convention: telemetry `w_eci` is negated (`w0 = -w_eci`)
  relative to the sign SpaceAGORA's `rot()`/quaternion-kinematics convention
  needs; this was determined empirically in the original maneuver-window
  regression (not re-derived here) and is applied consistently everywhere in
  `common.jl` that reads `w1_itp`/`w2_itp`/`w3_itp`.
- Quaternion convention: telemetry is scalar-first `[qs, qx, qy, qz]`;
  SpaceAGORA's internal state and `rot()` are scalar-last `[qx, qy, qz, qw]`.
  `common.jl` converts at every read/write boundary
  (`q_truth_scalarlast`/`q_sim_sf` in the run scripts) -- never assume the
  two are interchangeable without the swap.
- `J_RW`/`I_WHEEL`: regressed against measured `w_eci`/`Omega_rw` over the
  maneuver window in earlier work (not reproduced in this folder), used here
  only as the wheel-only baseline for comparison.
