# Study preflight checks

Checks completed on 2026-10-18 using the new integration. The production
simulation library was not modified. Test code, results, and reproducible
commands are retained here.

## Checklist outcome

| Requested check | Outcome and scope |
|---|---|
| 1. Laser on/off and zero power | Passed for both SMA and inclination schedulers in short J2 runs |
| 2. BP/c force and own-mass response | Passed, including unequal endpoint masses |
| 3. Initial geometry and LOS | All 105 planned initial geometries are noncoincident; overlap rejection tested; analytic LOS bound and runtime clearance checks passed |
| 4. Event edge cases | Passed short, multiple, simultaneous, grazing, helper-switch, and clipped start/end fixtures |
| 5. Numerical accuracy | Refinement experiments completed; sensitivity remains significant for some quantities. Production RTN solver-state tests are blocked because those states do not yet exist |
| 6. Recording and range segments | Passed in test harness, including fresh-process reads, complete segment partitions, and recorder noninterference |
| 7. Mechanical balances | Passed force/torque, work-energy, J2 axial angular momentum, and no-J2 full angular momentum tests |

Passed checks can be removed from the plan's outstanding checklist. Their
implementation still needs to be wired into the eventual production grid runner;
these results are not a claim that the whole production pipeline exists.

## Test matrix

Two selected 100-orbit cases:

| Case | Helpers | Target | Scheduler |
|---|---|---|---|
| `sma850` | 20 at 1000 km, inclination 0 deg | 850 km, inclination 1 deg, nu=1 deg | `gve_sma` |
| `inc1050` | 20 at 1000 km, inclination 0 deg | 1050 km, inclination 1 deg, nu=1 deg | `gve_inc` |

Each case used:

- Step caps 10, 5, and 2.5 s with 10 s saving.
- Saving intervals 10, 5, and 1 s with the 10 s step cap.
- A separate 10 s control with event recording disabled.

This is 12 runs, including the two controls. A thirteenth 100-orbit
`inc1050` run disabled J2 to test full-vector angular momentum conservation.
An initial two-orbit smoke run is also retained.

Settings: central gravity plus J2 except the explicitly labelled no-J2 run,
227 kg each, P=10 kW, B=100, range=200 km, one terminal, target-helper
candidate links, Tsit5, relative/absolute tolerances 1e-12. All 210 pairs were
monitored geometrically. Every simulation ran sequentially in a separate process.
All force is propagated by the existing continuous-force engine.

## Physics, force, and geometry

- Measured force: **0.003335640951981521 N**, consistent with BP/c, not 2BP/c.
- Pair forces are equal and opposite; pair torque residual was zero in the
  tested fixture (required below 1e-8 N m).
- With one endpoint twice as massive, its acceleration is half as large.
- Short 120 s laser-on/off J2 runs produced target position differences of
  **0.105673 m** (`gve_sma`) and **0.105673 m** (`gve_inc`).
- Zero-power runs had exactly zero integrated laser impulse and agreed with
  laser-model-disabled gravity trajectories within 1e-5 in Cartesian components.
- The equatorial inclination-scheduler fixture produced no laser impulse.
- All combinations of 5 target altitudes, 3 inclinations, and 7 helper counts
  were checked at the required 1 deg target phase. The minimum initial pair
  separation was **25,754.542 m**.
- A same-altitude, zero-phase overlapping case was explicitly rejected by the
  test configuration guard. This guard must be retained in the production runner.

### LOS

For any point on a line segment no longer than 200 km, the reverse triangle
inequality gives distance from Earth's center at least the endpoint radius
minus 200 km. At the lowest prescribed initial altitude, 850 km, that gives
at least **650 km above the equatorial-radius sphere**. Thus Earth cannot
occlude a 200 km link in that altitude domain.

The full runs additionally checked actual spacecraft radii at every accepted
step against `Earth_equatorial_radius + 200 km`; they remained well above it.
This validates the stated operating domain and selected trajectories, not
arbitrary future orbits approaching Earth. Pointing/acquisition is outside scope.

## Events and range-audit segments

The test recorder was extended with a segment-level auxiliary output:

- One row per contiguous within-range or outside-range part of actual firing.
- Run ID, firing ID, pair IDs, nullable encounter ID, start/end time, duration,
  and within-range flag.
- Full coverage of every firing interval, with no gaps or overlaps.
- Duration sums match actual firing, per-encounter overlaps, and outside-range
  totals. Outside-range portions are not silently removed.

Fixtures include a 0.5 ms encounter wholly between step endpoints, two
simultaneous target contacts, three separate encounters within one step,
grazing contact, helper switching, open initial/final encounters, wholly
outside-range firing, and no firing.

The new empty-case tests found a test-recorder bug: an empty encounter analysis
could produce an untyped laser-duration column. It now explicitly uses
Float64, and empty results survive summation and Feather round-trips.
Invalid interval times, duplicate IDs, and overlapping intervals for the same
pair are rejected explicitly.

## Step-cap sensitivity

These are **measured changes**, not certified error bounds. The changes include
the effect of scheduler updates occurring at different accepted-step times.
Osculating element changes include J2 and are not laser-on minus laser-off.

| Quantity | SMA case 10 s | 5 s | 2.5 s |
|---|---:|---:|---:|
| Total SMA change (m) | 205.971288 | 204.967854 | 203.853379 |
| Target R integral estimate (m/s) | -0.292637252 | -0.293386280 | -0.293283013 |
| Target T integral estimate (m/s) | 0.103791926 | 0.103275321 | 0.102700967 |
| Target N integral estimate (m/s) | 0.020407574 | 0.020285721 | 0.020466945 |
| Encounters | 61 | 61 | 61 |
| Actual firing (s) | 23965.218255 | 24010.811056 | 23993.607802 |

| Quantity | Inclination case 10 s | 5 s | 2.5 s |
|---|---:|---:|---:|
| Total inclination change (deg) | 0.001078629775 | 0.001078554020 | 0.001077951498 |
| Total SMA change (m) | 19.749804 | 21.570383 | 22.358195 |
| Target R integral estimate (m/s) | 0.284702891 | 0.284960805 | 0.285123711 |
| Target T integral estimate (m/s) | 0.008557323 | 0.009457570 | 0.009846869 |
| Target N integral estimate (m/s) | -0.006005305 | -0.005575745 | -0.005631290 |
| Encounters | 21 | 21 | 21 |
| Actual firing (s) | 48220.887480 | 48258.476325 | 48258.230491 |

All encounter pairs and chronological occurrences matched across refinements.

| Largest absolute change | SMA: 10 to 5 s | SMA: 5 to 2.5 s | Inclination: 10 to 5 s | Inclination: 5 to 2.5 s |
|---|---:|---:|---:|---:|
| Encounter entry (s) | 5.466390 | 5.454848 | 21.281818 | 2.421256 |
| Encounter exit (s) | 3.454400 | 4.074010 | 18.962991 | 3.147018 |
| Encounter duration (s) | 4.178042 | 4.159810 | 17.694763 | 0.725762 |
| Final target position (m) | 590.373616 | 606.001206 | 948.349849 | 203.752124 |

The SMA gain changes by 2.117909 m from 10 to 2.5 s, approximately 1.03%
of its 10 s value. The inclination gain changes by 6.782772e-7 deg,
approximately 0.063%. Some small net RTN components vary more in relative terms
and do not converge monotonically. Use the absolute values above rather than
calling all maneuver quantities equally accurate.

No study-wide error budget has been specified. Therefore the numerical-accuracy
item remains open; the accepted 10 s scheduler is not silently replaced.
Long-term encounter differences can exceed one scheduler step because trajectories
and later schedules diverge cumulatively.

## Saving-interval sensitivity and recorder interference

For both cases, saving at 10, 5, or 1 s produced:

- Exactly identical accepted-step schedules and geometric/firing records.
- Exactly identical native trajectory and diagnostic values at shared sample times.
- Exactly identical final orbital changes, RTN estimates, and firing totals.

Sampled kinematic extrema can still change. With 10-to-1 s saving:

- SMA-case minimum relative speed changed by **-0.003498 m/s**.
- Inclination-case minimum relative speed changed by **-0.001620 m/s**.
- Inclination-case minimum/maximum range rate changed by approximately
  **-0.001214 / +0.002159 m/s**.

These extrema include recorded encounter boundaries as well as regular samples,
and remain sampled estimates.

Removing the event recorder from otherwise identical diagnostic runs produced
exactly identical schedules, native saved states, and diagnostics in both cases.
The earlier N20 tests also compared against runs without the passive diagnostic
callback; see [that report](ENCOUNTER_RECORDING_VALIDATION_N20.md).

An independent reconstruction from regular 10 s position/velocity samples
agreed with current-step dense-interpolant roots within **1.189259e-6 s**
for the SMA case and **1.983717e-7 s** for the inclination case, below the
explicit 1e-3 s comparison threshold. Counts and pairs also matched.
This establishes numerical recording consistency for these trajectories,
not microsecond physical accuracy or subsecond scheduler switching.

## Work, energy, and angular momentum

The test diagnostics integrate actual selected-link forces using composite
Simpson quadrature on solver dense interpolation, including out-of-range firing.
Work is integrated separately for every spacecraft as integral(F dot v) dt.
RTN uses the target's changing frame; gravity is excluded.

**These are not the planned RTN solver states.** Their quadrature can be tested,
but the nonexistent production solver-state implementation cannot be certified.
Along-track sign separation is quadrature-based, not root-exact at every sign change.

Energy includes central gravity and, when enabled, the J2 potential:

```
U/m = -mu/r + mu*J2*R_equatorial^2/(2*r^3) * (3*(z/r)^2 - 1)
```

The tests check system mechanical energy change minus summed optical work
at every accepted step. The explicit acceptance threshold is
`20 J + 1e-4 * abs(final total optical work)`.
Angular-momentum tolerance is 1e-9 relative to the initial full-vector norm.

| Mechanical check | Maximum discrepancy |
|---|---:|
| J2 SMA case, 10 s: energy minus work | 2.067577 J |
| J2 inclination case, 10 s: energy minus work | 2.158363 J |
| J2 SMA case, 2.5 s: energy minus work | 0.008345 J |
| J2 inclination case, 2.5 s: energy minus work | 0.005713 J |
| No-J2 100-orbit case: energy minus work | 2.174376 J |
| J2 axial angular momentum, largest relative discrepancy | 8.363161e-12 |
| No-J2 full-vector angular momentum, relative discrepancy | 8.445713e-12 |

All passed. Mass-weighted integrated laser impulses also balance below
1e-8 kg m/s. Total satellite linear momentum is not used as a conserved
quantity because Earth supplies external gravitational force.

## Evidence and reproduction

- [Preflight runner](preflight_study_checks.jl)
- [Fresh-process verifier](verify_study_preflight.jl)
- [Shared test recorder](EncounterTimingRecorder.jl)
- [All results and logs](results/study_preflight/)
- [Sensitivity results](results/study_preflight/sensitivity_runs.csv)
- [Pairwise refinement differences](results/study_preflight/sensitivity_differences.csv)

The fixtures passed **298 assertions**, all 13 full runs passed their checks,
and fresh-process persisted-data verification passed **118,190 assertions**.
The large verification count includes per-sample range-label checks; it is not
118,190 independent physical scenarios.

From `Tests`:

```sh
julia --startup-file=no --threads=1 preflight_study_checks.jl fixtures
for case in sma850 inc1050; do
  for dt in 10 5 2.5; do
    julia --startup-file=no --threads=1 preflight_study_checks.jl run "$case" "$dt" 10
  done
  for interval in 5 1; do
    julia --startup-file=no --threads=1 preflight_study_checks.jl run "$case" 10 "$interval"
  done
  julia --startup-file=no --threads=1 preflight_study_checks.jl run "$case" 10 10 control
done
julia --startup-file=no --threads=1 preflight_study_checks.jl run inc1050 10 10 noJ2
julia --startup-file=no --threads=1 verify_study_preflight.jl
```

Simulation commands refuse to overwrite existing result directories. The verifier
can be repeated directly. The separate production-wide stale-include/routing
test-suite issues were outside this checklist execution and remain unresolved.
