# N20 encounter and firing recorder validation

## Scope and conclusion

Both requested cases were run with the new integration: 20 helpers at 1000 km,
target at 1050 km, target inclination **0 or 5 deg**, and target initial anomaly
**0 deg**. Duration was copied exactly from the reference:
318567.0415322443 s (50 initial prototype-reference target periods).
These are validation cases, not changes to the planned 100-orbit grid or its
required 1 deg initial target phase.

**The test recorder faithfully recorded these two new simulations.**
It did not change any saved trajectory value or scheduler entry, and an
independent reconstruction agreed with its encounter times within 0.5 microseconds.
The old and new physical runs are not identical, so their firing times need
not match exactly.

The new integration's default output does not automatically include these
interval files. They were produced by the
[test recorder](EncounterTimingRecorder.jl), attached by the
[validation runner](validate_encounter_recording.jl).
No production source or reference output was modified.

## References and setup

The supplied deliverable CSVs were checked against the corresponding saved
SpaceAGORA Feather bundles and matched exactly. Thus the primary comparisons
below are against **SpaceAGORA records**, not the prototype images/videos also
present in those deliverable folders.

- [0 deg deliverable](../../Summary_Table_for_Other_Labs/6_Deliverables/N20_h1000km_t1050km_ih0_it0_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/)
- [5 deg deliverable](../../Summary_Table_for_Other_Labs/6_Deliverables/N20_h1000km_t1050km_ih0_it5_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s/)

All 21 initial Cartesian states were copied from each saved reference.
Helpers retain IDs 1-20 and the target ID is 21; pair IDs require no remapping.
Settings: 227 kg each, 10 kW, B=100, 200 km range, one terminal,
target-helper candidate links, `gve_sma`, central gravity plus J2, no drag,
Tsit5, relative/absolute tolerances 1e-12, and a 10 s maximum step.
The epoch is 2026-01-01 00:00:00, matching the reference runner.
The new runs use the equatorial J2 radius, also used by the comparison project's
current gravity code. The new force uses c=299792458 m/s.

All 210 spacecraft pairs were checked for geometric encounters, although these
cases produced only target-helper encounters. At these altitudes and the 200 km
range, Earth occultation is nonlimiting; no pointing/acquisition model is assumed.
There are 31,858 regular trajectory rows per run, including the final endpoint.

## Encounter and firing totals

All durations below are seconds. Geometric duration is summed across encounters.

| Quantity | 0 deg reference | 0 deg new | 5 deg reference | 5 deg new |
|---|---:|---:|---:|---:|
| Geometric encounters | 11 | 11 | 18 | 18 |
| Total encounter duration | 53872.347596 | 53873.764974 | 8082.239277 | 8083.660773 |
| Recorded firing intervals | 11 | 10 | 10 | 9 |
| Recorded/actual firing duration | 26088.846925 | 26157.356350 | 3850.484340 | 3918.931234 |
| Firing within geometric encounters | 26088.846925 | 26157.356350 | 3850.484340 | 3886.676984 |
| Recorded firing outside range | 0 | 0 | 0 | 32.254250 |

All encounter pairs and chronological occurrences matched; none were unmatched.
The reference contains one initial 0.000379693605 s firing interval in each
case, absent in the new simulation. This accounts for the one-interval count
difference; it is retained as an unmatched row, not silently discarded.

For the 5 deg new run:

```
3918.931234 s actual firing
  = 3886.676984 s within range + 32.254250 s outside range
```

The largest outside-range contribution for one firing interval is 7.735132 s,
within the accepted 10 s step cap. In the coplanar case, no firing extended
beyond a geometric encounter.

### Largest differences per matched encounter

| Absolute difference | 0 deg | 5 deg |
|---|---:|---:|
| Entry time (s) | 6.905823 | 0.646620 |
| Exit time (s) | 4.710235 | 0.626283 |
| Encounter duration (s) | 2.730566 | 1.272903 |
| In-range laser-on duration (s) | 12.244544 | 7.277147 |

The old comparison implementation applies discrete velocity kicks and uses
continuous boundary callbacks. The new implementation applies continuous force
and holds links between accepted steps. This test deliberately does not add
boundary-triggered solver stops or clip the applied firing.
Scheduler update times, force integration, and accumulated trajectory differences
therefore differ. Old/new discrepancies are not estimates of recorder precision,
and cannot all be attributed to a single delayed shutoff.

## How recording was verified

1. **No interference with propagation.** Each case was also run without the
   passive recorder, with otherwise identical configuration and saved fields.
   Every saved column and the complete accepted-step schedule were exactly equal.
2. **Passive geometric roots.** A callback observes completed accepted steps.
   Cubic-Hermite relative-motion bounds and stationary points identify potential
   crossings, including an entry and exit within one step. Range roots and
   boundary states are then evaluated on the solver's current-step dense
   interpolation. No event-triggered timestep is introduced.
3. **Independent saved-data reconstruction.** Encounters were reconstructed
   separately with cubic Hermite interpolation on the regular 10 s position/
   velocity output, rather than accepted steps. Counts and pairs agreed; maximum
   boundary-time differences were 4.592e-7 s (0 deg) and 2.258e-7 s (5 deg),
   passing the explicit 1e-3 s comparison threshold.
4. **Boundary residuals.** Maximum absolute range residuals at recorded roots
   were 4.404e-8 m and 2.246e-7 m, passing the 1e-3 m threshold.
5. **Firing versus actual force.** Each step's solver-integrated laser impulse
   was checked against the previously selected pair and a Simpson estimate of
   BP/(cm) times its changing direction. Maximum vector discrepancies were
   1.081e-14 and 2.869e-12 m/s, below the 1e-9 m/s test threshold. Unselected
   spacecraft/laser-off steps had zero recorded impulse increments.
6. **Persisted interval accounting.** Firing boundaries lie on accepted-step
   times, summed firing duration equals schedule-integrated active time, and
   within-range plus outside-range time equals full firing time. Reopened
   Feather tables reproduce the saved analysis.
7. **Analytical fixtures.** 17 checks cover known crossings, an entire encounter
   between step endpoints (including a 0.5 ms encounter), grazing/no encounter,
   delayed shutoff, helper switching, and clipped start/end intervals.

Both full recorded runs passed 20 checks each; both no-recorder controls passed.
The fresh-process comparisons passed 23 checks each. The 5 deg recorded run
included three rejected solver trials; they did not become recorder samples.

These tests validate the recorder on these cases and fixtures, not absolute
physical accuracy at the microsecond level or completeness for every possible
grazing trajectory. The reference file is not a precision ground truth.
Further grid-wide and step-refinement tests remain appropriate.

## Outputs

Results are under [N20_encounter_recording_validation](results/N20_encounter_recording_validation/):

- `i0/recorded/` and `i5/recorded/`: native trajectory, geometry intervals,
  full actual firing intervals, per-encounter in-range firing, and per-firing
  within/outside-range time totals, in Feather and CSV.
- Each firing record's `start_encounter_id` identifies the encounter open when
  firing began; it is not a containment promise. Zero means no encounter was open.
- `firing_range_audit` is a **per-firing duration summary** in this test, not the
  planned production file with a row for every within/outside-range segment.
- `i0/baseline/` and `i5/baseline/`: no-recorder control trajectories/schedules.
- Each inclination directory also has `encounter_comparison.csv`,
  `laser_interval_comparison.csv`, `independent_10s_reconstruction.csv`,
  and `comparison.toml`.
- Root-level logs preserve test outcomes and peak-memory measurements.

Encounter and firing comparisons use unordered pair plus chronological
occurrence, not run-local IDs. Missing occurrences have blank CSV differences.
No speed/range-rate extrema, optical-work analysis, or RTN production diagnostics
were added by this timing test.

## RAM and reproduction

Each simulation ran sequentially in its **own Julia process**. The recorded
cases peaked at approximately **2.06 GiB RSS each**; their RAM did not add
together. The processes retained only solver endpoints, not full solver histories.
This is not a memory benchmark for 301 spacecraft over 100 orbits.

In a long-lived Julia process, keeping solutions, tables, callback histories,
or caches reachable can accumulate memory. Garbage collection cannot free
objects still referenced, and freed memory may remain reserved by the allocator.
For the grid, separate processes provide a straightforward per-run memory boundary.
OS file cache may remain after a process exits but is generally reclaimable.

From `Tests`, with Julia 1.12.1:

```sh
julia --startup-file=no --threads=1 validate_encounter_recording.jl fixtures
julia --startup-file=no --threads=1 validate_encounter_recording.jl run 0
julia --startup-file=no --threads=1 validate_encounter_recording.jl baseline 0
julia --startup-file=no --threads=1 validate_encounter_recording.jl compare 0
julia --startup-file=no --threads=1 validate_encounter_recording.jl run 5
julia --startup-file=no --threads=1 validate_encounter_recording.jl baseline 5
julia --startup-file=no --threads=1 validate_encounter_recording.jl compare 5
```

The runner activates the new integration's environment. Runs refuse to overwrite
existing result subdirectories; comparisons can be repeated directly.
