# SMA Scheduler Comparison: 50 Helpers over 50 Target Orbital Periods

## 1. Main result

**All three simulations completed successfully and chose the same sequence
of physical helpers for the sustained contacts: slots 50, 49, ..., 25.**
Their target SMA gains were close, but not identical:

| Implementation | Final target SMA (m) | Target SMA gain (m) | Final eccentricity |
|---|---:|---:|---:|
| Kuang prototype | 7,429,621.021573 | **+1,484.421573** | 0.00003950356224 |
| 4 / v1 | 7,429,616.714260 | **+1,480.114260** | 0.00003966913611 |
| 9 integration | 7,429,617.030079 | **+1,480.430079** | 0.00003971232555 |

The largest gain difference is **4.307313 m**, or **0.290710% of the mean
gain**. The v1 and 9 gains differ by only **0.315819 m**.

This demonstrates similar behavior for this equal-power, single-target
scenario. It does not demonstrate that the three schedulers are identical
for arbitrary link networks or unequal link strengths.

## 2. Configuration

Only helper count and duration were changed from the
[20-helper, five-orbit comparison](SMA_SCHEDULER_COMPARISON.md).

| Setting | Value |
|---|---|
| Helpers | 50, initially circular and equatorial at 1000 km altitude |
| Target | 1, initially circular and equatorial at 1050 km altitude |
| Helper phasing | Slots 1 through 50 at 0, 7.2, ..., 352.8 degrees |
| Initial target phase | 0 degrees |
| Spacecraft mass | 227 kg each |
| Allowed links | Target-helper only; no helper-helper links |
| Target terminals in 9 | One; hence at most one active target-helper link |
| Link range | 200 km |
| Input power / magnification | 10,000 W / 100 for every candidate |
| v1 force factors | `eta=1`, `beta=1`, matching the other force magnitudes |
| Gravity | Central point mass; no J2 or drag |
| Earth reference radius | 6,378,136.6 m |
| Propagation gravitational parameter | 3.98600436233e14 m^3/s^2 |
| Initial target SMA | 7,428,136.6 m |
| Initial target orbital period | 6,371.340360499 s |
| Duration | **318,567.018024937 simulation seconds**, approximately 88.49 hours |

"50 orbits" means 50 times the **initial target orbital period**, not an
event-based count of 50 completed revolutions after the orbit changes.
All runs ended at exactly the requested duration and returned `Success`.

For 9, the forbidden-pair set is empty and spacecraft are active. The
prototype additionally enables line-of-sight checking with a 100 km
atmospheric boundary. These eligibility rules were retained, not replaced
by one common scheduler.

## 3. Helper selections and contact history

Raw spacecraft indices differ between implementations:

- Prototype: helpers are spacecraft 1 through 50; target is spacecraft 51.
- v1 and 9: target is spacecraft 1; helpers are spacecraft 2 through 51.
- Physical helper slot $j$ means the helper initially phased at
  $7.2(j-1)$ degrees, regardless of its raw spacecraft index.

After mapping indices into physical slots, all three produced:

```text
50 -> 49 -> 48 -> 47 -> 46 -> 45 -> 44 -> 43 -> 42 -> 41 -> 40 -> 39
   -> 38 -> 37 -> 36 -> 35 -> 34 -> 33 -> 32 -> 31 -> 30 -> 29 -> 28
   -> 27 -> 26 -> 25
```

There are inactive gaps between these contacts; the arrows do not represent
immediate switches between active helpers.

Each run contains **26 non-startup contact segments**: 25 ended before the
simulation finished, and the 26th was still active at the final time.
All end with physical helper slot **25**, initially at **172.8 degrees**.
Thus, the raw final helper index of 25 in the prototype and 26 in v1/9
represents the **same spacecraft**, not a disagreement.

### Selected contact timings

All times below are **simulation seconds**, not wall-clock runtime.
Prototype times are sampled observations; see the measurement caveats.

| Contact | Prototype | v1 | 9 |
|---|---:|---:|---:|
| Slot 50: first observed activation | 9,909.954 | 9,918.455 | 9,911.251 |
| Slot 50: first observed deactivation | 12,512.953 | 12,518.455 | 12,521.251 |
| Slot 38: observed activation | 158,950.878 | 158,958.455 | 158,961.251 |
| Slot 38: observed deactivation | 161,508.697 | 161,518.455 | 161,511.251 |
| Slot 25: final observed activation | 318,111.028 | 318,128.455 | 318,131.251 |
| Simulation cutoff while slot 25 remains active | 318,567.018 | 318,567.018 | 318,567.018 |

The final cutoff is **not a deactivation event**.
The [complete contact table](results/50_helpers_50_orbits/contacts.csv)
records every observed active segment and flags contacts truncated by the
simulation end.

### v1 startup selection

v1 also selected physical helper slot 1 at initialization and deselected it
at **0.0000124087073 s**, about **12.4 microseconds of simulation time**.
Its activation counter therefore reports 27 rather than 26.

That initially radial force direction has zero ideal instantaneous SMA rate.
The tiny startup selection is consistent with numerical roundoff near the
strict positive-score boundary; it is not an additional sustained contact.
The analysis labels a segment as a startup transient only when it begins
at zero and ends before one second. No such segment was recorded in the
other two runs.

## 4. Quantitative schedule comparison

| Pair | Absolute SMA-gain difference (m) | Largest matched-contact start difference (s) | Largest matched-contact end difference (s) | Time-weighted recorded-schedule agreement |
|---|---:|---:|---:|---:|
| Prototype vs v1 | 4.307313 | 18.176 | 16.184 | 99.8566% |
| Prototype vs 9 | 3.991494 | 20.223 | 18.980 | 99.8526% |
| v1 vs 9 | 0.315819 | 7.204 | 7.204 | 99.9373% |

The sequence comparison excludes the v1 startup transient. The time-weighted
comparison includes the entire run, including startup and inactive periods.
It holds each recorded helper selection until the next timestamp and
integrates agreement over the union of both timestamp sets.
Consequently, it is an agreement statistic for the **recorded schedules**,
not an exact reconstruction of the prototype's force history.

| Diagnostic | Prototype | v1 | 9 |
|---|---:|---:|---:|
| Schedule samples / updates | 3,757 saved-state rescorings | 31,866 accepted-step updates | 31,866 accepted-step updates |
| Samples / updates with an active link | 1,092 | 6,466 | 6,470 |
| Estimated active time (s) | 64,654.987 | 64,638.563 | 64,685.767 |
| Estimated duty fraction | 20.2956% | 20.2904% | 20.3052% |
| Median schedule sample interval (s) | 119.597 | 10.000 | 10.000 |
| Maximum schedule sample interval (s) | 124.512 | 10.000 | 10.000 |
| Measured simulation-call runtime (wall-clock s) | 17.113 | 3.291 | 11.683 |

Active time is calculated by summing the recorded active intervals, with the
last interval truncated at the simulation end. **Do not compare raw active
sample counts as duty cycles**, because the sample spacing differs.

Runtime measurements include work and any compilation inside each timed
simulation call, but exclude earlier package loading and later CSV analysis.
The three processes were launched concurrently with one Julia thread and
one OpenBLAS thread each. These values are execution observations, **not a
controlled scheduler-performance benchmark**.

## 5. Why the behavior is similar but not identical

### Same helper choices in this configuration

All candidates have the same force magnitude and the same target mass.
Therefore, at an identical state with identical eligibility, multiplying
the directional GVE score by force/mass does not change its ranking.

The prototype and v1 choose one helper for one target. 9 selects a
maximum-total-score terminal matching, but here all allowed edges share the
same target terminal, so it also selects at most one helper.
Its ability to schedule multiple independent links is not exercised.

Although the helper count increased, their initial adjacent spacing is still
about 927 km, substantially larger than the 200 km link radius. This remains
a sparse-contact example, not a strong test of choosing among many
simultaneously eligible helpers. Matching and unequal-power behavior should
not be generalized from these results.

### Different numerical and scheduling paths

- **Prototype:** its original `run_open_cavity_multi` uses Vern9 with
  `reltol=1e-12`, `abstol=1e-12`, and no explicit 10 s step cap. It selects
  links during RHS force evaluation. The saved schedule is reconstructed
  by rescoring saved states, not logged at every RHS evaluation.
- **v1:** the actual v1 scheduler/force source is loaded through the existing
  wrapper, with a central-gravity ODE. It uses Tsit5, `reltol=1e-9`,
  `abstol=1e-11`, and a 10 s maximum step, updating selection at accepted
  steps.
- **9:** runs through native `run_simulation`, using Tsit5 and the same
  requested orbit tolerances and step cap as the v1 runner. It records
  accepted-step link history.
- **Scoring constant:** v1's GVE source uses 3.986004418e14 m^3/s^2,
  slightly different from the common propagation constant.

The repeated differences in contact timing and evolving orbital states are
consistent with the final SMA differences. This comparison does not isolate
each cause or include a timestep-convergence study, so the 4.307 m spread
should not be attributed solely to the score formula. A `Success` return
also does not by itself establish a numerical error bound.

For the mathematical relationship between the scoring formulas, see
[SMA Scheduler Scoring](SMA_SCHEDULER_SCORING.md).

## 6. Reproduce and inspect

From the workspace root:

```sh
GKSwstype=100 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project="4_SpaceAGORA.jl-main_v1" SciTech_Paper/run_1_prototype.jl 50 50
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project="4_SpaceAGORA.jl-main_v1" SciTech_Paper/run_4_v1.jl 50 50
GKSwstype=100 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project="9_SpaceAGORA.jl-main_full integration" SciTech_Paper/run_9_full.jl 50 50
julia --project="4_SpaceAGORA.jl-main_v1" SciTech_Paper/compare_sma_runs.jl 50 50
```

Executed using Julia **1.12.1**. The tracked repository revision at execution
was `324e7f5a47e043b40203bd9d6b94d1bf3f951859`; the comparison runners and
report are additional workspace files.

The [analysis script](compare_sma_runs.jl) verified successful return codes,
requested helper/orbit counts, common initial SMA and duration, finite
summary metrics, final-minus-initial SMA consistency, CSV schema, strictly
increasing timestamps covering the full duration, valid helper indices,
sample counts, and recorded final selections. It then normalized helper
indices, extracted contacts, and calculated the statistics above.

| Code | Runner | Summary | Recorded schedule |
|---|---|---|---|
| Prototype | [run_1_prototype.jl](run_1_prototype.jl) | [summary.toml](results/50_helpers_50_orbits/prototype/summary.toml) | [schedule.csv](results/50_helpers_50_orbits/prototype/schedule.csv) |
| v1 | [run_4_v1.jl](run_4_v1.jl) | [summary.toml](results/50_helpers_50_orbits/4_v1/summary.toml) | [schedule.csv](results/50_helpers_50_orbits/4_v1/schedule.csv) |
| 9 | [run_9_full.jl](run_9_full.jl) | [summary.toml](results/50_helpers_50_orbits/9_full/summary.toml) | [schedule.csv](results/50_helpers_50_orbits/9_full/schedule.csv) |

Derived outputs: [comparison.toml](results/50_helpers_50_orbits/comparison.toml)
and [contacts.csv](results/50_helpers_50_orbits/contacts.csv).
The earlier five-orbit result files were not overwritten.
